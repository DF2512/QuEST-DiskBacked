#include "diskbackedstate.h"
#include "gatescheduler.h"
#include "run.h"
#include "gates.h"
#include "chunkmanager.h"
#include <vector>
#include <iostream>
#include <omp.h>
#include <atomic>
#include <condition_variable>
#include <chrono>
#include <thread>
#include <mutex>
#include <optional>

// Shared memory pool counter
std::mutex memMtx;
std::condition_variable memCv;
int inFlightBlocks = 0;
constexpr int maxBlocksInMemory = 3;

// ─────────────────────────────────────────────────────────────
// Apply all gates from a subcircuit to a single block buffer
// ─────────────────────────────────────────────────────────────
void applySubCircuitToBlock(const SubCircuit& sub, std::vector<qcomp>& buffer, int qubitsPerBlock) {
    for (const GateOp& gate : sub.gates) {
        int target = sub.permutation[gate.target];
        int control = (gate.control != -1) ? sub.permutation[gate.control] : -1;

        if (gate.type == GateType::Swap) {
            if (gate.bufferFn) {
                gate.bufferFn(buffer, qubitsPerBlock);
            }
            continue;
        }

        qindex stride = 1ULL << target;
        qindex blockSize = 2 * stride;
        if (buffer.size() < blockSize) continue;

        qindex numGroups = buffer.size() / blockSize;

        #pragma omp parallel for schedule(static)
        for (qindex g = 0; g < numGroups; ++g) {
            qindex base = g * blockSize;
            for (qindex i = 0; i < stride; ++i) {
                qindex idx0 = base + i;
                qindex idx1 = idx0 + stride;

                if (gate.type == GateType::Single) {
                    gate.fn(buffer[idx0], buffer[idx1]);
                } else if (gate.type == GateType::Controlled) {
                    if (((idx0 >> control) & 1) == 1)
                        gate.fn(buffer[idx0], buffer[idx1]);
                }
            }
        }
    }
}


void applySubCircuitToBlockDebug(const SubCircuit& sub, std::vector<qcomp>& buffer, int qubitsPerBlock) {
    std::cout << "[Debug] Entering applySubCircuitToBlockDebug, buffer size: " << buffer.size() << ", qubitsPerBlock: " << qubitsPerBlock << "\n";
    
    std::cout << "[Debug] SubCircuit has " << sub.gates.size() << " gates\n";
    
    if (sub.gates.empty()) {
        std::cout << "[Debug] No gates to apply\n";
        return;
    }
    
    std::cout << "[Debug] About to start gate loop\n";
    
    for (size_t gateIdx = 0; gateIdx < sub.gates.size(); ++gateIdx) {
        std::cout << "[Debug] Processing gate " << gateIdx << "\n";
        
        const GateOp& gate = sub.gates[gateIdx];
        std::cout << "[Debug] Gate " << gateIdx << " type: " << static_cast<int>(gate.type) << "\n";

        int logicalTarget = gate.target;
        int physicalTarget = sub.permutation[gate.target];
        int logicalControl = gate.control;
        int physicalControl = (gate.control != -1) ? sub.permutation[gate.control] : -1;

        if (gate.type == GateType::Swap) {
            std::cout << "[Debug] Gate #" << gateIdx << ": type=" << static_cast<int>(gate.type)
                      << ", target=" << logicalTarget << ", control=" << logicalControl << "\n";
            std::cout << "Gate #" << gateIdx << " | SWAP(" << logicalTarget << ", " << logicalControl << ")\n";
            if (gate.bufferFn) {
                gate.bufferFn(buffer, qubitsPerBlock);
            }
            continue;
        }

        qindex stride = 1ULL << physicalTarget;
        qindex blockSize = 2 * stride;
        if (buffer.size() < blockSize) continue;

        qindex numGroups = buffer.size() / blockSize;

        std::cout << "[Debug] Gate #" << gateIdx << ": type=" << static_cast<int>(gate.type)
                  << ", target=" << logicalTarget << ", control=" << logicalControl << "\n";

        if (logicalControl != -1) {
            std::cout << "[Debug] Gate #" << gateIdx << ": type=" << static_cast<int>(gate.type)
                      << ", target=" << logicalTarget << ", control=" << logicalControl << "\n";
        }
        std::cout << std::endl;

        qindex debug_idx0 = 0;
        qindex debug_idx1 = stride;

        qcomp before0 = buffer[debug_idx0];
        qcomp before1 = buffer[debug_idx1];

        std::cout << "  Before gate: index " << debug_idx0 << " = " << before0
                  << ", index " << debug_idx1 << " = " << before1 << std::endl;

        #pragma omp parallel for schedule(static)
        for (qindex g = 0; g < numGroups; ++g) {
            qindex base = g * blockSize;
            for (qindex i = 0; i < stride; ++i) {
                qindex idx0 = base + i;
                qindex idx1 = idx0 + stride;

                if (gate.type == GateType::Single) {
                    gate.fn(buffer[idx0], buffer[idx1]);
                } else if (gate.type == GateType::Controlled) {
                    if (((idx0 >> physicalControl) & 1) == 1)
                        gate.fn(buffer[idx0], buffer[idx1]);
                }
            }
        }

        qcomp after0 = buffer[debug_idx0];
        qcomp after1 = buffer[debug_idx1];
        std::cout << "  After gate:  index " << debug_idx0 << " = " << after0
                  << ", index " << debug_idx1 << " = " << after1 << std::endl;
    }
}

std::vector<GateOp> scheduleSwaps(const std::vector<std::pair<int, int>>& swaps) {
    std::cout << "[Debug] scheduleSwaps called with " << swaps.size() << " swaps\n";
    GateScheduler sched;
    for (const auto& [a, b] : swaps) {
        std::cout << "[Debug] Scheduling SWAP(" << a << ", " << b << ")\n";
        addSwap(sched, a, b);
    }
    auto result = sched.getSchedule();
    std::cout << "[Debug] scheduleSwaps returning " << result.size() << " gates\n";
    return result;
}


void runPipeline(GateScheduler& scheduler, DiskBackedState& state, bool verbose) {
    int numQubits = state.getNumQubits();
    int numLocalQubits = state.getNumQubitsPerBlock();
    int numBlocks = state.getNumBlocks();
    int qubitsPerBlock = numLocalQubits;
    PermutationTracker& tracker = state.getPermutationTracker();

    // Partition into subcircuits
    std::vector<SubCircuit> subcircuits = scheduler.partitionIntoSubcircuits(numQubits, numLocalQubits, state, verbose);

    // Generate transitions
    std::vector<Transition> transitions = tracker.generateTransitions(subcircuits);

    std::cout << "[Pipeline] Generated " << transitions.size() << " transitions\n";
    for (size_t i = 0; i < transitions.size(); ++i) {
        std::cout << "[Pipeline] Transition " << i << ":\n";
        std::cout << "  swap1 size: " << transitions[i].swap1.size() << "\n";
        std::cout << "  swap2 size: " << transitions[i].swap2.size() << "\n";
        std::cout << "  swapLevels size: " << transitions[i].swapLevels.size() << "\n";
    }

    for (size_t i = 0; i < subcircuits.size(); ++i) {
        std::cout << "\n[Pipeline] Processing subcircuit " << i << "\n";
        std::cout << "[Pipeline] Current chunkMap: ";
        for (int j = 0; j < tracker.chunkMap.size(); ++j) {
            std::cout << tracker.chunkMap[j] << " ";
        }
        std::cout << "\n";

        std::vector<std::vector<int>> blockChunkMapping = tracker.getBlockChunkMapping();

        // Timing vectors
        std::vector<double> readerTimes(numBlocks, 0.0);
        std::vector<double> processorTimes(numBlocks, 0.0);
        std::vector<double> writerTimes(numBlocks, 0.0);

        // Set up triple buffer queues
        ThreadSafeQueue<BlockData> readQueue(maxBlocksInMemory);
        ThreadSafeQueue<BlockData> processQueue(maxBlocksInMemory);

        // --- Reader Thread ---
        std::thread reader([&]() {
            for (int blockIdx = 0; blockIdx < numBlocks; ++blockIdx) {
                auto t0 = std::chrono::high_resolution_clock::now();
                // Memory management
                std::unique_lock<std::mutex> lock(memMtx);
                memCv.wait(lock, []() { return inFlightBlocks < maxBlocksInMemory; });
                inFlightBlocks++;

                BlockData block;
                block.blockIdx = blockIdx;
                block.chunkIndices = blockChunkMapping[blockIdx];
                state.loadBlock(block.blockIdx, block.chunkIndices, block.buffer);

                // Apply swap2 if not the first subcircuit
                if (i > 0 && !transitions[i-1].swap2.empty()) {
                    std::cout << "[Reader] Applying swap2 to block " << blockIdx << "\n";
                    SubCircuit swap2;
                    swap2.gates = scheduleSwaps(transitions[i-1].swap2);
                    swap2.permutation = transitions[i-1].interimTarget;
                    if (blockIdx == 0)
                        applySubCircuitToBlockDebug(swap2, block.buffer, qubitsPerBlock);
                    else
                        applySubCircuitToBlock(swap2, block.buffer, qubitsPerBlock);
                }
                auto t1 = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> elapsed = t1 - t0;
                readerTimes[blockIdx] = elapsed.count();
                readQueue.push(std::move(block));
            }
            readQueue.setFinished();
        });

        // --- Processor Thread ---
        std::thread processor([&]() {
            BlockData block;
            while (readQueue.pop(block)) {
                auto t0 = std::chrono::high_resolution_clock::now();
                if (block.blockIdx == 0)
                    applySubCircuitToBlockDebug(subcircuits[i], block.buffer, qubitsPerBlock);
                else
                    applySubCircuitToBlock(subcircuits[i], block.buffer, qubitsPerBlock);
                auto t1 = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> elapsed = t1 - t0;
                processorTimes[block.blockIdx] = elapsed.count();
                processQueue.push(std::move(block));
            }
            processQueue.setFinished();
        });

        // --- Writer Thread ---
        std::thread writer([&]() {
            BlockData block;
            while (processQueue.pop(block)) {
                auto t0 = std::chrono::high_resolution_clock::now();
                // Apply swap1 if not the last subcircuit
                if (i < transitions.size() && !transitions[i].swap1.empty()) {
                    SubCircuit swap1;
                    swap1.gates = scheduleSwaps(transitions[i].swap1);
                    swap1.permutation = subcircuits[i].permutation;
                    if (block.blockIdx == 0)
                        applySubCircuitToBlockDebug(swap1, block.buffer, qubitsPerBlock);
                    else
                        applySubCircuitToBlock(swap1, block.buffer, qubitsPerBlock);
                }
                state.saveBlock(block.blockIdx, block.chunkIndices, block.buffer);

                // Memory management
                {
                    std::lock_guard<std::mutex> lock(memMtx);
                    inFlightBlocks--;
                }
                memCv.notify_one();
                auto t1 = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> elapsed = t1 - t0;
                writerTimes[block.blockIdx] = elapsed.count();
            }
        });

        // Join threads
        reader.join();
        processor.join();
        writer.join();

        // Print per-block and average times
        auto printTimes = [](const std::string& label, const std::vector<double>& times) {
            double sum = 0.0;
            std::cout << "[Timing] " << label << " times per block:\n";
            for (size_t i = 0; i < times.size(); ++i) {
                std::cout << "  Block " << i << ": " << times[i] << " s\n";
                sum += times[i];
            }
            std::cout << "  Average: " << (sum / times.size()) << " s\n";
        };
        printTimes("Reader", readerTimes);
        printTimes("Processor", processorTimes);
        printTimes("Writer", writerTimes);

        // --- After all blocks processed, update chunkMap if not last subcircuit ---
        if (i < transitions.size()) {
            std::cout << "[Pipeline] Applying area swap shuffle for transition " << i << "\n";
            for (int level : transitions[i].swapLevels) {
                std::cout << "[Pipeline] Area swap level: " << level << "\n";
                int maxLevel = *std::max_element(transitions[i].swapLevels.begin(), transitions[i].swapLevels.end());
                std::vector<int> newChunkMap = tracker.areaSwapShuffle(tracker.chunkMap, level, maxLevel);
                tracker.chunkMap = newChunkMap;
            }
            std::cout << "[Pipeline] New chunkMap: ";
            for (int j = 0; j < tracker.chunkMap.size(); ++j) {
                std::cout << tracker.chunkMap[j] << " ";
            }
            std::cout << "\n";
        }
    }
}

// Pipeline function without triple buffering, for benchmarking reasons
/*
void runPipeline(GateScheduler& scheduler, DiskBackedState& state, bool verbose) {
    int numQubits = state.getNumQubits();
    int numLocalQubits = state.getNumQubitsPerBlock();
    int numBlocks = state.getNumBlocks();
    int qubitsPerBlock = numLocalQubits;
    PermutationTracker& tracker = state.getPermutationTracker();

    // Partition into subcircuits
    std::vector<SubCircuit> subcircuits = scheduler.partitionIntoSubcircuits(numQubits, numLocalQubits, state, verbose);

    // Generate transitions
    std::vector<Transition> transitions = tracker.generateTransitions(subcircuits);

    std::cout << "[PipelineNoTripleBuffer] Generated " << transitions.size() << " transitions\n";
    for (size_t i = 0; i < transitions.size(); ++i) {
        std::cout << "[PipelineNoTripleBuffer] Transition " << i << ":\n";
        std::cout << "  swap1 size: " << transitions[i].swap1.size() << "\n";
        std::cout << "  swap2 size: " << transitions[i].swap2.size() << "\n";
        std::cout << "  swapLevels size: " << transitions[i].swapLevels.size() << "\n";
    }

    for (size_t i = 0; i < subcircuits.size(); ++i) {
        std::cout << "\n[PipelineNoTripleBuffer] Processing subcircuit " << i << "\n";
        std::cout << "[PipelineNoTripleBuffer] Current chunkMap: ";
        for (int j = 0; j < tracker.chunkMap.size(); ++j) {
            std::cout << tracker.chunkMap[j] << " ";
        }
        std::cout << "\n";

        std::vector<std::vector<int>> blockChunkMapping = tracker.getBlockChunkMapping();

        // Timing vectors
        std::vector<double> readerTimes(numBlocks, 0.0);
        std::vector<double> processorTimes(numBlocks, 0.0);
        std::vector<double> writerTimes(numBlocks, 0.0);

        for (int blockIdx = 0; blockIdx < numBlocks; ++blockIdx) {
            auto t0 = std::chrono::high_resolution_clock::now();
            BlockData block;
            block.blockIdx = blockIdx;
            block.chunkIndices = blockChunkMapping[blockIdx];
            state.loadBlock(block.blockIdx, block.chunkIndices, block.buffer);
            auto t1 = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsedRead = t1 - t0;
            readerTimes[blockIdx] = elapsedRead.count();

            auto t2 = std::chrono::high_resolution_clock::now();
            // Apply swap2 if not the first subcircuit
            if (i > 0 && !transitions[i-1].swap2.empty()) {
                std::cout << "[NoTripleBuffer] Applying swap2 to block " << blockIdx << "\n";
                SubCircuit swap2;
                swap2.gates = scheduleSwaps(transitions[i-1].swap2);
                swap2.permutation = transitions[i-1].interimTarget;
                if (blockIdx == 0)
                    applySubCircuitToBlockDebug(swap2, block.buffer, qubitsPerBlock);
                else
                    applySubCircuitToBlock(swap2, block.buffer, qubitsPerBlock);
            }

            if (block.blockIdx == 0)
                applySubCircuitToBlockDebug(subcircuits[i], block.buffer, qubitsPerBlock);
            else
                applySubCircuitToBlock(subcircuits[i], block.buffer, qubitsPerBlock);
            auto t3 = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsedProc = t3 - t2;
            processorTimes[blockIdx] = elapsedProc.count();

            // Apply swap1 if not the last subcircuit
            auto t4 = std::chrono::high_resolution_clock::now();
            if (i < transitions.size() && !transitions[i].swap1.empty()) {
                SubCircuit swap1;
                swap1.gates = scheduleSwaps(transitions[i].swap1);
                swap1.permutation = subcircuits[i].permutation;
                if (block.blockIdx == 0)
                    applySubCircuitToBlockDebug(swap1, block.buffer, qubitsPerBlock);
                else
                    applySubCircuitToBlock(swap1, block.buffer, qubitsPerBlock);
            }
            state.saveBlock(block.blockIdx, block.chunkIndices, block.buffer);
            auto t5 = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsedWrite = t5 - t4;
            writerTimes[blockIdx] = elapsedWrite.count();
        }

        // Print per-block and average times
        auto printTimes = [](const std::string& label, const std::vector<double>& times) {
            double sum = 0.0;
            std::cout << "[Timing] " << label << " times per block:\n";
            for (size_t i = 0; i < times.size(); ++i) {
                std::cout << "  Block " << i << ": " << times[i] << " s\n";
                sum += times[i];
            }
            std::cout << "  Average: " << (sum / times.size()) << " s\n";
        };
        printTimes("Reader", readerTimes);
        printTimes("Processor", processorTimes);
        printTimes("Writer", writerTimes);

        // --- After all blocks processed, update chunkMap if not last subcircuit ---
        if (i < transitions.size()) {
            std::cout << "[PipelineNoTripleBuffer] Applying area swap shuffle for transition " << i << "\n";
            for (int level : transitions[i].swapLevels) {
                std::cout << "[PipelineNoTripleBuffer] Area swap level: " << level << "\n";
                int maxLevel = *std::max_element(transitions[i].swapLevels.begin(), transitions[i].swapLevels.end());
                std::vector<int> newChunkMap = tracker.areaSwapShuffle(tracker.chunkMap, level, maxLevel);
                tracker.chunkMap = newChunkMap;
            }
            std::cout << "[PipelineNoTripleBuffer] New chunkMap: ";
            for (int j = 0; j < tracker.chunkMap.size(); ++j) {
                std::cout << tracker.chunkMap[j] << " ";
            }
            std::cout << "\n";
        }
    }
}
*/