#include "diskbackedstate.h"
#include "gatescheduler.h"
#include "run.h"
#include "chunkmanager.h"
#include "qureg.h"
#include "hook.h"
#include <vector>
#include <iostream>
#include <omp.h>
#include <atomic>
#include <condition_variable>
#include <chrono>
#include <thread>
#include <mutex>
#include <optional>
#include <numeric>
#include <algorithm>

// Shared memory pool counter
std::mutex memMtx;
std::condition_variable memCv;
int inFlightBlocks = 0;


// Helper: Given a logical qubit index, find its physical position in the permutation
static int logicalToPhysical(const std::vector<int>& permutation, int logicalQubit) {
    auto it = std::find(permutation.begin(), permutation.end(), logicalQubit);
    if (it == permutation.end()) {
        throw std::runtime_error("logicalToPhysical: Logical qubit not found in permutation");
    }
    return static_cast<int>(std::distance(permutation.begin(), it));
}

// ─────────────────────────────────────────────────────────────
// Apply all gates from a subcircuit to a single block buffer
// ─────────────────────────────────────────────────────────────
void applySubCircuitToBlock(const SubCircuit& sub, std::vector<qcomp>& buffer, int qubitsPerBlock) {
    Qureg tempQureg = createTempQureg(buffer, qubitsPerBlock);
    for (const GateOp& gate : sub.gates) {
        GateOp mappedGate = gate;
        mappedGate.target = logicalToPhysical(sub.permutation, gate.target);
        if (gate.control != -1)
            mappedGate.control = logicalToPhysical(sub.permutation, gate.control);
        GateScheduler::applyGateOpToQureg(mappedGate, tempQureg);
    }
}

void applySubCircuitToBlockDebug(const SubCircuit& sub, std::vector<qcomp>& buffer, int qubitsPerBlock) {
    std::cout << "[Debug] Entering applySubCircuitToBlockDebug, buffer size: " << buffer.size() << ", qubitsPerBlock: " << qubitsPerBlock << "\n";
    std::cout << "[Debug] SubCircuit has " << sub.gates.size() << " gates\n";
    if (sub.gates.empty()) {
        std::cout << "[Debug] No gates to apply\n";
        return;
    }
    std::cout << "[Debug] Subcircuit permutation (physical->logical): ";
    for (size_t i = 0; i < sub.permutation.size(); ++i) {
        std::cout << i << "->" << sub.permutation[i] << " ";
    }
    std::cout << "\n";
    // Print the block's logical qubits (physical pos -> logical qubit)
    std::vector<int> blockLogicalQubits(qubitsPerBlock, -1);
    for (int phys = 0; phys < qubitsPerBlock; ++phys) {
        blockLogicalQubits[phys] = sub.permutation[phys];
    }
    std::cout << "[Debug] Block physical positions (0.." << (qubitsPerBlock-1) << ") contain logical qubits: ";
    for (int i = 0; i < qubitsPerBlock; ++i) {
        std::cout << i << "(logical " << blockLogicalQubits[i] << ") ";
    }
    std::cout << "\n";
    std::cout << "[Debug] About to start gate loop\n";
    Qureg tempQureg = createTempQureg(buffer, qubitsPerBlock);
    for (size_t gateIdx = 0; gateIdx < sub.gates.size(); ++gateIdx) {
        const GateOp& gate = sub.gates[gateIdx];
        GateOp mappedGate = gate;
        mappedGate.target = logicalToPhysical(sub.permutation, gate.target);
        if (gate.control != -1)
            mappedGate.control = logicalToPhysical(sub.permutation, gate.control);
        // Check if mapped target/control are within block
        bool targetInBlock = mappedGate.target >= 0 && mappedGate.target < qubitsPerBlock;
        bool controlInBlock = (gate.control == -1) || (mappedGate.control >= 0 && mappedGate.control < qubitsPerBlock);
        if (!targetInBlock || !controlInBlock) {
            std::cout << "[ERROR] Gate refers to a qubit not present in this block!\n";
            std::cout << "  Gate: type " << static_cast<int>(gate.type)
                      << ", logical target: " << gate.target << ", mapped target: " << mappedGate.target
                      << ", logical control: " << gate.control << ", mapped control: " << mappedGate.control << "\n";
            std::cout << "  Block logical qubits: ";
            for (int i = 0; i < qubitsPerBlock; ++i) {
                std::cout << i << "(logical " << blockLogicalQubits[i] << ") ";
            }
            std::cout << "\n  Subcircuit permutation: ";
            for (size_t i = 0; i < sub.permutation.size(); ++i) {
                std::cout << i << "->" << sub.permutation[i] << " ";
            }
            std::cout << "\n";
        }
        GateScheduler::applyGateOpToQureg(mappedGate, tempQureg);
    }
}

std::vector<GateOp> scheduleSwaps(const std::vector<std::pair<int, int>>& swaps) {
    //std::cout << "[Debug] scheduleSwaps called with " << swaps.size() << " swaps\n";
    GateScheduler sched;
    for (const auto& [a, b] : swaps) {
        //std::cout << "[Debug] Scheduling SWAP(" << a << ", " << b << ")\n";
        sched.addSwap(a, b);
    }
    auto result = sched.getSchedule();
    //std::cout << "[Debug] scheduleSwaps returning " << result.size() << " gates\n";
    // Print the swaps for clarity
    //for (const auto& op : result) {
        //std::cout << "[Debug] SWAP gate: target " << op.target << ", control " << op.control << "\n";
    //}
    return result;
}


void runCircuit(GateScheduler& scheduler, DiskBackedState& state, bool verbose) {
    int numQubits = state.getNumQubits();
    int numLocalQubits = state.getNumQubitsPerBlock();
    int numBlocks = state.getNumBlocks();
    int qubitsPerBlock = numLocalQubits;
    int maxBlocksInMemory = state.getMaxBlocksInMemory();
    PermutationTracker& tracker = state.getPermutationTracker();

    std::vector<SubCircuit> subcircuits = scheduler.partitionIntoSubcircuits(numQubits, numLocalQubits, state, verbose);
    std::vector<std::vector<int>> blockChunkMapping = tracker.getBlockChunkMapping();
    std::vector<Transition> transitions = tracker.generateTransitions(subcircuits);

    if (tracker.currentPermutation != subcircuits[0].permutation) {
        std::cerr << "[Error] Initial permutation does not match first subcircuit's permutation!\n";
        return;
    }

    for (size_t i = 0; i < subcircuits.size(); ++i) {
        tracker.currentPermutation = subcircuits[i].permutation;
        std::vector<std::vector<int>> blockChunkMapping = tracker.getBlockChunkMapping();

        bool noGates = subcircuits[i].gates.empty();
        bool noSwap2 = (i == 0) || transitions[i-1].swap2.empty();
        bool noSwap1 = (i >= transitions.size()) || transitions[i].swap1.empty();

        if (!(noGates && noSwap1 && noSwap2)) {
            ThreadSafeQueue<BlockData> readQueue(maxBlocksInMemory);
            ThreadSafeQueue<BlockData> processQueue(maxBlocksInMemory);

            std::thread reader([&]() {
                for (int blockIdx = 0; blockIdx < numBlocks; ++blockIdx) {
                    std::unique_lock<std::mutex> lock(memMtx);
                    memCv.wait(lock, []() { return inFlightBlocks < maxBlocksInMemory; });
                    inFlightBlocks++;

                    BlockData block;
                    block.blockIdx = blockIdx;
                    block.chunkIndices = blockChunkMapping[blockIdx];
                    block.bufferIndex = blockIdx % maxBlocksInMemory;

                    void* alignedBuf = state.getAlignedBuffer(block.bufferIndex);
                    state.loadBlock(block.blockIdx, block.chunkIndices, alignedBuf, block.buffer);

                    if (i > 0 && !transitions[i-1].swap2.empty()) {
                        SubCircuit swap2;
                        swap2.gates = scheduleSwaps(transitions[i-1].swap2);
                        swap2.permutation.resize(qubitsPerBlock);
                        std::iota(swap2.permutation.begin(), swap2.permutation.end(), 0);
                        applySubCircuitToBlock(swap2, block.buffer, qubitsPerBlock);
                    }

                    readQueue.push(std::move(block));
                }
                readQueue.setFinished();
            });

            std::thread processor([&]() {
                BlockData block;
                while (readQueue.pop(block)) {
                    if (block.blockIdx == 0)
                        applySubCircuitToBlockDebug(subcircuits[i], block.buffer, qubitsPerBlock);
                    else
                        applySubCircuitToBlock(subcircuits[i], block.buffer, qubitsPerBlock);
                    processQueue.push(std::move(block));
                }
                processQueue.setFinished();
            });

            std::thread writer([&]() {
                BlockData block;
                while (processQueue.pop(block)) {
                    if (i < transitions.size() && !transitions[i].swap1.empty()) {
                        SubCircuit swap1;
                        swap1.gates = scheduleSwaps(transitions[i].swap1);
                        swap1.permutation.resize(qubitsPerBlock);
                        std::iota(swap1.permutation.begin(), swap1.permutation.end(), 0);
                        if (block.blockIdx == 0)
                            applySubCircuitToBlockDebug(swap1, block.buffer, qubitsPerBlock);
                        else
                            applySubCircuitToBlock(swap1, block.buffer, qubitsPerBlock);
                    }

                    void* alignedBuf = state.getAlignedBuffer(block.bufferIndex);
                    state.saveBlock(block.blockIdx, block.chunkIndices, alignedBuf, block.buffer);

                    {
                        std::lock_guard<std::mutex> lock(memMtx);
                        inFlightBlocks--;
                    }
                    memCv.notify_one();
                }
            });

            reader.join();
            processor.join();
            writer.join();
        }

        if (i < transitions.size()) {
            for (int level : transitions[i].swapLevels) {
                int maxLevel = *std::max_element(transitions[i].swapLevels.begin(), transitions[i].swapLevels.end());
                std::vector<int> newChunkMap = tracker.areaSwapShuffle(tracker.chunkMap, level, maxLevel);
                tracker.chunkMap = newChunkMap;
            }
        }
    }
}


