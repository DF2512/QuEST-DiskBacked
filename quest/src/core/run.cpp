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

// Operate in-place on the alignedBuf (no copying)
void applySubCircuitToBlock(const SubCircuit& sub, void* alignedBuf, int qubitsPerBlock) {
    // Interpret memory directly as qcomp amplitudes for the qubits in this block.
    auto* amps = static_cast<qcomp*>(alignedBuf);

    // Create a Qureg that wraps the existing memory
    Qureg tempQureg = createTempQureg(amps, qubitsPerBlock);

    // Apply each gate mapped to the local physical indices
    for (const GateOp& gate : sub.gates) {
        GateOp mappedGate = gate;
        mappedGate.target = logicalToPhysical(sub.permutation, gate.target);
        if (gate.control != -1)
            mappedGate.control = logicalToPhysical(sub.permutation, gate.control);
        GateScheduler::applyGateOpToQureg(mappedGate, tempQureg);
    }
}

std::vector<GateOp> scheduleSwaps(const std::vector<std::pair<int, int>>& swaps) {
    GateScheduler sched;
    for (const auto& [a, b] : swaps) {
        sched.addSwap(a, b);
    }
    auto result = sched.getSchedule();
    return result;
}

void runCircuit(GateScheduler& scheduler, DiskBackedState& state, bool verbose) {
    int numQubits = state.getNumQubits();
    int numLocalQubits = state.getNumQubitsPerBlock();
    int numBlocks = state.getNumBlocks();
    int qubitsPerBlock = numLocalQubits;
    int maxBlocksInMemory = state.getMaxBlocksInMemory();
    PermutationTracker& tracker = state.getPermutationTracker();

    std::vector<SubCircuit> subcircuits =
        scheduler.partitionIntoSubcircuits(numQubits, numLocalQubits, state, verbose);
    std::vector<Transition> transitions = tracker.generateTransitions(subcircuits);

    if (tracker.currentPermutation != subcircuits[0].permutation) {
        return;
    }

    std::atomic<int> blocksCompleted = 0;

    for (size_t i = 0; i < subcircuits.size(); ++i) {
        tracker.currentPermutation = subcircuits[i].permutation;
        std::vector<std::vector<int>> blockChunkMapping = tracker.getBlockChunkMapping();

        bool noGates = subcircuits[i].gates.empty();
        bool noSwap2 = (i == 0) || transitions[i - 1].swap2.empty();
        bool noSwap1 = (i >= transitions.size()) || transitions[i].swap1.empty();

        if (!(noGates && noSwap1 && noSwap2)) {
            blocksCompleted = 0;
            std::mutex doneMtx;
            std::condition_variable doneCv;

            std::mutex memMtx;
            std::condition_variable memCv;
            int inFlightBlocks = 0;

            std::vector<std::atomic<bool>> blockFinished(numBlocks);
            for (auto &f : blockFinished) f = false;

            std::vector<bool> bufferBusy(maxBlocksInMemory, false);
            std::vector<std::queue<int>> bufferQueues(maxBlocksInMemory);
            for (int blockIdx = 0; blockIdx < numBlocks; ++blockIdx) {
                int bufferIndex = blockIdx % maxBlocksInMemory;
                bufferQueues[bufferIndex].push(blockIdx);
            }

            for (int blockIdx = 0; blockIdx < numBlocks; ++blockIdx) {
                int bufferIndex = blockIdx % maxBlocksInMemory;

                std::unique_lock<std::mutex> lock(memMtx);
                memCv.wait(lock, [&]() {
                    return (inFlightBlocks < maxBlocksInMemory &&
                            !bufferBusy[bufferIndex] &&
                            !bufferQueues[bufferIndex].empty() &&
                            bufferQueues[bufferIndex].front() == blockIdx);
                });

                inFlightBlocks++;
                bufferBusy[bufferIndex] = true;
                bufferQueues[bufferIndex].pop();
                lock.unlock();

                std::vector<int> chunkIndices = blockChunkMapping[blockIdx];

                auto currentSubcircuit = subcircuits[i];
                auto currentSwap2 = (i > 0 ? transitions[i-1].swap2 : std::vector<std::pair<int,int>>{});
                auto currentSwap1 = (i < transitions.size() ? transitions[i].swap1 : std::vector<std::pair<int,int>>{});
                void* localBuf = state.getAlignedBuffer(bufferIndex);

                state.loadBlock(blockIdx, chunkIndices,
                    [blockIdx, bufferIndex, localBuf, chunkIndices,
                     currentSubcircuit, currentSwap2, currentSwap1,
                     &state, &blocksCompleted, &doneMtx, &doneCv,
                     &memMtx, &memCv, &bufferBusy, &inFlightBlocks,
                     &blockFinished, numBlocks, qubitsPerBlock]() {

                        if (!currentSwap2.empty()) {
                            SubCircuit swap2;
                            swap2.gates = scheduleSwaps(currentSwap2);
                            swap2.permutation.resize(qubitsPerBlock);
                            std::iota(swap2.permutation.begin(), swap2.permutation.end(), 0);
                            applySubCircuitToBlock(swap2, localBuf, qubitsPerBlock);
                        }

                        applySubCircuitToBlock(currentSubcircuit, localBuf, qubitsPerBlock);

                        if (!currentSwap1.empty()) {
                            SubCircuit swap1;
                            swap1.gates = scheduleSwaps(currentSwap1);
                            swap1.permutation.resize(qubitsPerBlock);
                            std::iota(swap1.permutation.begin(), swap1.permutation.end(), 0);
                            applySubCircuitToBlock(swap1, localBuf, qubitsPerBlock);
                        }

                        state.saveBlock(blockIdx, chunkIndices,
                            [blockIdx, bufferIndex,
                             &memMtx, &memCv, &bufferBusy, &inFlightBlocks,
                             &blocksCompleted, &doneMtx, &doneCv, &blockFinished, numBlocks]() {

                                {
                                    std::lock_guard<std::mutex> lock(memMtx);
                                    inFlightBlocks--;
                                    bufferBusy[bufferIndex] = false;
                                }
                                memCv.notify_all();

                                blockFinished[blockIdx] = true;
                                int done = ++blocksCompleted;
                                if (done == numBlocks) {
                                    std::lock_guard<std::mutex> lock(doneMtx);
                                    doneCv.notify_one();
                                }
                        });
                });
            }

            std::unique_lock<std::mutex> doneLock(doneMtx);
            doneCv.wait(doneLock, [&]() { return blocksCompleted == numBlocks; });
        }

        if (i < transitions.size()) {
            for (int level : transitions[i].swapLevels) {
                int maxLevel = *std::max_element(
                    transitions[i].swapLevels.begin(),
                    transitions[i].swapLevels.end());
                std::vector<int> newChunkMap =
                    tracker.areaSwapShuffle(tracker.chunkMap, level, maxLevel);
                tracker.chunkMap = newChunkMap;
            }
        }
    }
}
