/** @file
 * A minimum C++ example of running
 * QuEST with disk-backed state management.
 * 
 * @author Tyson Jones
*/

#include "quest.h"
#include "diskbackedstate.h"
#include "gatescheduler.h"
#include "run.h"
#include "chunkmanager.h"
#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include <random>
#include <ctime>
#include <complex>
#include <algorithm>
#include <iomanip>
#include <numeric>
#include <exception>
#include <omp.h>
#include <set> 
#include <limits>

struct RunLog {
    int runIdx;
    int numSubcircuits;
    bool anySwap1;
    bool anySwap2;
    bool anyChunkSwap;
    bool success;
    bool segfaultOrException;
    // Measurement outcome info
    int diskOutcome = -1;
    int regOutcome = -1;
    bool measurementMatch = false;
    // Restoration debug info
    int restorationSwap1Count = 0;
    int restorationSwap2Count = 0;
    int restorationAreaSwapCount = 0;
    size_t firstMismatchIdx = (size_t)-1;
    int diskQubit = -1;
    int regQubit = -1;
    int totalAreaShuffleMismatches = 0;
    bool bagsMatch = false;
    std::vector<int> permutationBeforeRestoration;
    // Final transition info
    std::vector<std::pair<int, int>> finalSwap1;
    std::vector<std::pair<int, int>> finalSwap2;
    std::vector<int> finalSwapLevels;
    std::vector<int> finalInterimTarget;
    // Power-of-2 amplitude mapping
    std::vector<int> powerOf2Mapping;
};

// Helper: For each power-of-2 index, find where that amplitude is in the disk-backed vector
std::vector<int> getPowerOf2AmplitudeMapping(const std::vector<qcomp>& ampsRegular, const std::vector<qcomp>& ampsDisk, double eps) {
    int n = static_cast<int>(std::log2(ampsRegular.size()));
    std::vector<int> mapping(n, -1);
    for (int k = 0; k < n; ++k) {
        size_t idx = 1ULL << k;
        if (idx >= ampsRegular.size()) break;
        qcomp val = ampsRegular[idx];
        // Find in disk-backed
        int found = -1;
        for (size_t j = 0; j < ampsDisk.size(); ++j) {
            if (std::abs(ampsDisk[j] - val) < eps) {
                found = static_cast<int>(j);
                break;
            }
        }
        mapping[k] = found;
    }
    return mapping;
}

int main() {
    initQuESTEnv();
    reportQuESTEnv();

    // === CONFIGURE NUMBER OF RUNS HERE ===
    const int numRuns = 50; 
    int numSuccesses = 0;
    std::vector<RunLog> logs;

    for (int run = 1; run <= numRuns; ++run) {
        std::cout << "\n====================\n";
        std::cout << "Run " << run << " of " << numRuns << "\n";
        std::cout << "====================\n";

        RunLog log = {};
        log.runIdx = run;
        log.numSubcircuits = 0;
        log.anySwap1 = false;
        log.anySwap2 = false;
        log.anyChunkSwap = false;
        log.success = false;
        log.segfaultOrException = false;

        const int numQubits = 20;
        const int numBlocks = 8;
        const int chunksPerBlock = 8;
        std::vector<std::string> diskRoots = {
            "C:/quantum_chunks0",
            "D:/quantum_chunks1"
        };

        // 1. Create and init Qureg
        Qureg qureg = createForcedQureg(numQubits);
        initRandomPureState(qureg);

        // 2. Copy amplitudes to vector
        std::vector<qcomp> amps(qureg.numAmps);
        getQuregAmps(amps.data(), qureg, 0, qureg.numAmps);

        // 3. Create disk-backed state and save amplitudes chunk by chunk
        DiskBackedState diskState(numQubits, numBlocks, chunksPerBlock, diskRoots);
        size_t ampsPerChunk = diskState.getAmpsPerChunk();
        size_t numChunks = diskState.getNumChunks();
        for (size_t chunk = 0; chunk < numChunks; ++chunk) {
            std::vector<qcomp> chunkBuf(amps.begin() + chunk * ampsPerChunk, amps.begin() + (chunk + 1) * ampsPerChunk);
            diskState.saveChunk(chunk, chunkBuf);
        }
        // Print and compare initial amplitudes between regular and disk-backed state
        std::vector<qcomp> ampsDiskInit(qureg.numAmps);
        for (size_t chunk = 0; chunk < numChunks; ++chunk) {
            std::vector<qcomp> chunkBuf;
            diskState.loadChunk(chunk, chunkBuf);
            std::copy(chunkBuf.begin(), chunkBuf.end(), ampsDiskInit.begin() + chunk * ampsPerChunk);
        }
        int mismatchCount = 0;
        double epsInit = 1e-12;
        for (size_t i = 0; i < qureg.numAmps; ++i) {
            if (std::abs(amps[i] - ampsDiskInit[i]) > epsInit) {
                if (mismatchCount < 10)
                    std::cout << "Initial mismatch at index " << i << ": " << amps[i] << " vs " << ampsDiskInit[i] << std::endl;
                ++mismatchCount;
            }
        }
        if (mismatchCount == 0) {
            std::cout << "SUCCESS: All initial amplitudes match between regular and disk-backed state." << std::endl;
        } else {
            std::cout << "FAILURE: " << mismatchCount << " initial amplitude mismatches found." << std::endl;
        }
        amps.clear();
        amps.shrink_to_fit();

        // 4. Build a schedule applying each gate to all qubits exactly once, in random order
        GateScheduler scheduler;
        std::vector<int> qubitOrder(numQubits);
        std::iota(qubitOrder.begin(), qubitOrder.end(), 0);
        std::mt19937 rng(static_cast<unsigned>(std::time(nullptr)) + run); // ensure different seed per run
        std::shuffle(qubitOrder.begin(), qubitOrder.end(), rng);
        // Uncomment ONE block below to test a single gate type at a time
        
        /**/
        // --- Hadamard ---
        for (int idx = 0; idx < numQubits; ++idx) {
            int q = qubitOrder[idx];
            scheduler.addHadamard(q);
        }
        
        // --- Phase ---
        for (int idx = 0; idx < numQubits; ++idx) {
            int q = qubitOrder[idx];
            scheduler.addPhase(q, M_PI/4.0); // or any fixed angle
        }
        
        // --- S ---
        for (int idx = 0; idx < numQubits; ++idx) {
            int q = qubitOrder[idx];
            scheduler.addS(q);
        }
        
        // --- T ---
        for (int idx = 0; idx < numQubits; ++idx) {
            int q = qubitOrder[idx];
            scheduler.addT(q);
        }
        
        // --- SqrtX ---
        for (int idx = 0; idx < numQubits; ++idx) {
            int q = qubitOrder[idx];
            scheduler.addSqrtX(q);
        }
        
        // --- SqrtY ---
        for (int idx = 0; idx < numQubits; ++idx) {
            int q = qubitOrder[idx];
            scheduler.addSqrtY(q);
        }
        
        // --- CNOT ---
        std::vector<int> cnotOrder(numQubits-1);
        std::iota(cnotOrder.begin(), cnotOrder.end(), 0);
        std::shuffle(cnotOrder.begin(), cnotOrder.end(), rng);
        for (int idx = 0; idx < numQubits-1; ++idx) {
            int q = cnotOrder[idx];
            scheduler.addCNOT(q, q+1);
        }
        
        // --- ControlledPhase ---
        std::vector<int> cpOrder(numQubits-1);
        std::iota(cpOrder.begin(), cpOrder.end(), 0);
        std::shuffle(cpOrder.begin(), cpOrder.end(), rng);
        for (int idx = 0; idx < numQubits-1; ++idx) {
            int q = cpOrder[idx];
            scheduler.addControlledPhase(q, q+1, M_PI/4.0);
        }
        
        // --- CZ ---
        std::vector<int> czOrder(numQubits-1);
        std::iota(czOrder.begin(), czOrder.end(), 0);
        std::shuffle(czOrder.begin(), czOrder.end(), rng);
        for (int idx = 0; idx < numQubits-1; ++idx) {
            int q = czOrder[idx];
            scheduler.addCZ(q, q+1);
        }
        
        // --- ControlledRK ---
        std::vector<int> crkOrder(numQubits-1);
        std::iota(crkOrder.begin(), crkOrder.end(), 0);
        std::shuffle(crkOrder.begin(), crkOrder.end(), rng);
        for (int idx = 0; idx < numQubits-1; ++idx) {
            int q = crkOrder[idx];
            scheduler.addControlledRK(q, q+1, 2);
        }
        
        // --- Swap ---
        std::vector<int> swapOrder(numQubits-1);
        std::iota(swapOrder.begin(), swapOrder.end(), 0);
        std::shuffle(swapOrder.begin(), swapOrder.end(), rng);
        for (int idx = 0; idx < numQubits-1; ++idx) {
            int q = swapOrder[idx];
            scheduler.addSwap(q, q+1);
        }
        
        // --- Quantum Supremacy Circuit ---
        //scheduler.addQSC(numQubits, 8); // 3 cycles of quantum supremacy circuit
        
        // --- Full QFT ---
        // scheduler.addFullQFT(numQubits);

        // 5. Apply the schedule to the regular Qureg using direct gate functions
        const auto& schedule = scheduler.getSchedule();
        for (const auto& op : schedule) {
            switch (op.type) {
                case GateType::Hadamard:
                    applyHadamard(qureg, op.target);
                    break;
                case GateType::Phase:
                    applyPhaseShift(qureg, op.target, op.angle);
                    break;
                case GateType::S:
                    applyS(qureg, op.target);
                    break;
                case GateType::T:
                    applyT(qureg, op.target);
                    break;
                case GateType::SqrtX:
                    applyRotateX(qureg, op.target, op.angle);
                    break;
                case GateType::SqrtY:
                    applyRotateY(qureg, op.target, op.angle);
                    break;
                case GateType::CNOT:
                    applyControlledPauliX(qureg, op.control, op.target);
                    break;
                case GateType::ControlledPhase:
                    applyTwoQubitPhaseShift(qureg, op.control, op.target, op.angle);
                    break;
                case GateType::CZ:
                    applyControlledPauliZ(qureg, op.control, op.target);
                    break;
                case GateType::ControlledRK:
                    applyTwoQubitPhaseShift(qureg, op.control, op.target, op.angle);
                    break;
                case GateType::Swap:
                    applySwap(qureg, op.target, op.control);
                    break;
                default:
                    std::cerr << "[Regular Qureg] Unknown gate type in schedule!\n";
                    break;
            }
        }

        // applyFullQuantumFourierTransform(qureg);


        // 6. Apply the schedule to the disk-backed state
        try {
            // Intercept subcircuit/transition info
            
            std::vector<SubCircuit> subcircuits = scheduler.partitionIntoSubcircuits(numQubits, diskState.getNumQubitsPerBlock(), diskState, false);
            log.numSubcircuits = (int)subcircuits.size();
            
            // Capture the permutation from the second-to-last subcircuit
            if (subcircuits.size() >= 2) {
                log.permutationBeforeRestoration = subcircuits[subcircuits.size() - 3].permutation;
            }
            
            // Generate transitions
            PermutationTracker& tracker = diskState.getPermutationTracker();
            std::vector<Transition> transitions = tracker.generateTransitions(subcircuits);
            for (const auto& t : transitions) {
                if (!t.swap1.empty()) log.anySwap1 = true;
                if (!t.swap2.empty()) log.anySwap2 = true;
                if (!t.swapLevels.empty()) log.anyChunkSwap = true;
            }
            
            // Capture final transition info (the last transition to identity)
            if (!transitions.empty()) {
                const auto& finalTransition = transitions.back();
                log.finalSwap1 = finalTransition.swap1;
                log.finalSwap2 = finalTransition.swap2;
                log.finalSwapLevels = finalTransition.swapLevels;
                log.finalInterimTarget = finalTransition.interimTarget;
            }
            
            runCircuit(scheduler, diskState, false);

            
            
        } catch (const std::exception& e) {
            std::cout << "[EXCEPTION] " << e.what() << std::endl;
            log.segfaultOrException = true;
        } catch (...) {
            std::cout << "[EXCEPTION] Unknown exception or segfault!" << std::endl;
            log.segfaultOrException = true;
        }
        // Measure qubit 0
        std::cout << "Measuring qubit 0 for disk-backed state" << std::endl;
        int diskOutcome = diskState.diskBacked_applyQubitMeasurement(0);
        std::cout << "Disk-backed state measurement outcome: " << diskOutcome << std::endl;
        std::cout << "Measuring qubit 0 for regular Qureg" << std::endl;
        int regOutcome = applyQubitMeasurement(qureg, 0);
        std::cout << "Regular Qureg measurement outcome: " << regOutcome << std::endl;
        // Store measurement outcomes in log
        log.diskOutcome = diskOutcome;
        log.regOutcome = regOutcome;
        log.measurementMatch = (diskOutcome == regOutcome);
        
        // 7. Extract amplitudes from both
        std::vector<qcomp> ampsRegular(qureg.numAmps);
        getQuregAmps(ampsRegular.data(), qureg, 0, qureg.numAmps);

        std::vector<qcomp> ampsDisk(qureg.numAmps);
        // Use the current chunk map from the permutation tracker
        const auto& chunkMap = diskState.getPermutationTracker().getCurrentChunkMap();
        for (size_t logicalChunk = 0; logicalChunk < numChunks; ++logicalChunk) {
            size_t physicalChunk = chunkMap[logicalChunk];
            std::vector<qcomp> chunkBuf;
            diskState.loadChunk(physicalChunk, chunkBuf);
            std::copy(chunkBuf.begin(), chunkBuf.end(), ampsDisk.begin() + logicalChunk * ampsPerChunk);
        }

        // 8. Compare amplitudes
        bool allMatch = true;
        double eps = 1e-12;
        size_t mismatchIdx = ampsRegular.size();
        for (size_t i = 0; i < ampsRegular.size(); ++i) {
            if (std::abs(ampsRegular[i] - ampsDisk[i]) > eps) {
                std::cout << "Mismatch at index " << i << ": " << ampsRegular[i] << " vs " << ampsDisk[i] << std::endl;
                allMatch = false;
                mismatchIdx = i;
                break;
            }
        }
        
        // Check if bags match (sets of amplitudes)
        std::set<std::pair<double, double>> regularBag, diskBag;
        for (size_t i = 0; i < ampsRegular.size(); ++i) {
            regularBag.insert({std::real(ampsRegular[i]), std::imag(ampsRegular[i])});
            diskBag.insert({std::real(ampsDisk[i]), std::imag(ampsDisk[i])});
        }
        log.bagsMatch = (regularBag == diskBag);
        
        if (!allMatch) {
            log.firstMismatchIdx = mismatchIdx;
            if (mismatchIdx > 0) {
                // Find the logical qubit for this amplitude index
                int diskQubit = -1, regQubit = -1;
                if (mismatchIdx > 0) {
                    diskQubit = static_cast<int>(std::log2(mismatchIdx));
                    // Find where the disk-backed value is in the regular amplitudes
                    qcomp diskVal = ampsDisk[mismatchIdx];
                    size_t foundIdx = ampsRegular.size();
                    #pragma omp parallel for
                    for (size_t j = 0; j < ampsRegular.size(); ++j) {
                        if (std::abs(ampsRegular[j] - diskVal) < eps) {
                            #pragma omp critical
                            {
                                if (foundIdx == ampsRegular.size() || j < foundIdx) foundIdx = j;
                            }
                        }
                    }
                    if (foundIdx < ampsRegular.size()) {
                        regQubit = static_cast<int>(std::log2(foundIdx));
                    }
                }
                log.diskQubit = diskQubit;
                log.regQubit = regQubit;
            }
        }
        
        // Store power-of-2 amplitude mapping for every run
        log.powerOf2Mapping = getPowerOf2AmplitudeMapping(ampsRegular, ampsDisk, eps);
        
        if (allMatch) {
            std::cout << "SUCCESS: All amplitudes match between regular and disk-backed state after all supported gates." << std::endl;
            ++numSuccesses;
            log.success = true;
        } else {
            // On failure, search ampsRegular for the disk-backed value at the mismatch index (using OpenMP)
            if (mismatchIdx < ampsDisk.size()) {
                qcomp diskVal = ampsDisk[mismatchIdx];
                size_t foundIdx = ampsRegular.size();
                #pragma omp parallel for
                for (size_t j = 0; j < ampsRegular.size(); ++j) {
                    if (std::abs(ampsRegular[j] - diskVal) < eps) {
                        #pragma omp critical
                        {
                            if (foundIdx == ampsRegular.size() || j < foundIdx) foundIdx = j;
                        }
                    }
                }
                if (foundIdx < ampsRegular.size()) {
                    std::cout << "[DEBUG] The value at disk-backed index " << mismatchIdx << " (" << diskVal << ") is found at regular index " << foundIdx << ".\n";
                } else {
                    std::cout << "[DEBUG] The value at disk-backed index " << mismatchIdx << " (" << diskVal << ") is NOT found in the regular amplitudes.\n";
                }
            }
            log.success = false;
            //std::exit(1); // DISABLED: do not exit on failure
        }

        destroyQureg(qureg);
        logs.push_back(log);
    }
    std::cout << "\n====================\n";
    std::cout << "Test summary: " << numSuccesses << " / " << numRuns << " runs succeeded.\n";
    std::cout << "====================\n";
    std::cout << "\nRun-by-run log:\n";
    std::cout << "Idx | Subcircuits | Swap1 | Swap2 | ChunkSwap | Success | Exception | RSwap1 | RSwap2 | RArea | FailIdx | DiskQ | RegQ | AreaMismatch | BagsMatch | DiskOut | RegOut | MeasMatch\n";
    for (const auto& log : logs) {
        std::cout << std::setw(3) << log.runIdx << " | "
                  << std::setw(11) << log.numSubcircuits << " | "
                  << std::setw(5) << (log.anySwap1 ? "Y" : "N") << " | "
                  << std::setw(5) << (log.anySwap2 ? "Y" : "N") << " | "
                  << std::setw(9) << (log.anyChunkSwap ? "Y" : "N") << " | "
                  << std::setw(7) << (log.success ? "Y" : "N") << " | "
                  << std::setw(9) << (log.segfaultOrException ? "Y" : "N") << " | "
                  << std::setw(6) << log.restorationSwap1Count << " | "
                  << std::setw(6) << log.restorationSwap2Count << " | "
                  << std::setw(5) << log.restorationAreaSwapCount << " | "
                  << std::setw(7) << (log.firstMismatchIdx == (size_t)-1 ? -1 : (int)log.firstMismatchIdx) << " | "
                  << std::setw(5) << log.diskQubit << " | "
                  << std::setw(5) << log.regQubit << " | "
                  << std::setw(12) << log.totalAreaShuffleMismatches << " | "
                  << std::setw(8) << (log.bagsMatch ? "Y" : "N") << " | "
                  << std::setw(7) << log.diskOutcome << " | "
                  << std::setw(7) << log.regOutcome << " | "
                  << std::setw(9) << (log.measurementMatch ? "Y" : "N") << "\n";
    }
    std::cout << "\n====================\n";
    std::cout << "Permutations before restoration:\n";
    std::cout << "====================\n";
    for (const auto& log : logs) {
        std::cout << "Run " << std::setw(2) << log.runIdx << " (" << (log.success ? "SUCCESS" : "FAILURE") << "): ";
        if (log.permutationBeforeRestoration.empty()) {
            std::cout << "No restoration needed (already identity)\n";
        } else {
            for (size_t i = 0; i < log.permutationBeforeRestoration.size(); ++i) {
                std::cout << log.permutationBeforeRestoration[i];
                if (i < log.permutationBeforeRestoration.size() - 1) std::cout << " ";
            }
            std::cout << "\n";
        }
    }
    
    std::cout << "\n====================\n";
    std::cout << "Final transition information (last transition to identity):\n";
    std::cout << "====================\n";
    for (size_t runIdx = 0; runIdx < logs.size(); ++runIdx) {
        const auto& log = logs[runIdx];
        std::cout << "Run " << std::setw(2) << log.runIdx << " (" << (log.success ? "SUCCESS" : "FAILURE") << "):\n";
        std::cout << "  Final Swap1 (" << log.finalSwap1.size() << " swaps): ";
        for (const auto& [i, j] : log.finalSwap1) {
            std::cout << "(" << i << "," << j << ") ";
        }
        std::cout << "\n";
        std::cout << "  Final Swap2 (" << log.finalSwap2.size() << " swaps): ";
        for (const auto& [i, j] : log.finalSwap2) {
            std::cout << "(" << i << "," << j << ") ";
        }
        std::cout << "\n";
        std::cout << "  Final SwapLevels (" << log.finalSwapLevels.size() << " levels): ";
        for (int level : log.finalSwapLevels) {
            std::cout << level << " ";
        }
        std::cout << "\n";
        
    }
    
    std::cout << "\n====================\n";
    std::cout << "Power-of-2 amplitude mapping (regular idx 2^k -> disk-backed idx):\n";
    std::cout << "====================\n";
    for (const auto& log : logs) {
        std::cout << "Run " << std::setw(2) << log.runIdx << " (" << (log.success ? "SUCCESS" : "FAILURE") << "): ";
        
        // First, save all the log2 values to a vector
        std::vector<int> log2Values(log.powerOf2Mapping.size(), -1);
        for (size_t k = 0; k < log.powerOf2Mapping.size(); ++k) {
            if (log.powerOf2Mapping[k] >= 0) {
                log2Values[k] = static_cast<int>(std::log2(log.powerOf2Mapping[k]));
            }
        }
        
        // Create correctMapping by inverting log2Values
        std::vector<int> correctMapping(log2Values.size(), -1);
        for (size_t i = 0; i < log2Values.size(); ++i) {
            if (log2Values[i] >= 0 && log2Values[i] < static_cast<int>(correctMapping.size())) {
                correctMapping[log2Values[i]] = static_cast<int>(i);
            }
        }
        
        // Print correctMapping instead of log2Values
        for (size_t k = 0; k < correctMapping.size(); ++k) {
            std::cout << correctMapping[k];
            if (k < correctMapping.size() - 1) std::cout << " ";
        }
        std::cout << "\n";
        
        // Compare with permutationBeforeRestoration
        if (!log.permutationBeforeRestoration.empty() && correctMapping.size() == log.permutationBeforeRestoration.size()) {
            bool match = true;
            for (size_t i = 0; i < correctMapping.size(); ++i) {
                if (correctMapping[i] != log.permutationBeforeRestoration[i]) {
                    match = false;
                    break;
                }
            }
            //std::cout << "  Permutation match: " << (match ? "YES" : "NO") << "\n";
        //} else {
            //std::cout << "  Permutation match: CANNOT_COMPARE\n";
        }
    }
    
    finalizeQuESTEnv();
    return 0;
}
