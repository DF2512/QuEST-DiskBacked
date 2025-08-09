#pragma once

#include <numeric>
#include <iostream>
#include <algorithm>
#include "diskbackedstate.h" 
#include <iostream>
#include <numeric>
#include <unordered_set>
#include <set>
#include <unordered_map>
#include <deque>
#include <random>
#include <ctime>
#include "gatescheduler.h"
#include "chunkmanager.h" 
#include "operations.h"
#include "qureg.h"

// ───  Apply a GateOp to a Qureg ─────────────
void GateScheduler::applyGateOpToQureg(const GateOp& op, Qureg& qureg) {
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
        case GateType::RotateX:
            applyRotateX(qureg, op.target, op.angle);
            break;
        case GateType::RotateY:
            applyRotateY(qureg, op.target, op.angle);
            break;
        case GateType::RotateZ:
            applyRotateZ(qureg, op.target, op.angle);
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
        case GateType::PauliX:
            applyPauliX(qureg, op.target);
            break;
        case GateType::PauliY:
            applyPauliY(qureg, op.target);
            break;
        case GateType::PauliZ:
            applyPauliZ(qureg, op.target);
            break;
        // ... add more gates as needed
        default:
            std::cerr << "[GateScheduler] Unknown gate type in dispatcher!\n";
            break;
    }
}

// ─── Helper: Check if a Gate is Local ─────────────
bool isGateLocal(const GateOp& gate, const std::vector<int>& layout, int numLocalQubits) {
    int n = layout.size();
    std::vector<int> logicalToPhysical(n);
    for (int i = 0; i < n; ++i)
        logicalToPhysical[layout[i]] = i;

    int targMapped = logicalToPhysical[gate.target];
    int ctrlMapped = (gate.control != -1) ? logicalToPhysical[gate.control] : -1;

    return targMapped < numLocalQubits && (ctrlMapped == -1 || ctrlMapped < numLocalQubits);
}

// ─── Helper: Lookahead Swap Optimisation ────────────────────────────────
std::vector<int> optimisePermutationWithLookahead(
    const std::vector<int>& currentPerm,
    const std::vector<GateOp>& upcomingGates,
    int numLocalQubits,
    int maxPermutableQubits,
    bool verbose = false
) {
    int n = static_cast<int>(currentPerm.size());

    std::vector<int> posInPerm(n);
    for (int i = 0; i < n; ++i)
        posInPerm[currentPerm[i]] = i;

    std::unordered_map<int, int> firstUse;
    for (int i = 0; i < static_cast<int>(upcomingGates.size()); ++i) {
        const GateOp& g = upcomingGates[i];
        if (!firstUse.count(g.target))
            firstUse[g.target] = i;
        if (g.control != -1 && !firstUse.count(g.control))
            firstUse[g.control] = i;
    }

    std::vector<std::pair<int, int>> localUsed, nonLocalUsed;
    for (int q = 0; q < n; ++q) {
        int pos = posInPerm[q];
        int usage = firstUse.count(q) ? firstUse[q] : INT_MAX;
        if (pos < numLocalQubits)
            localUsed.emplace_back(usage, q);
        else
            nonLocalUsed.emplace_back(usage, q);
    }

    std::sort(localUsed.begin(), localUsed.end(), std::greater<>());
    std::sort(nonLocalUsed.begin(), nonLocalUsed.end());

    int swaps = std::min({maxPermutableQubits, (int)localUsed.size(), (int)nonLocalUsed.size()});

    std::vector<int> toEvict, toImport;
    for (int i = 0; i < swaps; ++i) {
        toEvict.push_back(localUsed[i].second);
        toImport.push_back(nonLocalUsed[i].second);
    }

    if (verbose) {
        std::cout << "[Permute] Lookahead decided to evict: ";
        for (int q : toEvict) std::cout << q << " ";
        std::cout << "\n[Permute] Lookahead decided to import: ";
        for (int q : toImport) std::cout << q << " ";
        std::cout << "\n";
    }

    std::vector<int> newPerm = currentPerm;
    for (int i = 0; i < numLocalQubits; ++i) {
        if (std::find(toEvict.begin(), toEvict.end(), newPerm[i]) != toEvict.end()) {
            for (int j = numLocalQubits; j < n; ++j) {
                if (std::find(toImport.begin(), toImport.end(), newPerm[j]) != toImport.end()) {
                    std::swap(newPerm[i], newPerm[j]);
                    break;
                }
            }
        }
    }

    if (verbose) {
        std::cout << "[Permute] Final permutation: ";
        for (int q : newPerm) std::cout << q << " ";
        std::cout << "\n";
    }

    return newPerm;
}

// ─── Full Partitioning + Lookahead Optimisation ─────────────────────────
std::vector<SubCircuit> GateScheduler::partitionIntoSubcircuits(
    int numQubits,
    int numLocalQubits,
    DiskBackedState& state,
    bool verbose
) const {
    std::vector<SubCircuit> result;
    std::deque<GateOp> remaining(schedule.begin(), schedule.end());
    std::vector<int> perm(numQubits);
    std::iota(perm.begin(), perm.end(), 0);
    // PermutationTracker& tracker = state.getPermutationTracker(); // Remove this

    int subIndex = 0;

    while (!remaining.empty()) {
        std::vector<GateOp> localGates;

        if (verbose) {
            std::cout << "\n[Partition] Subcircuit #" << subIndex << " begins\n";
        }

        // 1. Drain local gates
        while (!remaining.empty() && isGateLocal(remaining.front(), perm, numLocalQubits)) {
            if (verbose) {
                const auto& g = remaining.front();
                std::cout << "  [LocalGate] Q" << g.target;
                if (g.control != -1) std::cout << " (ctrl Q" << g.control << ")";
                std::cout << "\n";
            }

            localGates.push_back(remaining.front());
            remaining.pop_front();
        }

        SubCircuit sub;
        sub.gates = localGates;
        sub.permutation = perm; // Store the current permutation

        result.push_back(sub);

        if (verbose) {
            std::cout << "[Partition] Subcircuit #" << subIndex
                      << " stored with " << localGates.size() << " gates\n";
        }

        ++subIndex;

        if (remaining.empty()) break;

        if (verbose) {
            const auto& g = remaining.front();
            std::cout << "[Partition] Non-local gate encountered on Q" << g.target;
            if (g.control != -1)
                std::cout << " (ctrl Q" << g.control << ")";
            std::cout << "\n";
        }

        std::vector<GateOp> nextGates;
        for (const auto& gate : remaining) {
            nextGates.push_back(gate);
            if (nextGates.size() >= 32) break;
        }
        
        
        // Update permutation for lookahead, but don't update the tracker
        perm = optimisePermutationWithLookahead(
            perm, nextGates,
            numLocalQubits,
            state.getMaxPermutableQubits(),
            verbose
        );
    }

    // NOTE: This is a temporary bandaid fix. There is a bug where if the last qubit is correct
    // on the penultimate subcircuit, it will not be restored correctly. Root cause is
    // unknown. More accurately, it seems that the issue happens when the swap doesnt involve the last qubit.
    // TODO: Fix this properly. 
    // Check if the final subcircuit ends with the highest qubit
    bool needsTwoStepRestoration = false;
    if (!result.empty()) {
        const auto& lastSub = result.back();
        if (lastSub.permutation.back() == numQubits - 1) {
            needsTwoStepRestoration = true;
        }
    }

    if (needsTwoStepRestoration) {
        // Add first restoration step: identity with last 2 qubits swapped
        SubCircuit firstRestoration;
        firstRestoration.gates = {}; // Empty gates
        firstRestoration.permutation.resize(numQubits);
        std::iota(firstRestoration.permutation.begin(), firstRestoration.permutation.end(), 0);
        // Swap the last two qubits
        std::swap(firstRestoration.permutation[numQubits - 1], firstRestoration.permutation[numQubits - 2]);
        result.push_back(firstRestoration);
        
        // Add second restoration step: true identity
        SubCircuit secondRestoration;
        secondRestoration.gates = {}; // Empty gates
        secondRestoration.permutation.resize(numQubits);
        std::iota(secondRestoration.permutation.begin(), secondRestoration.permutation.end(), 0);
        result.push_back(secondRestoration);
    } else {
        // Add final empty subcircuit with identity permutation
        SubCircuit finalSub;
        finalSub.gates = {}; // Empty gates
        finalSub.permutation.resize(numQubits);
        std::iota(finalSub.permutation.begin(), finalSub.permutation.end(), 0); // Identity permutation
        result.push_back(finalSub);
    }
    

    if (verbose) {
        std::cout << "\n[Partition] Total subcircuits: " << result.size() << "\n";
        for (size_t i = 0; i < result.size(); ++i) {
            const auto& sub = result[i];
            std::cout << "  Subcircuit #" << i << ": gates = " << sub.gates.size()
                      << ", permutation = [ ";
            for (int q : sub.permutation) std::cout << q << " ";
            std::cout << "]\n";
        }
    }

    return result;
}

// ─── Quantum Supremacy Circuit ─────────────────────────────
void GateScheduler::addQSC(int numQubits, int depth) {
    // Initialize random number generator
    std::mt19937 rng(std::time(nullptr));
    std::uniform_int_distribution<int> patternDist(0, 7);
    std::uniform_int_distribution<int> singleQubitDist(0, 2); // 0=SqrtX, 1=SqrtY, 2=T
    
    // Track previous gates for each qubit
    std::vector<GateType> previousGates(numQubits, GateType::Hadamard); // Initial Hadamard cycle
    std::vector<bool> hadCZInPreviousCycle(numQubits, false);
    
    // 1. Start with a cycle of Hadamard gates (0 clock cycle)
    for (int q = 0; q < numQubits; q++) {
        addHadamard(q);
        hadCZInPreviousCycle[q] = false; // No CZ gates in the initial cycle
    }
    
    // 2. Repeat for d clock cycles
    for (int cycle = 1; cycle <= depth; cycle++) {
        // Eight CZ patterns (similar to Fig. 6)
        std::vector<Pattern> patterns = {
            {{ {2,3}, {6,7}, {10,11}, {14,15}, {18,19}, {22,23}, {26,27}, {30,31}, {34,35} }},
            {{ {0,1}, {4,5}, {8,9}, {12,13}, {16,17}, {20,21}, {24,25}, {28,29}, {32,33} }},
            {{ {7,13}, {9,15}, {11,17}, {18,24}, {20,26}, {22,28} }},
            {{ {6,12}, {8,14}, {10,16}, {19,25}, {21,27}, {23,29} }},
            {{ {3,4}, {7,8}, {15,16}, {19,20}, {27,28}, {31,32} }},
            {{ {1,2}, {9,10}, {13,14}, {21,22}, {25,26}, {33,34} }},
            {{ {0,6}, {2,8}, {4,10}, {13,19}, {15,21}, {17,23}, {24,30}, {26,32}, {28,34} }},
            {{ {1,7}, {3,9}, {5,11}, {12,18}, {14,20}, {16,22}, {25,31}, {27,33}, {29,35} }}
        };
        
        // Select pattern randomly, ensuring it doesn't conflict with rules
        Pattern selectedPattern;
        bool validPattern = false;
        int attempts = 0;
        const int maxAttempts = 100;
        
        while (!validPattern && attempts < maxAttempts) {
            int patternIdx = patternDist(rng);
            selectedPattern = patterns[patternIdx];
            
            // Check if this pattern conflicts with the "no CNOT twice in a row" rule
            validPattern = true;
            for (const auto& pair : selectedPattern.pairs) {
                // Check if either qubit had a CZ gate in the previous cycle
                if (hadCZInPreviousCycle[pair.first] || hadCZInPreviousCycle[pair.second]) {
                    validPattern = false;
                    break;
                }
            }
            attempts++;
        }
        
        // If we couldn't find a valid pattern, use the first one (fallback)
        // A better fallback will be implemented later
        if (!validPattern) {
            selectedPattern = patterns[0];
        }
        
        // Apply CZ gates from the selected pattern
        std::vector<bool> hasCZThisCycle(numQubits, false);
        for (const auto& pair : selectedPattern.pairs) {
            if (pair.first < numQubits && pair.second < numQubits) {
                addCZ(pair.first, pair.second);
                hasCZThisCycle[pair.first] = true;
                hasCZThisCycle[pair.second] = true;
            }
        }
        
        // Apply single-qubit gates to qubits not occupied by CZ gates
        for (int q = 0; q < numQubits; q++) {
            if (!hasCZThisCycle[q]) {
                // Rule: Place a gate at qubit q only if this qubit is occupied by a CZ gate in the previous cycle
                if (hadCZInPreviousCycle[q]) {
                    // Rule: Place a T gate if there are no single-qubit gates in previous cycles except initial Hadamard
                    bool hasHadSingleQubitGate = false;
                    for (int prevCycle = 1; prevCycle < cycle; prevCycle++) {
                        if (previousGates[q] != GateType::Hadamard && previousGates[q] != GateType::CZ) {
                            hasHadSingleQubitGate = true;
                            break;
                        }
                    }
                    
                    GateType selectedGate;
                    if (!hasHadSingleQubitGate) {
                        // Must place T gate
                        selectedGate = GateType::T;
                    } else {
                        // Randomly select from {X1/2, Y1/2, T}
                        int gateChoice = singleQubitDist(rng);
                        switch (gateChoice) {
                            case 0: selectedGate = GateType::SqrtX; break;
                            case 1: selectedGate = GateType::SqrtY; break;
                            case 2: selectedGate = GateType::T; break;
                        }
                        
                        // Rule: Any gate at qubit q should be different from the gate at qubit q in the previous cycle
                        if (selectedGate == previousGates[q]) {
                            // Try a different gate
                            if (previousGates[q] == GateType::SqrtX) {
                                selectedGate = (singleQubitDist(rng) < 1) ? GateType::SqrtY : GateType::T;
                            } else if (previousGates[q] == GateType::SqrtY) {
                                selectedGate = (singleQubitDist(rng) < 1) ? GateType::SqrtX : GateType::T;
                            } else { 
                                selectedGate = (singleQubitDist(rng) < 1) ? GateType::SqrtX : GateType::SqrtY;
                            }
                        }
                    }
                    
                    // Apply the selected gate
                    switch (selectedGate) {
                        case GateType::SqrtX:
                            addSqrtX(q);
                            break;
                        case GateType::SqrtY:
                            addSqrtY(q);
                            break;
                        case GateType::T:
                            addT(q);
                            break;
                        default:
                            break;
                    }
                    
                    previousGates[q] = selectedGate;
                }
            }
        }
        
        // Update tracking for next cycle
        hadCZInPreviousCycle = hasCZThisCycle;
    }
    
    // 3. Repeat the Hadamards at the end
    for (int q = 0; q < numQubits; q++) {
        addHadamard(q);
    }
}


