#pragma once

#include "gates.h"
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
#include "gatescheduler.h"
#include "chunkmanager.h" 
#include "operations.h"
#include "qureg.h"

// ─── Dispatcher: Apply a GateOp to a Qureg ─────────────
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
        // ... add more cases as needed
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




