#pragma once

#include "types.h"
#include "qureg.h"
#include <vector>
#include <functional>
#include <utility>
class DiskBackedState;
#define _USE_MATH_DEFINES

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

// ─── Gate Types ────────────────────────────────────
enum class GateType {
    Hadamard,
    Phase,
    S,
    T,
    CNOT,
    ControlledPhase,
    CZ,
    ControlledRK,
    Swap,
    PauliX,
    PauliY,
    PauliZ,
    RotateX,
    RotateY,
    RotateZ,
    SqrtX,
    SqrtY
    // ... add more as needed
};

struct Pattern {
    std::vector<std::pair<int, int>> pairs;
};

// ─── A Single Gate Instruction ────────────────────────
struct GateOp {
    GateType type;
    int target;      // Target qubit (or first target for 2-qubit gates)
    int control;     // Control qubit (or second target for 2-qubit gates, -1 if not used)
    double angle;    // For phase/rotation gates (0.0 if not used)
    int k;           // For ControlledRK (0 if not used)
    // Add more fields as needed for other gate types
};

// ─── A Subcircuit: Local Gates + Permutation Layout ───
struct SubCircuit {
    std::vector<GateOp> gates;
    std::vector<int> permutation;
    std::vector<std::vector<int>> blockChunkMap;
};


// ─── Scheduler Class ──────────────────────────────────
class GateScheduler {
public:
    // Add a Hadamard gate
    void addHadamard(int target) {
        schedule.push_back(GateOp{GateType::Hadamard, target, -1, 0.0, 0});
    }
    // Add a Phase gate
    void addPhase(int target, double angle) {
        schedule.push_back(GateOp{GateType::Phase, target, -1, angle, 0});
    }
    // Add an S gate (π/2 phase)
    void addS(int target) {
        schedule.push_back(GateOp{GateType::S, target, -1, M_PI / 2.0, 0});
    }
    // Add a T gate (π/4 phase)
    void addT(int target) {
        schedule.push_back(GateOp{GateType::T, target, -1, M_PI / 4.0 , 0});
    }
    // Add a Rotate-X gate
    void addRotateX(int target, double angle) {
        schedule.push_back(GateOp{GateType::RotateX, target, -1, angle, 0});
    }
    // Add a Rotate-Y gate
    void addRotateY(int target, double angle) {
        schedule.push_back(GateOp{GateType::RotateY, target, -1, angle, 0});
    }
    // Add a Rotate-Z gate
    void addRotateZ(int target, double angle) {
        schedule.push_back(GateOp{GateType::RotateZ, target, -1, angle, 0});
    }
    // Add a Sqrt-X gate
    void addSqrtX(int target) {
        schedule.push_back(GateOp{GateType::SqrtX, target, -1, M_PI / 2.0, 0});
    }
    // Add a Sqrt-Y gate
    void addSqrtY(int target) {
        schedule.push_back(GateOp{GateType::SqrtY, target, -1, M_PI / 2.0, 0});
    }
    // Add a CNOT gate
    void addCNOT(int control, int target) {
        schedule.push_back(GateOp{GateType::CNOT, target, control, 0.0, 0});
    }
    // Add a Controlled-Phase gate
    void addControlledPhase(int control, int target, double angle) {
        schedule.push_back(GateOp{GateType::ControlledPhase, target, control, angle, 0});
    }
    // Add a CZ gate (controlled-Z)
    void addCZ(int control, int target) {
        schedule.push_back(GateOp{GateType::CZ, target, control, 0.0, 0});
    }
    // Add a Controlled-RK gate
    void addControlledRK(int control, int target, int k) {
        schedule.push_back(GateOp{GateType::ControlledRK, target, control, 2 * M_PI / pow(2, k), k});
    }
    // Add a SWAP gate
    void addSwap(int qubitA, int qubitB) {
        schedule.push_back(GateOp{GateType::Swap, qubitA, qubitB, 0.0, 0});
    }
    // Add a Pauli-X gate
    void addPauliX(int target) {
        schedule.push_back(GateOp{GateType::PauliX, target, -1, 0.0, 0});
    }
    // Add a Pauli-Y gate
    void addPauliY(int target) {
        schedule.push_back(GateOp{GateType::PauliY, target, -1, 0.0, 0});
    }
    // Add a Pauli-Z gate
    void addPauliZ(int target) {
        schedule.push_back(GateOp{GateType::PauliZ, target, -1, 0.0, 0});
    }
    // Add a full QFT
    void addFullQFT(int numQubits) {
        // Apply QFT to all qubits in reverse order
        for (int n = numQubits - 1; n >= 0; n--) {
            addHadamard(n);
            for (int m = 0; m < n; m++) {
                double arg = M_PI / (1 << (m + 1)); // 2^(m+1)
                addControlledPhase(n, n - m - 1, arg);
            }
        }
        // Apply SWAP gates to reverse the order
        int mid = numQubits / 2;
        for (int n = 0; n < mid; n++) {
            addSwap(n, numQubits - 1 - n);
        }
    }
    // Add a Quantum Supremacy Circuit
    void addQSC(int numQubits, int depth); 
    

    // Get the full schedule
    const std::vector<GateOp>& getSchedule() const { return schedule; }

    // Partition schedule into subcircuits (unchanged)
    std::vector<SubCircuit> partitionIntoSubcircuits(int numQubits, int numLocalQubits,
                                                     DiskBackedState& state,
                                                     bool verbose = false) const;

    // Static dispatcher: apply a GateOp to a Qureg
    static void applyGateOpToQureg(const GateOp& op, Qureg& qureg);

private:
    std::vector<GateOp> schedule;
};

