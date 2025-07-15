#pragma once

#include "types.h"
#include <vector>
#include <functional>
class DiskBackedState;

// ─── Supported Gate Types ─────────────────────────────
enum class GateType {
    Single,
    Controlled,
    Swap
};

// ─── A Single Gate Instruction ────────────────────────
struct GateOp {
    GateType type;
    int target;
    int control;  // -1 if not applicable
    std::function<void(qcomp&, qcomp&)> fn;
    std::function<void(std::vector<qcomp>&, int)> bufferFn;  
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
    void addSingleQubitGate(int target, std::function<void(qcomp&, qcomp&)> fn);
    void addControlledGate(int control, int target, std::function<void(qcomp&, qcomp&)> fn);
    void addSwapGate(int qubitA, int qubitB, std::function<void(std::vector<qcomp>&, int)> fn);


    const std::vector<GateOp>& getSchedule() const;

    
    std::vector<SubCircuit> partitionIntoSubcircuits(int numQubits, int numLocalQubits,
                                                 DiskBackedState& state,
                                                 bool verbose = false) const;

private:
    std::vector<GateOp> schedule;
};

