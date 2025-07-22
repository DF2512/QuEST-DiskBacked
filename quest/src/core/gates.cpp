/*
#include "quest/include/gatescheduler.h"
#include "types.h"
#define _USE_MATH_DEFINES
#include <cmath>
#include <random>
#include <complex>
#include <vector>
#include <iostream>
#include <functional>
#include <omp.h>
#include "gates.h"

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

// Quantum Fourier Transform (QFT)
void addQFT(GateScheduler& sched, int numQubits) {
    for (int target = 0; target < numQubits; ++target) {
        addHadamard(sched, target);
        for (int control = target + 1; control < numQubits; ++control) {
            int k = control - target + 1;
            addControlledRK(sched, control, target, k);
        }
    }

    // Optionally reverse qubit order with SWAPs
    //for (int i = 0; i < numQubits / 2; ++i) {
       //addSwap(sched, i, numQubits - i - 1);
    //}
}

// Quantum Supremacy Circuit (QSC)
void addQSC(GateScheduler& sched, int numQubits, int depth) {
    std::mt19937 rng(42);  // Seeded RNG for reproducibility
    std::uniform_int_distribution<int> gateDist(0, 4);

    // Apply initial Hadamard layer
    for (int q = 0; q < numQubits; ++q)
        addHadamard(sched, q);

    // Middle random layers
    for (int d = 0; d < depth; ++d) {
        // 1-qubit layer with random gates
        for (int q = 0; q < numQubits; ++q) {
            switch (gateDist(rng)) {
                case 0: addHadamard(sched, q); break;
                case 1: addT(sched, q); break;
                case 2: addSqrtX(sched, q); break;
                case 3: addSqrtY(sched, q); break;
                case 4: addPhase(sched, q, M_PI / 8.0); break;
            }
        }

        // 2-qubit entangling layer: CZ between nearest neighbors
        for (int q = 0; q < numQubits - 1; q += 2) {
            addCZ(sched, q, q + 1);
        }
        for (int q = 1; q < numQubits - 1; q += 2) {
            addCZ(sched, q, q + 1);
        }
    }

    // Apply final Hadamard layer
    for (int q = 0; q < numQubits; ++q)
        addHadamard(sched, q);
}


// Hadamard gate
void addHadamard(GateScheduler& sched, int qubit) {
    sched.addSingleQubitGate(qubit, [](qcomp& a0, qcomp& a1) {
        const qreal invSqrt2 = 1.0 / std::sqrt(2.0);
        qcomp new0 = invSqrt2 * (a0 + a1);
        qcomp new1 = invSqrt2 * (a0 - a1);
        a0 = new0;
        a1 = new1;
    });
}

// Phase gate: applies exp(i * theta) to |1⟩
void addPhase(GateScheduler& sched, int qubit, qreal theta) {
    sched.addSingleQubitGate(qubit, [theta](qcomp& a0, qcomp& a1) {
        a1 *= std::exp(qcomp(0, theta));
    });
}

// T gate is a π/4 phase gate
void addT(GateScheduler& sched, int qubit) {
    addPhase(sched, qubit, M_PI / 4.0);
}

// Sqrt-X gate
void addSqrtX(GateScheduler& sched, int qubit) {
    sched.addSingleQubitGate(qubit, [](qcomp& a0, qcomp& a1) {
        const qreal inv2 = 0.5;
        qcomp i(0, 1);
        qcomp new0 = inv2 * ((1.0 + i) * a0 + (1.0 - i) * a1);
        qcomp new1 = inv2 * ((1.0 - i) * a0 + (1.0 + i) * a1);
        a0 = new0;
        a1 = new1;
    });
}

// Sqrt-Y gate
void addSqrtY(GateScheduler& sched, int qubit) {
    sched.addSingleQubitGate(qubit, [](qcomp& a0, qcomp& a1) {
        const qreal inv2 = 0.5;
        qcomp i(0, 1);
        qcomp new0 = inv2 * ((1.0 + i) * a0 + (-1.0 - i) * a1);
        qcomp new1 = inv2 * ((1.0 + i) * a0 + (1.0 + i) * a1);
        a0 = new0;
        a1 = new1;
    });
}

// Controlled-phase gate: adds exp(i*theta) if control=1 and target=1
void addControlledPhase(GateScheduler& sched, int control, int target, qreal theta) {
    sched.addControlledGate(control, target, [theta](qcomp& a0, qcomp& a1) {
        a1 *= std::exp(qcomp(0, theta));
    });
}

// Controlled-Z is just controlled phase with theta = pi
void addCZ(GateScheduler& sched, int control, int target) {
    addControlledPhase(sched, control, target, M_PI);
}

void addControlledRK(GateScheduler& scheduler, int control, int target, int k) {
    const double angle = 2.0 * M_PI / (1ULL << k);  // θ = 2π / 2^k
    const qcomp phase = std::polar(1.0, angle);     // exp(iθ)

    scheduler.addControlledGate(control, target, [=](qcomp& a0, qcomp& a1) {
        // Controlled gate: a1 is the |11⟩ amplitude
        a1 *= phase;
    });
}

void addSwap(GateScheduler& scheduler, int qubitA, int qubitB) {
    if (qubitA == qubitB) return;

    scheduler.addSwapGate(qubitA, qubitB, [=](std::vector<qcomp>& buffer, int ) {
        std::cout << "[Debug] SWAP gate called with qubits " << qubitA << ", " << qubitB << "\n";
        std::cout << "[Debug] Buffer size in SWAP: " << buffer.size() << "\n";
        
        const qindex mask = (1ULL << qubitA) | (1ULL << qubitB);
        const qindex size = buffer.size();
        
        std::cout << "[Debug] SWAP mask: " << mask << ", size: " << size << "\n";
        
        std::atomic<bool> hasPrinted = false;
        std::atomic<qindex> first_i(-1);
        std::atomic<qindex> first_j(-1);

        #pragma omp parallel for schedule(static)
        for (qindex i = 0; i < size; ++i) {
            bool bitA = (i >> qubitA) & 1;
            bool bitB = (i >> qubitB) & 1;

            if (bitA != bitB) {
                qindex j = i ^ mask;
                if (i < j) {
                    bool expected = false;
                    if (hasPrinted.compare_exchange_strong(expected, true)) {
                        first_i = i;
                        first_j = j;
                    }

                    std::swap(buffer[i], buffer[j]);
                }
            }
        }

        if (hasPrinted.load()) {
            qindex i = first_i.load();
            qindex j = first_j.load();

            std::cout << "  [SWAP DEBUG] First swap on qubits " << qubitA << " and " << qubitB << "\n";
            std::cout << "    Indices: " << i << " <-> " << j << "\n";
            std::cout << "    After:  buffer[" << i << "] = " << buffer[i]
                      << ", buffer[" << j << "] = " << buffer[j] << "\n";
        }
    });
}

// This is the same function as the above, but without the debug prints
/*
void addSwap(GateScheduler& scheduler, int qubitA, int qubitB) {
    if (qubitA == qubitB) return;

    scheduler.addSwapGate(qubitA, qubitB, [=](std::vector<qcomp>& buffer, int ) {
        const qindex mask = (1ULL << qubitA) | (1ULL << qubitB);
        const qindex size = buffer.size();

        #pragma omp parallel for schedule(static)
        for (qindex i = 0; i < size; ++i) {
            bool bitA = (i >> qubitA) & 1;
            bool bitB = (i >> qubitB) & 1;

            if (bitA != bitB) {
                qindex j = i ^ mask;
                if (i < j) std::swap(buffer[i], buffer[j]);
            }
        }
    });
}

*/