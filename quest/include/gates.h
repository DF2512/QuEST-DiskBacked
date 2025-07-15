#pragma once

#include "gatescheduler.h"
#include "types.h"

void addQFT(GateScheduler& sched, int numQubits);
void addQSC(GateScheduler& sched, int numQubits, int depth);
void addHadamard(GateScheduler&, int);
void addPhase(GateScheduler&, int, qreal);
void addT(GateScheduler&, int);
void addSqrtX(GateScheduler&, int);
void addSqrtY(GateScheduler&, int);
void addControlledPhase(GateScheduler&, int, int, qreal);
void addCZ(GateScheduler&, int, int);
void addSwap(GateScheduler&, int, int);
void addControlledRK(GateScheduler&, int control, int target, int k);
