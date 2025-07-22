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
#include <vector>
#include <string>
#include <iostream>

int main(void) {

    initQuESTEnv();
    reportQuESTEnv();

    std::vector<std::string> diskRoots = {
        "C:/quantum_chunks0",
        "D:/quantum_chunks1"
    };
    
    DiskBackedState state(28, 8, 8, diskRoots);
    state.diskBacked_initRandomPureState();

    qreal totalProb = state.diskBacked_calcTotalProbability();
    reportScalar("Total Probability", totalProb);

    GateScheduler scheduler;
    for (int i = 0; i <= 27; ++i) {
        scheduler.addHadamard(i);
    }
    
    runCircuit(scheduler, state, true);

    qreal finalProb = state.diskBacked_calcTotalProbability();
    reportScalar("Final Probability", finalProb);
    finalizeQuESTEnv();


    return 0;
}

