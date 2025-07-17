/** @file
 * A minimum C++ example of running
 * QuEST with disk-backed state management.
 * 
 * @author Tyson Jones
*/

#include "quest.h"
#include "diskbackedstate.h"
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

    double totalProb = state.diskBacked_calcTotalProbability();

    std::cout << "Total probability: " << totalProb << "\n";

    finalizeQuESTEnv();

    return 0;
}

