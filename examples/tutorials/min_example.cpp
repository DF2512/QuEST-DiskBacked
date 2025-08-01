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
#include <chrono>
#include <fstream>
#include <filesystem>

struct RunData {
    int runIdx;
    int numQubits;
    int numBlocks;
    int chunksPerBlock;
    double elapsed;
    double measureElapsed;
    double totalElapsed;
    double initialProb;
    double finalProb;
    bool probCheckPassed;
    int measurementOutcome;
};

int main() {
    initQuESTEnv();
    reportQuESTEnv();

    // Parameters
    std::vector<int> qubitSizes = {26,27,28,29};
    std::vector<int> numBlocksList = {4,8,16,32};
    std::vector<int> chunksPerBlockList = {4,8,16,32};
    const int numRuns = 1; 
    std::vector<std::string> diskRoots = {
        "C:/quantum_chunks0",
        "D:/quantum_chunks1"
    };
    
    std::vector<RunData> runData;
    
    // Find available log file name
    std::string logFileName;
    int logIndex = 0;
    do {
        logFileName = "logs" + std::to_string(logIndex) + ".txt";
        logIndex++;
    } while (std::filesystem::exists(logFileName));
    
    // Initialize log file with header
    std::ofstream logFile(logFileName);
    if (logFile.is_open()) {
        logFile << "QuEST Disk-Backed State Performance Log" << std::endl;
        logFile << "======================================" << std::endl;
        logFile << "Parameters: numRuns=" << numRuns << " per qubit size" << std::endl;
        logFile << "Qubit sizes tested: ";
        for (size_t i = 0; i < qubitSizes.size(); ++i) {
            logFile << qubitSizes[i];
            if (i < qubitSizes.size() - 1) logFile << ", ";
        }
        logFile << std::endl;
        logFile << std::endl;
        
        logFile << "Qubits | Blocks | Chunks | Run | Elapsed(s) | Measure(s) | Total(s) | InitialProb | FinalProb | ProbCheck | Outcome" << std::endl;
        logFile << "-------|--------|--------|-----|------------|------------|----------|-------------|-----------|-----------|---------" << std::endl;
        logFile.flush(); // Ensure header is written immediately
    } else {
        std::cerr << "Failed to open log file: " << logFileName << std::endl;
        return 1;
    }
    
    std::cout << "Running " << numRuns << " iterations for each qubit size..." << std::endl;
    std::cout << "Results will be saved to: " << logFileName << std::endl;
    
    for (size_t qubitIndex = 0; qubitIndex < qubitSizes.size(); ++qubitIndex) {
        int numQubits = qubitSizes[qubitIndex];
        int numBlocks = numBlocksList[qubitIndex];
        int chunksPerBlock = chunksPerBlockList[qubitIndex];
        
        std::cout << "\n==========================================" << std::endl;
        std::cout << "Testing " << numQubits << " qubits (numBlocks=" << numBlocks 
                  << ", chunksPerBlock=" << chunksPerBlock << ")" << std::endl;
        std::cout << "==========================================" << std::endl;
        
        for (int run = 1; run <= numRuns; ++run) {
            std::cout << "Run " << run << "/" << numRuns << " for " << numQubits << " qubits" << std::endl;
            
            RunData data;
            data.runIdx = run;
            data.numQubits = numQubits;
            data.numBlocks = numBlocks;
            data.chunksPerBlock = chunksPerBlock;
            
            // Begin timer
            auto start = std::chrono::high_resolution_clock::now();

            // Create disk-backed state, initialise with random pure state
            DiskBackedState diskState(numQubits, numBlocks, chunksPerBlock, diskRoots);
            diskState.diskBacked_initRandomPureState();
            //Qureg qureg = createForcedQureg(numQubits);
            //initRandomPureState(qureg);

            // Create a schedule, add either a QFT or a QSC with chosen parameters
            GateScheduler scheduler;
            scheduler.addFullQFT(numQubits);
            // scheduler.addQSC(numQubits, 8); // Adjust depth as needed
            //applyFullQuantumFourierTransform(qureg);

            // Apply the schedule to the disk-backed state
            runCircuit(scheduler, diskState, false);

            // Report time to initialise and run the circuit
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end - start;
            data.elapsed = elapsed.count();

            // Check total probability, start new timer
            auto measureStart = std::chrono::high_resolution_clock::now();

            qreal totalProb = diskState.diskBacked_calcTotalProbability();
            //qreal totalProb = calcTotalProb(qureg);
            data.initialProb = totalProb;

            // Measure a qubit
            int outcome = diskState.diskBacked_applyQubitMeasurement(0); // Choose a qubit to measure
            //int outcome = applyQubitMeasurement(qureg, 0); 
            data.measurementOutcome = outcome;

            // Check probability again to ensure it remains 1
            qreal finalProb = diskState.diskBacked_calcTotalProbability();
            //qreal finalProb = calcTotalProb(qureg);
            data.finalProb = finalProb;
            data.probCheckPassed = (std::abs(finalProb - 1.0) < 1e-10);

            // Check time taken for measurement and probability checks
            auto measureEnd = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> measureElapsed = measureEnd - measureStart;
            data.measureElapsed = measureElapsed.count();

            // Check total runtime
            std::chrono::duration<double> totalElapsed = measureEnd - start;
            data.totalElapsed = totalElapsed.count();
            
            runData.push_back(data);
            
            std::cout << "  Elapsed: " << data.elapsed << "s, Measure: " << data.measureElapsed 
                      << "s, Total: " << data.totalElapsed << "s, Prob: " << data.finalProb 
                      << ", Outcome: " << data.measurementOutcome << std::endl;
            
            // Write this run's data to log file immediately
            logFile << std::setw(6) << data.numQubits << " | "
                    << std::setw(6) << data.numBlocks << " | "
                    << std::setw(6) << data.chunksPerBlock << " | "
                    << std::setw(3) << data.runIdx << " | "
                    << std::fixed << std::setprecision(6) << std::setw(10) << data.elapsed << " | "
                    << std::setw(10) << data.measureElapsed << " | "
                    << std::setw(8) << data.totalElapsed << " | "
                    << std::setw(11) << data.initialProb << " | "
                    << std::setw(9) << data.finalProb << " | "
                    << std::setw(9) << (data.probCheckPassed ? "PASS" : "FAIL") << " | "
                    << std::setw(7) << data.measurementOutcome << std::endl;
            logFile.flush(); // Ensure data is written immediately
        }
    }
    
    // Calculate and write summary statistics
    logFile << std::endl;
    logFile << "Summary Statistics by Qubit Size:" << std::endl;
    logFile << "=================================" << std::endl;
    
    for (int qubits : qubitSizes) {
        double avgElapsed = 0.0, avgMeasure = 0.0, avgTotal = 0.0;
        int probPassCount = 0;
        int count = 0;
        
        for (const auto& data : runData) {
            if (data.numQubits == qubits) {
                avgElapsed += data.elapsed;
                avgMeasure += data.measureElapsed;
                avgTotal += data.totalElapsed;
                if (data.probCheckPassed) probPassCount++;
                count++;
            }
        }
        
        if (count > 0) {
            avgElapsed /= count;
            avgMeasure /= count;
            avgTotal /= count;
            
            logFile << "Qubits " << qubits << ":" << std::endl;
            logFile << "  Average Elapsed Time: " << std::fixed << std::setprecision(6) << avgElapsed << "s" << std::endl;
            logFile << "  Average Measure Time: " << avgMeasure << "s" << std::endl;
            logFile << "  Average Total Time: " << avgTotal << "s" << std::endl;
            logFile << "  Probability Check Pass Rate: " << probPassCount << "/" << count 
                    << " (" << (100.0 * probPassCount / count) << "%)" << std::endl;
            logFile << std::endl;
        }
    }
    
    logFile.close();
    std::cout << "Final results saved to: " << logFileName << std::endl;
    
    finalizeQuESTEnv();
    return 0;
   
}