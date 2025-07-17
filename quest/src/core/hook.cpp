#include "quest/include/hook.h"
#include <vector>

// Forward declarations for original QuEST routines


// Given a buffer and metadata, create a Qureg struct that points to it
Qureg createTempQureg(std::vector<qcomp>& buffer, int qubits) {
    Qureg tempQureg = {};
    tempQureg.cpuAmps = buffer.data();
    tempQureg.numQubits = qubits;
    tempQureg.numAmpsPerNode = buffer.size();
    tempQureg.isDensityMatrix = 0;
    tempQureg.isDistributed = 0;
    tempQureg.isGpuAccelerated = 0;
    tempQureg.isMultithreaded = 0;
    // ...set other fields as needed
    return tempQureg;
}

