#include "quest/include/hook.h"
#include "quest/include/environment.h"
#include <vector>

// Given a buffer and metadata, create a Qureg that points to it
Qureg createTempQureg(std::vector<qcomp>& buffer, int qubits) {
    QuESTEnv env = getQuESTEnv();
    
    Qureg tempQureg = {};
    tempQureg.cpuAmps = buffer.data();
    tempQureg.numQubits = qubits;
    tempQureg.numAmps = buffer.size();
    tempQureg.logNumAmps = log2(buffer.size());
    tempQureg.numAmpsPerNode = buffer.size();
    tempQureg.logNumAmpsPerNode = qubits;
    tempQureg.isDensityMatrix = 0;
    tempQureg.isDistributed = env.isDistributed;
    tempQureg.isGpuAccelerated = env.isGpuAccelerated;
    tempQureg.isMultithreaded = env.isMultithreaded;
    tempQureg.numNodes = 1;
    return tempQureg;
}

