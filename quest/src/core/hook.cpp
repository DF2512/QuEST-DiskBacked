#include "quest/include/hook.h"
#include "quest/include/environment.h"
#include <vector>

// Given a buffer and metadata, create a Qureg that points to it
Qureg createTempQureg(void* rawBuf, int qubits) {
    QuESTEnv env = getQuESTEnv();

    Qureg tempQureg = {};
    tempQureg.cpuAmps = static_cast<qcomp*>(rawBuf);
    tempQureg.numQubits = qubits;
    tempQureg.numAmps = 1ULL << qubits;
    tempQureg.logNumAmps = qubits;
    tempQureg.numAmpsPerNode = tempQureg.numAmps;
    tempQureg.logNumAmpsPerNode = qubits;
    tempQureg.isDensityMatrix = 0;
    tempQureg.isDistributed = env.isDistributed;
    tempQureg.isGpuAccelerated = env.isGpuAccelerated;
    tempQureg.isMultithreaded = env.isMultithreaded;
    tempQureg.numNodes = 1;
    return tempQureg;
}

