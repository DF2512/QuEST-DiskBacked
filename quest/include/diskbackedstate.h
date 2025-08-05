#pragma once

#include <string>
#include <vector>
#include <complex>
#include <liburing.h>
#include "types.h"
#include "chunkmanager.h"
#include "calculations.h"

class DiskBackedState {
public:
    DiskBackedState(int numQubits, int numBlocks, int chunksPerBlock,
                    const std::vector<std::string>& diskRoots, int maxBlocksInMemory);

    int getNumQubits() const;
    int getNumBlocks() const;
    int getNumQubitsPerBlock() const;
    int getNumQubitsPerChunk() const;
    int getMaxPermutableQubits() const;
    int getChunksPerBlock() const;
    int getMaxBlocksInMemory() const;
    size_t getNumAmplitudes() const;
    size_t getNumChunks() const;
    size_t getAmpsPerChunk() const;
    void loadChunk(size_t chunkIndex, void* alignedBuf, std::vector<qcomp>& out) const;
    void saveChunk(size_t chunkIndex, void* alignedBuf, const std::vector<qcomp>& in) const;
    void loadBlock(int blockIdx, const std::vector<int>& chunkIndices, void* alignedBuf, std::vector<qcomp>& buffer) const;
    void saveBlock(int blockIdx, const std::vector<int>& chunkIndices, void* alignedBuf, const std::vector<qcomp>& buffer) const;
    void diskBacked_initRandomPureState();
    void deleteAllChunkFiles();
    void* getAlignedBuffer(int idx) const;
    double computeTotalProbability() const;
    qreal diskBacked_calcTotalProbability() const;
    int diskBacked_applyQubitMeasurement(int qubit);
    void diskBacked_initZeroState();
    ~DiskBackedState();
    PermutationTracker& getPermutationTracker(); 

private:
    int numQubits; // total number of qubits in the state vector
    int numBlocks; // number of blocks
    int chunksPerBlock; // number of chunks per block
    int numQubitsPerChunk; // number of qubits per chunk
    int numQubitsPerBlock; // number of qubits per block
    int maxPermutableQubits; // maximum number of non-local qubits that can be swapped in a single permutation
    size_t numAmplitudes; // total number of amplitudes in the state vector
    size_t numChunks; // total number of chunks
    size_t ampsPerChunk; // number of amplitudes per chunk
    mutable io_uring ring;
    mutable bool ioInitialised = false;
    mutable std::vector<void*> alignedBufferPool;
    mutable std::vector<int> chunkFDs;
    mutable bool filesRegistered = false;
    void ensureIoUringInitialised() const;
    const int maxBlocksInMemory; // maximum number of blocks that can be in memory at once

    std::vector<std::string> diskRoots;
    std::vector<std::string> chunkPaths;

    PermutationTracker permTracker; // tracks current permutation and chunk map

    void generateChunkPaths(); // file initialisation
};