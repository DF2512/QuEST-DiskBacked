#pragma once

#include <string>
#include <vector>
#include <complex>
#include <types.h>
#include "chunkmanager.h"


using qreal = double;
using qcomp = std::complex<qreal>;

class DiskBackedState {
public:
    DiskBackedState(int numQubits, int numBlocks, int chunksPerBlock,
                    const std::vector<std::string>& diskRoots);

    int getNumQubits() const;
    int getNumBlocks() const;
    int getNumQubitsPerBlock() const;
    int getMaxPermutableQubits() const;
    int getChunksPerBlock() const;
    size_t getNumAmplitudes() const;
    size_t getNumChunks() const;
    size_t getAmpsPerChunk() const;
    void loadChunk(size_t chunkIndex, std::vector<qcomp>& buffer) const;
    void saveChunk(size_t chunkIndex, const std::vector<qcomp>& buffer) const;
    void loadBlock(int blockIdx, const std::vector<int>& chunkIndices, std::vector<qcomp>& buffer) const;
    void saveBlock(int blockIdx, const std::vector<int>& chunkIndices, const std::vector<qcomp>& buffer) const;
    void initialiseRandomState();
    void deleteAllChunkFiles();
    double computeTotalProbability() const;
    ~DiskBackedState();

    PermutationTracker& getPermutationTracker(); // NEW
    //const PermutationTracker& getPermutationTracker() const;
    

private:
    int numQubits;
    int numBlocks;
    int chunksPerBlock;
    int numQubitsPerChunk;
    int numQubitsPerBlock;
    int maxPermutableQubits;
    size_t numAmplitudes;
    size_t numChunks;
    size_t ampsPerChunk;

    std::vector<std::string> diskRoots;
    std::vector<std::string> chunkPaths;

    PermutationTracker permTracker; // tracks current permutation and chunk map

    void generateChunkPaths(); // internal helper
};