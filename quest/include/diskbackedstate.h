#pragma once

#include <vector>
#include <string>
#include <filesystem>
#include <atomic>
#include <thread>
#include <complex>
#include <functional>
#include <liburing.h>
#include <mutex>
#include <condition_variable>
#include <deque>


#include "types.h"
#include "chunkmanager.h"
#include "calculations.h"

class DiskBackedState {
public:
    DiskBackedState(int numQubits, int numBlocks, int chunksPerBlock,
                    const std::vector<std::string>& diskRoots, int maxBlocksInMemory);

    // --- State info accessors ---
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

    // --- Buffer helpers ---
    int findBufferIndex(void* buf) const;
    void* getAlignedBuffer(int blockBufIndex) const;
    void* getChunkBuffer(int blockBufIndex, int localChunkIndex) const;

    // --- Disk I/O (single chunk, blocking) ---
    void loadChunk(size_t chunkIndex, void* alignedBuf) const;
    void saveChunk(size_t chunkIndex, void* alignedBuf) const;

    // --- Disk I/O (entire block, async) ---
    void loadBlock(int blockIdx,
                   const std::vector<int>& chunkIndices,
                   std::function<void()> onComplete) const;

    void saveBlock(int blockIdx,
                   const std::vector<int>& chunkIndices,
                   std::function<void()> onComplete) const;

    // --- State initialisation ---
    void diskBacked_initRandomPureState();
    void diskBacked_initZeroState();
    void diskBacked_initPlusState();
    void diskBacked_initBlankState();

    // --- State analysis ---
    double computeTotalProbability() const;
    qreal diskBacked_calcTotalProbability() const;
    int diskBacked_applyQubitMeasurement(int qubit);

    // --- Management ---
    void ioCompletionLoop() const;
    void stopIoWorker() const;
    void deleteAllChunkFiles();
    PermutationTracker& getPermutationTracker();

    ~DiskBackedState();

private:
    // --- Quantum layout ---
    int numQubits;                // total number of qubits in the state vector
    int numBlocks;                 // number of blocks
    int chunksPerBlock;            // number of chunks per block
    int numQubitsPerChunk;         // number of qubits per chunk
    int numQubitsPerBlock;         // number of qubits per block
    int maxPermutableQubits;       // maximum number of non-local qubits that can be swapped
    size_t numAmplitudes;          // total amplitudes in the state vector
    size_t numChunks;              // total number of chunks
    size_t ampsPerChunk;           // number of amplitudes per chunk

    // --- I/O state ---
    const int maxBlocksInMemory;   // max number of blocks in memory at once
    std::vector<int> chunkFDs;     // file descriptors for each chunk file
    std::vector<std::string> diskRoots;
    std::vector<std::string> chunkPaths;
    std::vector<struct iovec> registeredIovecs;

    mutable io_uring ring;
    mutable std::atomic<bool> ioInitialised = false;
    mutable std::atomic<bool> filesRegistered = false;
    mutable std::atomic<int> totalSubmitted{0};
    mutable std::atomic<bool> ioRunning = false;
    mutable std::thread ioCompletionThread;

    mutable std::mutex completionQueueMtx;
    mutable std::condition_variable completionQueueCv;
    mutable std::deque<std::function<void()>> completionWorkQueue;
    mutable std::thread completionWorker;
    mutable std::atomic<bool> completionWorkerRunning{false};

    // --- Buffers ---
    mutable std::vector<void*> blockBuffers;                     // base pointer for each in-memory block
    mutable std::vector<std::vector<void*>> chunkViews;           // per-block chunk pointers

    // --- Permutation tracking ---
    PermutationTracker permTracker;

    // --- Internal helpers ---
    void ensureIoUringInitialised() const;
    void generateChunkPaths();
};

// Context used for async block I/O
struct BlockIOContext {
    int blockIdx;
    std::vector<int> chunkIndices;
    void* userData; // placeholder for future, unused in current .cpp
    std::atomic<int> remainingChunks;
    std::function<void()> onComplete;
};
