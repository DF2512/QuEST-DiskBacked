#include "diskbackedstate.h"
#include "quest/include/qureg.h"
#include "quest/include/initialisations.h"
#include "quest/include/hook.h"
#include "quest/include/types.h"
#include <filesystem>
#include <sys/resource.h>
#include <iostream>
#include <memory>
#include <cmath>
#include <stdexcept>
#include <random>
#include <cstdio>
#include <unistd.h>
#include <cstring>
#include <omp.h>
#include "quest/src/core/randomiser.hpp"
#include "quest/src/core/localiser.hpp"

// Constructor for DiskBackedState
DiskBackedState::DiskBackedState(int numQubits_, int numBlocks_, int chunksPerBlock_,
                                 const std::vector<std::string>& diskRoots_, int maxBlocksInMemory_)
    : numQubits(numQubits_),
      numBlocks(numBlocks_),
      chunksPerBlock(chunksPerBlock_),
      diskRoots(diskRoots_),
      maxBlocksInMemory(maxBlocksInMemory_),
      permTracker(numQubits_, numBlocks_, chunksPerBlock_)
{
    numAmplitudes = 1ULL << numQubits;
    numChunks = numBlocks * chunksPerBlock;
    ampsPerChunk = numAmplitudes / numChunks;
    numQubitsPerChunk = log2(ampsPerChunk);
    numQubitsPerBlock = log2(chunksPerBlock * ampsPerChunk);
    maxPermutableQubits = log2(chunksPerBlock);

    if (numAmplitudes % numChunks != 0) {
        throw std::runtime_error("Amplitudes must divide evenly across chunks");
    }
    struct io_uring_params params;
    memset(&params, 0, sizeof(params));
    params.flags |= IORING_SETUP_SQPOLL;
    params.sq_thread_idle = 2000;

    generateChunkPaths();

    chunkFDs.resize(numChunks);
    for (size_t i = 0; i < numChunks; ++i) {
        int fd = open(chunkPaths[i].c_str(), O_RDWR | O_CREAT | O_DIRECT, 0644);
        if (fd < 0) {
            throw std::runtime_error(
                "Failed to open file: " + chunkPaths[i] +
                " errno=" + std::to_string(errno) + " (" + strerror(errno) + ")"
            );
        }
        chunkFDs[i] = fd;
    }

    if (io_uring_queue_init_params((chunksPerBlock * maxBlocksInMemory * 2) + 1, &ring, &params) < 0) {
        throw std::runtime_error("Failed to initialize io_uring");
    }

    const size_t chunkBytes = ampsPerChunk * sizeof(qcomp);
    const size_t blockBytes = chunkBytes * chunksPerBlock;
    const size_t totalBuffers = maxBlocksInMemory * chunksPerBlock;

    registeredIovecs.resize(totalBuffers);
    blockBuffers.resize(maxBlocksInMemory);
    chunkViews.resize(maxBlocksInMemory);

    for (int blockIdx = 0; blockIdx < maxBlocksInMemory; ++blockIdx) {
        void* blockBuf = nullptr;
        if (posix_memalign(&blockBuf, 4096, blockBytes) != 0) {
            throw std::runtime_error("Block buffer allocation failed");
        }

        blockBuffers[blockIdx] = blockBuf;
        chunkViews[blockIdx].resize(chunksPerBlock);

        for (int chunkIdx = 0; chunkIdx < chunksPerBlock; ++chunkIdx) {
            void* chunkPtr = static_cast<char*>(blockBuf) + chunkIdx * chunkBytes;
            size_t bufIndex = blockIdx * chunksPerBlock + chunkIdx;
            registeredIovecs[bufIndex].iov_base = chunkPtr;
            registeredIovecs[bufIndex].iov_len = chunkBytes;
            chunkViews[blockIdx][chunkIdx] = chunkPtr;
        }
    }

    if (io_uring_register_buffers(&ring, registeredIovecs.data(),
                                  static_cast<unsigned int>(registeredIovecs.size())) < 0) {
        int err = errno;
        throw std::runtime_error("io_uring_register_buffers failed");
    }
    ioInitialised = true;
}

// Buffer accessor
void* DiskBackedState::getAlignedBuffer(int blockBufIndex) const {
    return blockBuffers.at(blockBufIndex);
}

// Registered chunk buffer accessor
void* DiskBackedState::getChunkBuffer(int blockBufIndex, int localChunkIndex) const {
    return chunkViews.at(blockBufIndex).at(localChunkIndex);
}

// File generation helper
void DiskBackedState::generateChunkPaths() {
    chunkPaths.resize(numChunks);
    for (size_t i = 0; i < numChunks; ++i) {
        size_t diskIndex = i % diskRoots.size();
        std::filesystem::create_directories(diskRoots[diskIndex]);
        chunkPaths[i] = diskRoots[diskIndex] + "/chunk_" + std::to_string(i) + ".dat";
    }
}

PermutationTracker& DiskBackedState::getPermutationTracker() {
    return permTracker;
}


// Accessors
int DiskBackedState::getChunksPerBlock() const { return chunksPerBlock; }
int DiskBackedState::getMaxPermutableQubits() const {return maxPermutableQubits; }
int DiskBackedState::getNumQubitsPerBlock() const { return numQubitsPerBlock; }
int DiskBackedState::getNumQubitsPerChunk() const { return numQubitsPerChunk; }
int DiskBackedState::getNumBlocks() const { return numBlocks; }
int DiskBackedState::getNumQubits() const { return numQubits; }
int DiskBackedState::getMaxBlocksInMemory() const { return maxBlocksInMemory; }
size_t DiskBackedState::getNumAmplitudes() const { return numAmplitudes; }
size_t DiskBackedState::getNumChunks() const { return numChunks; }
size_t DiskBackedState::getAmpsPerChunk() const { return ampsPerChunk; }

void DiskBackedState::ensureIoUringInitialised() const {
    if (!ioInitialised.load()) {
        bool expected = false;
        if (ioInitialised.compare_exchange_strong(expected, true)) {
            if (io_uring_queue_init((chunksPerBlock * maxBlocksInMemory * 2) + 1, &ring, 0) < 0) {
                ioInitialised.store(false);
                int err = errno;
                throw std::runtime_error("Failed to initialize io_uring");
            }

            if (io_uring_register_buffers(&ring, registeredIovecs.data(),
                                          static_cast<unsigned int>(registeredIovecs.size())) < 0) {
                int err = errno;
                io_uring_queue_exit(&ring);
                ioInitialised.store(false);
                throw std::runtime_error("io_uring_register_buffers failed");
            }

            ioInitialised.store(true);
        } else {
            while (!ioInitialised.load()) {
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
            }
        }
    }

    if (!filesRegistered.load()) {
        bool expectedFiles = false;
        if (filesRegistered.compare_exchange_strong(expectedFiles, true)) {
            if (io_uring_register_files(&ring, chunkFDs.data(), chunkFDs.size()) < 0) {
                int err = errno;
                filesRegistered.store(false);
                throw std::runtime_error("io_uring_register_files failed");
            }
        } else {
            while (!filesRegistered.load()) {
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
            }
        }
    }

    bool expectedRunning = false;
    if (ioRunning.compare_exchange_strong(expectedRunning, true)) {
        ioCompletionThread = std::thread([this]() {
            this->ioCompletionLoop();
        });
    }

    if (!completionWorkerRunning.load()) {
        bool expectedWorker = false;
        if (completionWorkerRunning.compare_exchange_strong(expectedWorker, true)) {
            completionWorker = std::thread([this]() {
                while (completionWorkerRunning.load()) {
                    std::function<void()> task;
                    {
                        std::unique_lock<std::mutex> lk(this->completionQueueMtx);
                        this->completionQueueCv.wait(lk, [this]() {
                            return !this->completionWorkQueue.empty() || !completionWorkerRunning.load();
                        });
                        if (!completionWorkerRunning.load() && this->completionWorkQueue.empty())
                            break;
                        if (!this->completionWorkQueue.empty()) {
                            task = std::move(this->completionWorkQueue.front());
                            this->completionWorkQueue.pop_front();
                        }
                    }
                    if (task) {
                        try {
                            task();
                        } catch (const std::exception& e) {
                        } catch (...) {
                        }
                    }
                }
            });
        } else {
            while (!completionWorkerRunning.load()) {
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
            }
        }
    }
}

int DiskBackedState::findBufferIndex(void* buf) const {
    for (size_t i = 0; i < registeredIovecs.size(); ++i) {
        if (registeredIovecs[i].iov_base == buf)
            return static_cast<int>(i);
    }
    return -1;
}

void printRingStatus(int submitted,  int ring_size) {
    printf("\r[IO Ring] In queue: %d / %d      ", submitted, ring_size);
    fflush(stdout);
}

// Access a single file (for debugging)
void DiskBackedState::loadChunk(size_t chunkIndex, void* alignedBuf) const {
    ensureIoUringInitialised();

    if (chunkIndex >= numChunks)
        throw std::runtime_error("loadChunk: chunkIndex out of bounds");

    const size_t chunkBytes = ampsPerChunk * sizeof(qcomp);
    int fdIndex = static_cast<int>(chunkIndex);
    int bufIndex = findBufferIndex(alignedBuf);

    if (bufIndex < 0)
        throw std::runtime_error("loadChunk: buffer not registered");

    struct io_uring_sqe* sqe = io_uring_get_sqe(&ring);
    if (!sqe)
        throw std::runtime_error("loadChunk: failed to get SQE");

    io_uring_prep_read_fixed(sqe, fdIndex, alignedBuf, chunkBytes, 0, bufIndex);
    sqe->flags |= IOSQE_FIXED_FILE;
    io_uring_sqe_set_data(sqe, nullptr);

    if (io_uring_submit(&ring) < 0)
        throw std::runtime_error("loadChunk: submit failed");

    struct io_uring_cqe* cqe;
    if (io_uring_wait_cqe(&ring, &cqe) < 0)
        throw std::runtime_error("loadChunk: failed to wait for CQE");

    if (cqe->res < 0 || static_cast<size_t>(cqe->res) != chunkBytes) {
        int err = -cqe->res;
        io_uring_cqe_seen(&ring, cqe);
        throw std::runtime_error("loadChunk: failed IO read (" + std::string(strerror(err)) + ")");
    }

    io_uring_cqe_seen(&ring, cqe);
}

// Save a single file (for debugging)
void DiskBackedState::saveChunk(size_t chunkIndex, void* alignedBuf) const {
    ensureIoUringInitialised();

    if (chunkIndex >= numChunks)
        throw std::runtime_error("saveChunk: chunkIndex out of bounds");

    const size_t chunkBytes = ampsPerChunk * sizeof(qcomp);
    int fdIndex = static_cast<int>(chunkIndex);
    int bufIndex = findBufferIndex(alignedBuf);

    if (bufIndex < 0)
        throw std::runtime_error("saveChunk: buffer not registered");

    struct io_uring_sqe* sqe = io_uring_get_sqe(&ring);
    if (!sqe)
        throw std::runtime_error("saveChunk: failed to get SQE");

    io_uring_prep_write_fixed(sqe, fdIndex, alignedBuf, chunkBytes, 0, bufIndex);
    sqe->flags |= IOSQE_FIXED_FILE;
    io_uring_sqe_set_data(sqe, nullptr);

    if (io_uring_submit(&ring) < 0)
        throw std::runtime_error("saveChunk: submit failed");

    struct io_uring_cqe* cqe;
    if (io_uring_wait_cqe(&ring, &cqe) < 0)
        throw std::runtime_error("saveChunk: failed to wait for CQE");

    if (cqe->res < 0 || static_cast<size_t>(cqe->res) != chunkBytes) {
        int err = -cqe->res;
        io_uring_cqe_seen(&ring, cqe);
        throw std::runtime_error("saveChunk: failed IO write (" + std::string(strerror(err)) + ")");
    }

    io_uring_cqe_seen(&ring, cqe);
}

// Block reading
void DiskBackedState::loadBlock(
    int blockIdx,
    const std::vector<int>& chunkIndices,
    std::function<void()> onComplete
) const {
    if (chunkIndices.size() != static_cast<size_t>(chunksPerBlock))
        throw std::runtime_error("loadBlock: chunkIndices size mismatch");

    ensureIoUringInitialised();
    const size_t chunkBytes = ampsPerChunk * sizeof(qcomp);

    auto* rawCtx = new BlockIOContext{
        blockIdx,
        chunkIndices,
        nullptr,
        0,
        std::move(onComplete)
    };

    auto* sptrHeap = new std::shared_ptr<BlockIOContext>(std::shared_ptr<BlockIOContext>(rawCtx));

    int memSlot = blockIdx % maxBlocksInMemory;
    int submittedChunks = 0;

    for (int i = 0; i < chunksPerBlock; ++i) {
        struct io_uring_sqe* sqe = io_uring_get_sqe(&ring);
        if (!sqe) {
            delete sptrHeap;
            throw std::runtime_error("loadBlock: failed to get SQE");
        }

        int fileIndex = chunkIndices[i];
        int bufIndex  = memSlot * chunksPerBlock + i;

        io_uring_prep_read_fixed(
            sqe,
            fileIndex,
            chunkViews[memSlot][i],
            chunkBytes,
            0,
            bufIndex
        );
        sqe->flags |= IOSQE_FIXED_FILE;

        io_uring_sqe_set_data(sqe, sptrHeap);
        submittedChunks++;
    }

    rawCtx->remainingChunks.store(submittedChunks);
    totalSubmitted.fetch_add(submittedChunks);

    int rc = io_uring_submit(&ring);
    if (rc < 0) {
        totalSubmitted.fetch_sub(submittedChunks);
        delete sptrHeap;
        throw std::runtime_error("loadBlock: submit failed");
    }

    if (rc == 0) {
        totalSubmitted.fetch_sub(submittedChunks);
        delete sptrHeap;
        throw std::runtime_error("loadBlock: submit accepted 0 SQEs");
    }

    if (rc < submittedChunks) {
        int notSubmitted = submittedChunks - rc;
        totalSubmitted.fetch_sub(notSubmitted);
        rawCtx->remainingChunks.store(rc);
    }
}

// Block writing
void DiskBackedState::saveBlock(
    int blockIdx,
    const std::vector<int>& chunkIndices,
    std::function<void()> onComplete
) const {
    if (chunkIndices.size() != static_cast<size_t>(chunksPerBlock))
        throw std::runtime_error("saveBlock: chunkIndices size mismatch");

    ensureIoUringInitialised();
    const size_t chunkBytes = ampsPerChunk * sizeof(qcomp);

    auto* rawCtx = new BlockIOContext{
        blockIdx,
        chunkIndices,
        nullptr,
        0,
        std::move(onComplete)
    };

    auto* sptrHeap = new std::shared_ptr<BlockIOContext>(std::shared_ptr<BlockIOContext>(rawCtx));

    int memSlot = blockIdx % maxBlocksInMemory;
    int submittedChunks = 0;

    for (int i = 0; i < chunksPerBlock; ++i) {
        struct io_uring_sqe* sqe = io_uring_get_sqe(&ring);
        if (!sqe) {
            delete sptrHeap;
            throw std::runtime_error("saveBlock: failed to get SQE");
        }

        int fileIndex = chunkIndices[i];
        int bufIndex  = memSlot * chunksPerBlock + i;

        io_uring_prep_write_fixed(
            sqe,
            fileIndex,
            chunkViews[memSlot][i],
            chunkBytes,
            0,
            bufIndex
        );
        sqe->flags |= IOSQE_FIXED_FILE;

        io_uring_sqe_set_data(sqe, sptrHeap);
        submittedChunks++;
    }

    rawCtx->remainingChunks.store(submittedChunks);
    totalSubmitted.fetch_add(submittedChunks);

    int rc = io_uring_submit(&ring);
    if (rc < 0) {
        totalSubmitted.fetch_sub(submittedChunks);
        delete sptrHeap;
        throw std::runtime_error("saveBlock: submit failed");
    }

    if (rc == 0) {
        totalSubmitted.fetch_sub(submittedChunks);
        delete sptrHeap;
        throw std::runtime_error("saveBlock: submit accepted 0 SQEs");
    }

    if (rc < submittedChunks) {
        int notSubmitted = submittedChunks - rc;
        totalSubmitted.fetch_sub(notSubmitted);
        rawCtx->remainingChunks.store(rc);
    }
}

// io_uring completion queue thread
void DiskBackedState::ioCompletionLoop() const {
    using namespace std::chrono;
    const milliseconds idleSleep{1};
    const milliseconds noDataSleep{2};

    while (ioRunning.load()) {
        try {
            io_uring_cqe* cqe = nullptr;
            int rc = io_uring_peek_cqe(&ring, &cqe);

            if (rc == -EAGAIN) {
                std::this_thread::sleep_for(idleSleep);
                continue;
            } else if (rc < 0) {
                std::this_thread::sleep_for(idleSleep);
                continue;
            }

            if (!cqe) {
                std::this_thread::sleep_for(idleSleep);
                continue;
            }

            void* data = io_uring_cqe_get_data(cqe);
            int res = cqe->res;

            if (!data) {
                std::this_thread::sleep_for(noDataSleep);
                continue;
            }

            // Remove the CQE from the CQ (we own it)
            io_uring_cqe_seen(&ring, cqe);

            // Update global submitted counter
            int prevTotal = totalSubmitted.fetch_sub(1);
            int newTotal = prevTotal - 1;
            if (newTotal < 0) {
                continue;
            }

            // Interpret userdata as pointer-to-shared_ptr<BlockIOContext>
            auto* sptrHeap = static_cast<std::shared_ptr<BlockIOContext>*>(data);
            if (!sptrHeap) {
                continue;
            }

            std::shared_ptr<BlockIOContext> ctx_sp = *sptrHeap;
            BlockIOContext* ctx = ctx_sp.get();
            if (!ctx) {
                continue;
            }

            // Sanity-checks
            int ctxBlock = ctx->blockIdx;
            size_t ctxChunkCount = ctx->chunkIndices.size();
            if (ctxBlock < 0 || ctxBlock >= numBlocks) {
                continue;
            }
            if (ctxChunkCount == 0 || ctxChunkCount > static_cast<size_t>(chunksPerBlock) * 1000) {
                continue;
            }

            if (res < 0) {
                std::terminate();
            }

            int prev = ctx->remainingChunks.fetch_sub(1);
            int rem = prev - 1;

            if (rem < 0) {
                continue;
            }

            if (rem == 0) {
                std::function<void()> cb;
                try {
                    cb = ctx->onComplete;
                } catch (...) {
                    cb = nullptr;
                }

                if (cb) {
                    {
                        std::lock_guard<std::mutex> lk(this->completionQueueMtx);
                        this->completionWorkQueue.push_back(std::move(cb));
                    }
                    this->completionQueueCv.notify_one();
                }

                delete sptrHeap;
            }

           // printRingStatus((int)totalSubmitted.load(), (chunksPerBlock * maxBlocksInMemory * 2) + 1);

        } catch (const std::exception &e) {
            std::this_thread::sleep_for(milliseconds(10));
        } catch (...) {
            std::this_thread::sleep_for(milliseconds(10));
        }
    }
}


void DiskBackedState::stopIoWorker() const {
    // Stop the completion worker gracefully
    completionWorkerRunning.store(false);
    {
        std::lock_guard<std::mutex> lk(completionQueueMtx);
        // leave the queue as-is; worker will drain remaining items if desired.
    }
    completionQueueCv.notify_all();
    if (completionWorker.joinable()) {
        completionWorker.join();
    }
}

void DiskBackedState::deleteAllChunkFiles() {
    for (const auto& file : chunkPaths) {
        std::remove(file.c_str());
    }
}

DiskBackedState::~DiskBackedState() {
    // Stop accepting new IO
    ioRunning = false;

    if (ioInitialised.load()) {
        // Wake up the completion loop with a no-op SQE
        io_uring_sqe* sqe = io_uring_get_sqe(&ring);
        if (sqe) {
            io_uring_prep_nop(sqe);
        }
        io_uring_submit(&ring);

        // Wait for reaper thread to finish
        if (ioCompletionThread.joinable()) {
            ioCompletionThread.join();
        }

        // Now stop and join the worker thread (drains all enqueued callbacks)
        stopIoWorker();

        // Unregister files if they were registered
        if (filesRegistered.load()) {
            io_uring_unregister_files(&ring);
            filesRegistered.store(false);
        }

        // Unregister buffers
        io_uring_unregister_buffers(&ring);

        // Exit the ring
        io_uring_queue_exit(&ring);
        ioInitialised.store(false);
    }

    // Close all FDs
    for (int fd : chunkFDs) {
        if (fd >= 0) {
            close(fd);
        }
    }
    deleteAllChunkFiles();

    // Free memory buffers
    for (void* blockBuf : blockBuffers) {
        free(blockBuf);
    }
}


void DiskBackedState::diskBacked_initRandomPureState() {
    // --- Preconditions & derived sizes (read-only) ---
    const int numBlocksLocal = getNumBlocks();
    const int qubitsPerBlock = getNumQubitsPerBlock();
    const int maxBlocks = getMaxBlocksInMemory();

    const size_t chunkBytes = ampsPerChunk * sizeof(qcomp);
    const size_t blockAmps = static_cast<size_t>(chunksPerBlock) * ampsPerChunk;
    const size_t blockBytes = blockAmps * sizeof(qcomp);

    const double normFactor = 1.0 / std::sqrt(static_cast<double>(numBlocksLocal));

    if (blockBytes == 0) {
        return;
    }

    // --- Synchronisation primitives (match runCircuit() style) ---
    std::mutex memMtx;
    std::condition_variable memCv;
    int inFlightBlocks = 0;

    std::vector<bool> bufferBusy(maxBlocks, false);

    // Ensure buffers exist
    for (int i = 0; i < maxBlocks; ++i) {
        void* p = getAlignedBuffer(i);
        if (p == nullptr) {
            throw std::runtime_error("Aligned buffer missing");
        }
    }

    // Main loop: process each block
    for (int blockIdx = 0; blockIdx < numBlocksLocal; ++blockIdx) {
        int bufferIndex = blockIdx % maxBlocks;

        // Wait until buffer is free and capacity available
        {
            std::unique_lock<std::mutex> lock(memMtx);
            memCv.wait(lock, [&]() {
                return (inFlightBlocks < maxBlocks) && (!bufferBusy[bufferIndex]);
            });
            inFlightBlocks++;
            bufferBusy[bufferIndex] = true;
        }

        // Zero the buffer
        void* blockBuf = getAlignedBuffer(bufferIndex);
        std::memset(blockBuf, 0, blockBytes);

        // Create a Qureg that uses the block buffer
        Qureg tempQureg = createTempQureg(blockBuf, qubitsPerBlock);
        if (tempQureg.cpuAmps == nullptr) {
            try { destroyQureg(tempQureg); } catch (...) {}
            {
                std::lock_guard<std::mutex> lk(memMtx);
                inFlightBlocks--;
                bufferBusy[bufferIndex] = false;
            }
            memCv.notify_all();
            throw std::runtime_error("createTempQureg produced null cpuAmps");
        }

        // Initialize random pure state
        initRandomPureState(tempQureg);

        // Normalise amplitudes
        for (size_t a = 0; a < blockAmps; ++a) {
            tempQureg.cpuAmps[a] *= normFactor;
        }

        // Build chunkIndices vector
        std::vector<int> chunkIndices;
        chunkIndices.reserve(chunksPerBlock);
        for (int local = 0; local < chunksPerBlock; ++local) {
            chunkIndices.push_back(blockIdx * chunksPerBlock + local);
        }

        // Save block asynchronously
        saveBlock(blockIdx, chunkIndices,
                  [&, bufferIndex, blockIdx]() {
                      {
                          std::lock_guard<std::mutex> lock(memMtx);
                          inFlightBlocks--;
                          bufferBusy[bufferIndex] = false;
                      }
                      memCv.notify_all();
                  });
    }

    // Wait for all in-flight saves to finish
    {
        std::unique_lock<std::mutex> lock(memMtx);
        memCv.wait(lock, [&]() {
            return inFlightBlocks == 0;
        });
    }
}

void DiskBackedState::diskBacked_initZeroState() {
    // --- Preconditions & derived sizes (read-only) ---
    const int numBlocksLocal = getNumBlocks();
    const int qubitsPerBlock = getNumQubitsPerBlock();
    const int maxBlocks = getMaxBlocksInMemory();

    const size_t blockAmps = static_cast<size_t>(chunksPerBlock) * ampsPerChunk;
    const size_t blockBytes = blockAmps * sizeof(qcomp);

    if (blockBytes == 0) {
        return;
    }

    // --- Synchronisation primitives (match runCircuit() style) ---
    std::mutex memMtx;
    std::condition_variable memCv;
    int inFlightBlocks = 0;

    std::vector<bool> bufferBusy(maxBlocks, false);

    // Ensure buffers exist
    for (int i = 0; i < maxBlocks; ++i) {
        void* p = getAlignedBuffer(i);
        if (p == nullptr) {
            throw std::runtime_error("Aligned buffer missing");
        }
    }

    // Main loop: process each block
    for (int blockIdx = 0; blockIdx < numBlocksLocal; ++blockIdx) {
        int bufferIndex = blockIdx % maxBlocks;

        // Wait until buffer is free and capacity available
        {
            std::unique_lock<std::mutex> lock(memMtx);
            memCv.wait(lock, [&]() {
                return (inFlightBlocks < maxBlocks) && (!bufferBusy[bufferIndex]);
            });
            inFlightBlocks++;
            bufferBusy[bufferIndex] = true;
        }

        // Zero the buffer
        void* blockBuf = getAlignedBuffer(bufferIndex);
        std::memset(blockBuf, 0, blockBytes);

        // Create a Qureg that uses the block buffer (non-owning)
        Qureg tempQureg = createTempQureg(blockBuf, qubitsPerBlock);
        if (tempQureg.cpuAmps == nullptr) {
            {
                std::lock_guard<std::mutex> lk(memMtx);
                inFlightBlocks--;
                bufferBusy[bufferIndex] = false;
            }
            memCv.notify_all();
            throw std::runtime_error("createTempQureg produced null cpuAmps");
        }

        // Initialise zero or blank state depending on block index
        if (blockIdx == 0) {
            initZeroState(tempQureg);
        } else {
            initBlankState(tempQureg);
        }

        // Build chunkIndices vector
        std::vector<int> chunkIndices;
        chunkIndices.reserve(chunksPerBlock);
        for (int local = 0; local < chunksPerBlock; ++local) {
            chunkIndices.push_back(blockIdx * chunksPerBlock + local);
        }

        // Save block asynchronously
        saveBlock(blockIdx, chunkIndices,
                  [&, bufferIndex]() {
                      {
                          std::lock_guard<std::mutex> lock(memMtx);
                          inFlightBlocks--;
                          bufferBusy[bufferIndex] = false;
                      }
                      memCv.notify_all();
                  });
    }

    // Wait for all in-flight saves to finish
    {
        std::unique_lock<std::mutex> lock(memMtx);
        memCv.wait(lock, [&]() {
            return inFlightBlocks == 0;
        });
    }
}

void DiskBackedState::diskBacked_initPlusState() {
    // --- Preconditions & derived sizes (read-only) ---
    const int numBlocksLocal = getNumBlocks();
    const int qubitsPerBlock = getNumQubitsPerBlock();
    const int maxBlocks = getMaxBlocksInMemory();

    const size_t chunkBytes = ampsPerChunk * sizeof(qcomp);
    const size_t blockAmps = static_cast<size_t>(chunksPerBlock) * ampsPerChunk;
    const size_t blockBytes = blockAmps * sizeof(qcomp);

    const size_t numAmplitudes = static_cast<size_t>(numBlocksLocal) * blockAmps;
    const qcomp amp = 1.0 / std::sqrt(static_cast<double>(numAmplitudes));

    if (blockBytes == 0) {
        return;
    }

    // --- Synchronisation primitives ---
    std::mutex memMtx;
    std::condition_variable memCv;
    int inFlightBlocks = 0;

    std::vector<bool> bufferBusy(maxBlocks, false);

    // Ensure buffers exist
    for (int i = 0; i < maxBlocks; ++i) {
        void* p = getAlignedBuffer(i);
        if (p == nullptr) {
            throw std::runtime_error("Aligned buffer missing");
        }
    }

    // Main loop: process each block
    for (int blockIdx = 0; blockIdx < numBlocksLocal; ++blockIdx) {
        int bufferIndex = blockIdx % maxBlocks;

        // Wait until buffer is free and capacity available
        {
            std::unique_lock<std::mutex> lock(memMtx);
            memCv.wait(lock, [&]() {
                return (inFlightBlocks < maxBlocks) && (!bufferBusy[bufferIndex]);
            });
            inFlightBlocks++;
            bufferBusy[bufferIndex] = true;
        }

        // Zero the buffer
        void* blockBuf = getAlignedBuffer(bufferIndex);
        std::memset(blockBuf, 0, blockBytes);

        // Create a Qureg that uses the block buffer
        Qureg tempQureg = createTempQureg(blockBuf, qubitsPerBlock);
        if (tempQureg.cpuAmps == nullptr) {
            try { destroyQureg(tempQureg); } catch (...) {}
            {
                std::lock_guard<std::mutex> lk(memMtx);
                inFlightBlocks--;
                bufferBusy[bufferIndex] = false;
            }
            memCv.notify_all();
            throw std::runtime_error("createTempQureg produced null cpuAmps");
        }

        // Initialize plus state (uniform superposition)
        localiser_statevec_initUniformState(tempQureg, amp);

        // Build chunkIndices vector
        std::vector<int> chunkIndices;
        chunkIndices.reserve(chunksPerBlock);
        for (int local = 0; local < chunksPerBlock; ++local) {
            chunkIndices.push_back(blockIdx * chunksPerBlock + local);
        }

        // Save block asynchronously
        saveBlock(blockIdx, chunkIndices,
                  [&, bufferIndex, blockIdx]() {
                      {
                          std::lock_guard<std::mutex> lock(memMtx);
                          inFlightBlocks--;
                          bufferBusy[bufferIndex] = false;
                      }
                      memCv.notify_all();
                  });
    }

    // Wait for all in-flight saves to finish
    {
        std::unique_lock<std::mutex> lock(memMtx);
        memCv.wait(lock, [&]() {
            return inFlightBlocks == 0;
        });
    }
}

void DiskBackedState::diskBacked_initBlankState() {
    // --- Preconditions & derived sizes (read-only) ---
    const int numBlocksLocal = getNumBlocks();
    const int qubitsPerBlock = getNumQubitsPerBlock();
    const int maxBlocks = getMaxBlocksInMemory();

    const size_t blockAmps = static_cast<size_t>(chunksPerBlock) * ampsPerChunk;
    const size_t blockBytes = blockAmps * sizeof(qcomp);

    if (blockBytes == 0) {
        return;
    }

    // --- Synchronisation primitives ---
    std::mutex memMtx;
    std::condition_variable memCv;
    int inFlightBlocks = 0;

    std::vector<bool> bufferBusy(maxBlocks, false);

    // Ensure buffers exist
    for (int i = 0; i < maxBlocks; ++i) {
        void* p = getAlignedBuffer(i);
        if (p == nullptr) {
            throw std::runtime_error("Aligned buffer missing");
        }
    }

    // Main loop: process each block
    for (int blockIdx = 0; blockIdx < numBlocksLocal; ++blockIdx) {
        int bufferIndex = blockIdx % maxBlocks;

        // Wait until buffer is free and capacity available
        {
            std::unique_lock<std::mutex> lock(memMtx);
            memCv.wait(lock, [&]() {
                return (inFlightBlocks < maxBlocks) && (!bufferBusy[bufferIndex]);
            });
            inFlightBlocks++;
            bufferBusy[bufferIndex] = true;
        }

        // Zero the buffer
        void* blockBuf = getAlignedBuffer(bufferIndex);
        std::memset(blockBuf, 0, blockBytes);

        // Create a Qureg that uses the block buffer (non-owning)
        Qureg tempQureg = createTempQureg(blockBuf, qubitsPerBlock);
        if (tempQureg.cpuAmps == nullptr) {
            {
                std::lock_guard<std::mutex> lk(memMtx);
                inFlightBlocks--;
                bufferBusy[bufferIndex] = false;
            }
            memCv.notify_all();
            throw std::runtime_error("createTempQureg produced null cpuAmps");
        }

        // Initialise blank state
        initBlankState(tempQureg);

        // Build chunkIndices vector
        std::vector<int> chunkIndices;
        chunkIndices.reserve(chunksPerBlock);
        for (int local = 0; local < chunksPerBlock; ++local) {
            chunkIndices.push_back(blockIdx * chunksPerBlock + local);
        }

        // Save block asynchronously
        saveBlock(blockIdx, chunkIndices,
                  [&, bufferIndex]() {
                      {
                          std::lock_guard<std::mutex> lock(memMtx);
                          inFlightBlocks--;
                          bufferBusy[bufferIndex] = false;
                      }
                      memCv.notify_all();
                  });
    }

    // Wait for all in-flight saves to finish
    {
        std::unique_lock<std::mutex> lock(memMtx);
        memCv.wait(lock, [&]() {
            return inFlightBlocks == 0;
        });
    }
}

// TODO: optimise the probability function using initialisation format
qreal DiskBackedState::diskBacked_calcTotalProbability() const {
    void* alignedBuf = getChunkBuffer(0,0);
    qreal total = 0.0;
    for (size_t chunk = 0; chunk < numChunks; ++chunk) {
        loadChunk(chunk, alignedBuf);
        Qureg tempQureg = createTempQureg(alignedBuf, numQubitsPerChunk);
        total += calcTotalProb(tempQureg);
    }
    return total;
}

int DiskBackedState::diskBacked_applyQubitMeasurement(int qubit) {
    void* alignedBuf = getChunkBuffer(0,0);
    std::vector<qreal> probs(2, 0.0);

    if (qubit <= numQubitsPerChunk) {
        for (size_t chunk = 0; chunk < numChunks; ++chunk) {
            loadChunk(chunk, alignedBuf);
            Qureg tempQureg = createTempQureg(alignedBuf, numQubitsPerChunk);
            probs[0] += calcProbOfQubitOutcome(tempQureg, qubit, 0);
            probs[1] += calcProbOfQubitOutcome(tempQureg, qubit, 1);
        }
    } else {
        const auto& chunkMap = permTracker.getCurrentChunkMap();
        int qubitOffset = qubit - numQubitsPerChunk - 1;
        int regionSize = 1 << qubitOffset;

        for (size_t logicalChunk = 0; logicalChunk < numChunks; ++logicalChunk) {
            size_t physicalChunk = chunkMap[logicalChunk];
            loadChunk(physicalChunk, alignedBuf);
            Qureg tempQureg = createTempQureg(alignedBuf, numQubitsPerChunk);
            int region = logicalChunk / regionSize;
            int outcome = region % 2;
            probs[outcome] += calcTotalProb(tempQureg);
        }
    }

    int outcome = rand_getRandomSingleQubitOutcome(probs[0]);
    qreal correctProb = probs[outcome];
    qreal normalizationFactor = 1.0 / std::sqrt(correctProb);

    for (size_t logicalChunk = 0; logicalChunk < numChunks; ++logicalChunk) {
        size_t physicalChunk = permTracker.getCurrentChunkMap()[logicalChunk];
        loadChunk(physicalChunk, alignedBuf);

        qcomp* amps = static_cast<qcomp*>(alignedBuf);
        if (qubit <= numQubitsPerChunk) {
            for (size_t i = 0; i < ampsPerChunk; ++i) {
                bool isWrongOutcome = ((i >> qubit) & 1) != outcome;
                amps[i] = isWrongOutcome ? qcomp(0.0, 0.0) : amps[i] * normalizationFactor;
            }
        } else {
            int qubitOffset = qubit - numQubitsPerChunk - 1;
            int regionSize = 1 << qubitOffset;
            int region = logicalChunk / regionSize;
            int chunkOutcome = region % 2;
            for (size_t i = 0; i < ampsPerChunk; ++i) {
                if (chunkOutcome != outcome)
                    amps[i] = qcomp(0.0, 0.0);
                else
                    amps[i] *= normalizationFactor;
            }
        }
        saveChunk(physicalChunk, alignedBuf);
    }

    return outcome;
}

                                                                                       