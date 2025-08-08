#include "diskbackedstate.h"
#include "quest/include/qureg.h"
#include "quest/include/initialisations.h"
#include "quest/include/hook.h"
#include "quest/include/types.h"
#include <filesystem>
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <random>
#include <omp.h>
#include "quest/src/core/randomiser.hpp"

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
    numQubitsPerBlock = log2(chunksPerBlock *  ampsPerChunk);
    maxPermutableQubits = log2(chunksPerBlock);

    if (numAmplitudes % numChunks != 0) {
        throw std::runtime_error("Amplitudes must divide evenly across chunks");
    }

    generateChunkPaths();

    chunkFDs.resize(numChunks);
    for (size_t i = 0; i < numChunks; ++i) {
        int fd = open(chunkPaths[i].c_str(), O_RDWR | O_CREAT, 0644);
        if (fd < 0) {
                throw std::runtime_error(
        "Failed to open file: " + chunkPaths[i] +
        " errno=" + std::to_string(errno) + " (" + strerror(errno) + ")"
    );
        }
        chunkFDs[i] = fd;
    }

    const size_t totalBytes = ampsPerChunk * chunksPerBlock * sizeof(qcomp);
    for (int i = 0; i < maxBlocksInMemory; ++i) {
        void* buf = nullptr;
        if (posix_memalign(&buf, 4096, totalBytes) != 0)
            throw std::runtime_error("Buffer pool allocation failed");
        alignedBufferPool.push_back(buf);
    }
}

void* DiskBackedState::getAlignedBuffer(int idx) const {
    return alignedBufferPool.at(idx);
}

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
    if (!ioInitialised) {
        if (io_uring_queue_init((chunksPerBlock * maxBlocksInMemory * 2) + 1, &ring, 0) < 0) {
            throw std::runtime_error("Failed to initialize io_uring");
        }
        ioInitialised = true;

        // Start the IO completion thread
        ioRunning = true;
        ioCompletionThread = std::thread([this]() {
            this->ioCompletionLoop();
        });
    }

    if (!filesRegistered) {
        if (io_uring_register_files(&ring, chunkFDs.data(), chunkFDs.size()) < 0) {
            throw std::runtime_error("io_uring_register_files failed");
        }
        filesRegistered = true;
    }
}

void DiskBackedState::loadChunk(size_t chunkIndex, void* alignedBuf) const {
    if (chunkIndex >= numChunks)
        throw std::runtime_error("loadChunk: chunkIndex out of bounds");

    size_t chunkBytes = ampsPerChunk * sizeof(qcomp);
    int fd = chunkFDs[chunkIndex];

    ssize_t bytes = pread(fd, alignedBuf, chunkBytes, 0);
    if (bytes != (ssize_t)chunkBytes) {
        throw std::runtime_error("loadChunk: read failed, bytes=" + std::to_string(bytes));
    }
}

void DiskBackedState::saveChunk(size_t chunkIndex, void* alignedBuf) const {
    if (chunkIndex >= numChunks)
        throw std::runtime_error("saveChunk: chunkIndex out of bounds");

    size_t chunkBytes = ampsPerChunk * sizeof(qcomp);
    int fd = chunkFDs[chunkIndex];

    ssize_t bytes = pwrite(fd, alignedBuf, chunkBytes, 0);
    if (bytes != (ssize_t)chunkBytes) {
        throw std::runtime_error("saveChunk: write failed, bytes=" + std::to_string(bytes));
    }
}

struct IOContext {
    size_t chunkIndex;
    void* chunkBuf;
};

void printRingStatus(int submitted,  int ring_size) {
    printf("\r[IO Ring] In queue: %d / %d      ", submitted, ring_size);
    fflush(stdout);
}

void DiskBackedState::loadBlock(int blockIdx, const std::vector<int>& chunkIndices, void* alignedBuf, std::function<void()> onComplete) const {
    if (chunkIndices.size() != chunksPerBlock)
        throw std::runtime_error("loadBlock: chunkIndices size mismatch");

    ensureIoUringInitialised();
    const size_t chunkBytes = ampsPerChunk * sizeof(qcomp);

    auto* ctx = new BlockIOContext{
        blockIdx,
        chunkIndices,
        alignedBuf,
        0,
        std::move(onComplete)
    };

    int submittedChunks = 0;
    for (size_t i = 0; i < chunksPerBlock; ++i) {
        io_uring_sqe* sqe = io_uring_get_sqe(&ring);
        if (!sqe) {
            delete ctx;
            throw std::runtime_error("loadBlock: failed to get SQE");
        }

        void* chunkBuf = static_cast<char*>(alignedBuf) + i * chunkBytes;
        int fd = chunkFDs[chunkIndices[i]];

        io_uring_prep_read(sqe, fd, chunkBuf, chunkBytes, 0);
        io_uring_sqe_set_data(sqe, ctx);
        submittedChunks++;
    }

    if (submittedChunks == 0) {
        delete ctx;
        throw std::runtime_error("loadBlock: no SQEs submitted");
    }

    ctx->remainingChunks = submittedChunks;

    int submitted = io_uring_submit(&ring);
    if (submitted < 0) {
        delete ctx;
        throw std::runtime_error("loadBlock: submit failed");
    } else if (submitted < submittedChunks) {
        delete ctx;
        throw std::runtime_error("loadBlock: partial submit");
    }
    totalSubmitted += submitted;
}


void DiskBackedState::saveBlock(int blockIdx, const std::vector<int>& chunkIndices, void* alignedBuf, std::function<void()> onComplete) const {
    if (chunkIndices.size() != chunksPerBlock)
        throw std::runtime_error("saveBlock: chunkIndices size mismatch");

    ensureIoUringInitialised();
    const size_t chunkBytes = ampsPerChunk * sizeof(qcomp);

    auto* ctx = new BlockIOContext{
        blockIdx,
        chunkIndices,
        alignedBuf,
        0,
        std::move(onComplete)
    };

    int submittedChunks = 0;
    for (size_t i = 0; i < chunksPerBlock; ++i) {
        io_uring_sqe* sqe = io_uring_get_sqe(&ring);
        if (!sqe) {
            delete ctx;
            throw std::runtime_error("saveBlock: failed to get SQE");
        }

        void* chunkBuf = static_cast<char*>(alignedBuf) + i * chunkBytes;
        int fd = chunkFDs[chunkIndices[i]];

        io_uring_prep_write(sqe, fd, chunkBuf, chunkBytes, 0);
        io_uring_sqe_set_data(sqe, ctx);
        submittedChunks++;
    }

    if (submittedChunks == 0) {
        delete ctx;
        throw std::runtime_error("saveBlock: no SQEs submitted");
    }

    ctx->remainingChunks = submittedChunks;

    int submitted = io_uring_submit(&ring);
    if (submitted < 0) {
        delete ctx;
        throw std::runtime_error("saveBlock: submit failed");
    } else if (submitted < submittedChunks) {
        delete ctx;
        throw std::runtime_error("saveBlock: partial submit");
    }
    totalSubmitted += submitted;
}




void DiskBackedState::ioCompletionLoop() const {
    while (ioRunning) {
        io_uring_cqe* cqe;
        int ret = io_uring_wait_cqe(&ring, &cqe);
        if (ret < 0) continue;

        auto* ctx = static_cast<BlockIOContext*>(io_uring_cqe_get_data(cqe));
        io_uring_cqe_seen(&ring, cqe);

       --totalSubmitted;

        if (ctx) {
            if (cqe->res < 0) {
                fprintf(stderr, "\nIO error on block %d: %s\n", ctx->blockIdx, strerror(-cqe->res));
                std::terminate();
            }

            if (--ctx->remainingChunks == 0) {
                ctx->onComplete();
                delete ctx;
            }
        }

        printRingStatus(totalSubmitted, (chunksPerBlock * maxBlocksInMemory *4) + 1);
    }
}

DiskBackedState::~DiskBackedState() {
    ioRunning = false;

    if (ioInitialised) {
        // Wake up the completion loop with a no-op
        io_uring_sqe* sqe = io_uring_get_sqe(&ring);
        if (sqe) io_uring_prep_nop(sqe);
        io_uring_submit(&ring);

        if (ioCompletionThread.joinable())
            ioCompletionThread.join();

        io_uring_queue_exit(&ring);
    }

    for (int fd : chunkFDs) {
        if (fd >= 0) close(fd);
    }
    deleteAllChunkFiles();

    for (void* ptr : alignedBufferPool)
        free(ptr);
}


void DiskBackedState::deleteAllChunkFiles() {
    for (const auto& file : chunkPaths) {
        std::remove(file.c_str());
    }
}

// creates a Qureg for each chunk in the disk-backed state
void DiskBackedState::diskBacked_initRandomPureState() {
    for (size_t chunk = 0; chunk < numChunks; ++chunk) {
        Qureg chunkQureg = createForcedQureg(numQubitsPerChunk);
        initRandomPureState(chunkQureg);

        // initialising the pure states normalises to 1 for each chunk.
        // this then makes the total probability to be equal to the number
        // of chunks, so we need to renormalise each chunk to have a total
        // probability of 1 by dividing by the square root of the number of chunks.
        // TODO: is there a solution that preserves Haar randomness?
        double square_root = std::sqrt(numChunks);        
        for (qindex i = 0; i < ampsPerChunk; ++i) {
            chunkQureg.cpuAmps[i] /= square_root;
        }

        void* alignedBuf = getAlignedBuffer(chunk % maxBlocksInMemory);
        memcpy(alignedBuf, chunkQureg.cpuAmps, ampsPerChunk * sizeof(qcomp));
        saveChunk(chunk, alignedBuf);
        destroyQureg(chunkQureg); 
    }
}

void DiskBackedState::diskBacked_initZeroState() {
    for (int chunk = 0; chunk < numChunks; ++chunk) {
        void* alignedBuf = getAlignedBuffer(chunk % maxBlocksInMemory);
        qcomp* amps = static_cast<qcomp*>(alignedBuf);
        for (size_t i = 0; i < ampsPerChunk; ++i)
            amps[i] = 0.0;
        if (chunk == 0)
            amps[0] = 1.0;
        saveChunk(chunk, alignedBuf);
    }
}

void DiskBackedState::diskBacked_initPlusState() {
    qcomp amp = 1.0 / std::sqrt(numAmplitudes);
    for (int chunk = 0; chunk < numChunks; ++chunk) {
        void* alignedBuf = getAlignedBuffer(chunk % maxBlocksInMemory);
        qcomp* amps = static_cast<qcomp*>(alignedBuf);
        std::fill(amps, amps + ampsPerChunk, amp);
        saveChunk(chunk, alignedBuf);
    }
}

qreal DiskBackedState::diskBacked_calcTotalProbability() const {
    qreal total = 0.0;
    for (size_t i = 0; i < numChunks; ++i) {
        void* alignedBuf = getAlignedBuffer(i % maxBlocksInMemory);
        loadChunk(i, alignedBuf);
        Qureg tempQureg = createTempQureg(alignedBuf, numQubitsPerChunk);
        total += calcTotalProb(tempQureg);
    }
    return total;
}

int DiskBackedState::diskBacked_applyQubitMeasurement(int qubit) {
    std::vector<qreal> probs(2, 0.0);

    if (qubit <= numQubitsPerChunk) {
        for (size_t i = 0; i < numChunks; ++i) {
            void* alignedBuf = getAlignedBuffer(i % maxBlocksInMemory);
            loadChunk(i, alignedBuf);
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
            void* alignedBuf = getAlignedBuffer(logicalChunk % maxBlocksInMemory);
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
        void* alignedBuf = getAlignedBuffer(logicalChunk % maxBlocksInMemory);
        loadChunk(physicalChunk, alignedBuf);
        qcomp* amps = static_cast<qcomp*>(alignedBuf);

        if (qubit <= numQubitsPerChunk) {
            for (size_t i = 0; i < ampsPerChunk; ++i) {
                bool isWrongOutcome = ((i >> qubit) & 1) != outcome;
                if (isWrongOutcome) {
                    amps[i] = 0.0;
                } else {
                    amps[i] *= normalizationFactor;
                }
            }
        } else {
            int qubitOffset = qubit - numQubitsPerChunk - 1;
            int regionSize = 1 << qubitOffset;
            int region = logicalChunk / regionSize;
            int chunkOutcome = region % 2;

            if (chunkOutcome != outcome) {
                for (size_t i = 0; i < ampsPerChunk; ++i)
                    amps[i] = 0.0;
            } else {
                for (size_t i = 0; i < ampsPerChunk; ++i)
                    amps[i] *= normalizationFactor;
            }
        }

        saveChunk(physicalChunk, alignedBuf);
    }

    return outcome;
}

/*
void DiskBackedState::ensureIoUringInitialised() const {
    // Ensure ring exists and buffers are registered; tolerate races between threads.

    if (!ioInitialised.load()) {
        bool expected = false;
        if (ioInitialised.compare_exchange_strong(expected, true)) {
            // We won the race to initialise the ring and register buffers
            fprintf(stderr, "DEBUG: ensureIoUringInitialised: performing ring+buffer init\n");

            if (io_uring_queue_init((chunksPerBlock * maxBlocksInMemory * 2) + 1, &ring, 0) < 0) {
                ioInitialised.store(false);
                int err = errno;
                fprintf(stderr, "ERROR: io_uring_queue_init failed errno=%d (%s)\n", err, strerror(err));
                throw std::runtime_error("Failed to initialize io_uring");
            }

            // Register buffers on this ring (registeredIovecs should already be allocated by constructor)
            if (io_uring_register_buffers(&ring, registeredIovecs.data(),
                                          static_cast<unsigned int>(registeredIovecs.size())) < 0) {
                int err = errno;
                io_uring_queue_exit(&ring);
                ioInitialised.store(false);
                fprintf(stderr, "ERROR: io_uring_register_buffers failed errno=%d (%s)\n", err, strerror(err));
                throw std::runtime_error("io_uring_register_buffers failed");
            }

            fprintf(stderr, "DEBUG: ring and buffers registered successfully (ensure path)\n");
            ioInitialised.store(true);
        } else {
            // Another thread is initializing; wait until it's actually ready
            fprintf(stderr, "DEBUG: ensureIoUringInitialised: waiting for other init to finish\n");
            while (!ioInitialised.load()) {
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
            }
        }
    } else {
        // already initialised
        //fprintf(stderr, "DEBUG: ensureIoUringInitialised: ioInitialised already true\n");
    }

    // Ensure files are registered for fixed-file IO on this ring
    if (!filesRegistered.load()) {
        bool expectedFiles = false;
        if (filesRegistered.compare_exchange_strong(expectedFiles, true)) {
            fprintf(stderr, "DEBUG: ensureIoUringInitialised: registering files\n");
            if (io_uring_register_files(&ring, chunkFDs.data(), chunkFDs.size()) < 0) {
                int err = errno;
                filesRegistered.store(false);
                fprintf(stderr, "ERROR: io_uring_register_files failed errno=%d (%s)\n", err, strerror(err));
                throw std::runtime_error("io_uring_register_files failed");
            }
            fprintf(stderr, "DEBUG: files registered successfully\n");
        } else {
            fprintf(stderr, "DEBUG: ensureIoUringInitialised: waiting for files registration\n");
            while (!filesRegistered.load()) {
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
            }
        }
    } else {
        //fprintf(stderr, "DEBUG: ensureIoUringInitialised: filesRegistered already true\n");
    }

    // Ensure reaper thread is running (use atomic CAS)
    bool expectedRunning = false;
    if (ioRunning.compare_exchange_strong(expectedRunning, true)) {
        // we set ioRunning=true and start the thread
        fprintf(stderr, "DEBUG: ensureIoUringInitialised: starting ioCompletionThread\n");
        ioCompletionThread = std::thread([this]() {
            this->ioCompletionLoop();
        });
    } else {
        // someone else started it already
        //fprintf(stderr, "DEBUG: ensureIoUringInitialised: ioCompletionThread already running\n");
    }

// Ensure worker thread is running to process heavy callbacks
    if (!completionWorkerRunning.load()) {
        bool expectedWorker = false;
        if (completionWorkerRunning.compare_exchange_strong(expectedWorker, true)) {
            fprintf(stderr, "DEBUG: ensureIoUringInitialised: starting completionWorker\n");
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
                            fprintf(stderr, "Worker task threw: %s\n", e.what());
                        } catch (...) {
                            fprintf(stderr, "Worker task threw unknown exception\n");
                        }
                    }
                }
                fprintf(stderr, "DEBUG: completionWorker exiting\n");
            });
        } else {
            while (!completionWorkerRunning.load()) {
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
            }
        }
    } else {
        // already running
    }
} 

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
    io_uring_sqe_set_data(sqe, nullptr); // Explicitly mark as sync op

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

    io_uring_cqe_seen(&ring, cqe); // Consume our CQE
}

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
    io_uring_sqe_set_data(sqe, nullptr); // Explicitly mark as sync op

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

    io_uring_cqe_seen(&ring, cqe); // Consume our CQE
} 

void DiskBackedState::loadBlock(
    int blockIdx,
    const std::vector<int>& chunkIndices,
    std::function<void()> onComplete
) const {
    if (chunkIndices.size() != static_cast<size_t>(chunksPerBlock))
        throw std::runtime_error("loadBlock: chunkIndices size mismatch");

    ensureIoUringInitialised();
    const size_t chunkBytes = ampsPerChunk * sizeof(qcomp);

    // Allocate raw ctx and wrap in heap-allocated shared_ptr
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

    // Prepare SQEs
    for (int i = 0; i < chunksPerBlock; ++i) {
        struct io_uring_sqe* sqe = io_uring_get_sqe(&ring);
        if (!sqe) {
            delete sptrHeap;
            throw std::runtime_error("loadBlock: failed to get SQE");
        }

        int fileIndex = chunkIndices[i];                // index in registered files table
        int bufIndex  = memSlot * chunksPerBlock + i;   // index in registered buffers table

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

    // Set optimistic expected remaining chunks in raw ctx
    rawCtx->remainingChunks.store(submittedChunks);

    // *** IMPORTANT: increment outstanding count BEFORE calling submit to avoid races ***
    totalSubmitted.fetch_add(submittedChunks);

    // Submit
    int rc = io_uring_submit(&ring);
    if (rc < 0) {
        int err = -rc;
        // No SQEs accepted by kernel -> undo accounting and cleanup
        totalSubmitted.fetch_sub(submittedChunks);
        delete sptrHeap;
        fprintf(stderr, "ERROR: io_uring_submit returned %d (%s) in loadBlock for block %d\n", rc, strerror(err), blockIdx);
        throw std::runtime_error("loadBlock: submit failed");
    }

    if (rc == 0) {
        // Nothing accepted -> undo accounting and cleanup
        totalSubmitted.fetch_sub(submittedChunks);
        delete sptrHeap;
        fprintf(stderr, "ERROR: io_uring_submit accepted 0 SQEs in loadBlock for block %d\n", blockIdx);
        throw std::runtime_error("loadBlock: submit accepted 0 SQEs");
    }
 
if (rc < submittedChunks) {
        // Partial submit: kernel accepted fewer than we prepared.
        // Adjust accounting by removing the not-submitted count, and update ctx.
        int notSubmitted = submittedChunks - rc;
        totalSubmitted.fetch_sub(notSubmitted);
        rawCtx->remainingChunks.store(rc);
        fprintf(stderr, "WARNING: loadBlock partial submit: requested=%d submitted=%d for block %d. ctx retained until completions.\n",
                submittedChunks, rc, blockIdx);
    }

    // At this point our accounting totalSubmitted matches the number of SQEs the kernel accepted (rc).
    //fprintf(stderr, "DEBUG: loadBlock submitted %d SQEs for block %d (memSlot=%d). totalSubmitted=%d\n",
    //        rc, blockIdx, memSlot, totalSubmitted.load());
}

void DiskBackedState::saveBlock(
    int blockIdx,
    const std::vector<int>& chunkIndices,
    std::function<void()> onComplete
) const {
    if (chunkIndices.size() != static_cast<size_t>(chunksPerBlock))
        throw std::runtime_error("saveBlock: chunkIndices size mismatch");

    ensureIoUringInitialised();
    const size_t chunkBytes = ampsPerChunk * sizeof(qcomp);

    // Allocate raw ctx and wrap in heap-allocated shared_ptr
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

    // Prepare SQEs
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

    // Set optimistic expected remaining chunks in raw ctx
    rawCtx->remainingChunks.store(submittedChunks);

    // *** IMPORTANT: increment outstanding count BEFORE calling submit to avoid races ***
    totalSubmitted.fetch_add(submittedChunks);
 
// Submit
    int rc = io_uring_submit(&ring);
    if (rc < 0) {
        int err = -rc;
        totalSubmitted.fetch_sub(submittedChunks);
        delete sptrHeap;
        fprintf(stderr, "ERROR: io_uring_submit returned %d (%s) in saveBlock for block %d\n", rc, strerror(err), blockIdx);
        throw std::runtime_error("saveBlock: submit failed");
    }

    if (rc == 0) {
        totalSubmitted.fetch_sub(submittedChunks);
        delete sptrHeap;
        fprintf(stderr, "ERROR: io_uring_submit accepted 0 SQEs in saveBlock for block %d\n", blockIdx);
        throw std::runtime_error("saveBlock: submit accepted 0 SQEs");
    }

    if (rc < submittedChunks) {
        int notSubmitted = submittedChunks - rc;
        totalSubmitted.fetch_sub(notSubmitted);
        rawCtx->remainingChunks.store(rc);
       // fprintf(stderr, "WARNING: saveBlock partial submit: requested=%d submitted=%d for block %d. ctx retained until completions.\n",
               // submittedChunks, rc, blockIdx);
    }

    //fprintf(stderr, "DEBUG: saveBlock submitted %d SQEs for block %d (memSlot=%d). totalSubmitted=%d\n",
           // rc, blockIdx, memSlot, totalSubmitted.load());
}

void DiskBackedState::ioCompletionLoop() const {
    using namespace std::chrono;
    auto threadId = std::this_thread::get_id();
    //fprintf(stderr, "DEBUG: ioCompletionLoop started on thread %zu\n",
            std::hash<std::thread::id>{}(threadId));

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
                int err = -rc;
                fprintf(stderr, "DEBUG: io_uring_peek_cqe returned %d (%s). sleeping briefly\n", rc, strerror(err));
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
                // belongs to sync waiter — don't consume it here
                static int skipCount = 0;
                if ((++skipCount & 0x3FF) == 0) {
                    fprintf(stderr, "DEBUG: peeked CQE with NULL data (likely sync op). res=%d totalSubmitted=%d\n",
                            res, (int)totalSubmitted.load());
                }
                std::this_thread::sleep_for(noDataSleep);
                continue;
            }

            // This CQE corresponds to async BlockIO — consume it fully.
           // fprintf(stderr, "DEBUG: CQE (async) arrived: res=%d totalSubmitted=%d data=%p\n",
                   // res, (int)totalSubmitted.load(), data);

            // Remove the CQE from the CQ (we own it)
            io_uring_cqe_seen(&ring, cqe);

            // Update global submitted counter
            int prevTotal = totalSubmitted.fetch_sub(1);
            int newTotal = prevTotal - 1;
            if (newTotal < 0) {
                fprintf(stderr, "WARNING: totalSubmitted went negative: %d (prev=%d)\n", newTotal, prevTotal);
            }

            // Interpret userdata as pointer-to-shared_ptr<BlockIOContext>
            auto* sptrHeap = static_cast<std::shared_ptr<BlockIOContext>*>(data);
            if (!sptrHeap) {
                fprintf(stderr, "DEBUG: sptrHeap == nullptr (unexpected)\n");
                continue;
            }

            // Make a local owning copy of the shared_ptr so ctx stays alive while we process this CQE.
            std::shared_ptr<BlockIOContext> ctx_sp = *sptrHeap;
            BlockIOContext* ctx = ctx_sp.get();
            if (!ctx) {
                fprintf(stderr, "DEBUG: ctx_sp.get() returned nullptr (unexpected)\n");
                continue;
            }
 // Sanity-checks to detect corruption early
            int ctxBlock = ctx->blockIdx;
            size_t ctxChunkCount = ctx->chunkIndices.size();
            if (ctxBlock < 0 || ctxBlock >= numBlocks) {
                fprintf(stderr, "CORRUPTION DETECTED: ctx->blockIdx out of range: %d (numBlocks=%d). ctx=%p\n",
                        ctxBlock, numBlocks, ctx);
                // Don't delete sptrHeap here; continue to avoid double-free.
                continue;
            }
            if (ctxChunkCount == 0 || ctxChunkCount > static_cast<size_t>(chunksPerBlock) * 1000) {
                fprintf(stderr, "CORRUPTION DETECTED: ctx->chunkIndices.size() suspicious: %zu (chunksPerBlock=%d). ctx=%p\n",
                        ctxChunkCount, chunksPerBlock, ctx);
                continue;
            }

           // fprintf(stderr, "DEBUG: CQE for block %d (remainingChunks before update = %d). chunkCount=%zu\n",
             //       ctxBlock, ctx->remainingChunks.load(), ctxChunkCount);

            if (res < 0) {
                int err = -res;
               // fprintf(stderr, "ERROR: IO error on block %d: res=%d (%s)\n", ctxBlock, res, strerror(err));
                std::terminate();
            }

            // Atomically decrement remainingChunks and check if block-level IO is done.
            int prev = ctx->remainingChunks.fetch_sub(1);
            int rem = prev - 1;
            //fprintf(stderr, "DEBUG: After decrement remainingChunks=%d for block %d\n", rem, ctxBlock);

            if (rem < 0) {
                fprintf(stderr, "CORRUPTION WARNING: remainingChunks negative for block %d (prev=%d). ctx=%p\n",
                        ctxBlock, prev, ctx);
                continue;
            }

            if (rem == 0) {
                // Final completion for this ctx.
                // Extract callback (still via ctx_sp), schedule worker task, then free the heap-allocated shared_ptr.
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
                   // fprintf(stderr, "DEBUG: Enqueued onComplete for block %d to worker queue\n", ctxBlock);
                } else {
                    //fprintf(stderr, "DEBUG: No onComplete callback to enqueue for block %d\n", ctxBlock);
                }

                // Now safe to delete the heap-allocated shared_ptr: this removes the last owner we created at submit time.
                delete sptrHeap; // this does NOT delete ctx immediately if there are other shared_ptr copies; otherwise it will free ctx.
            }

            // Optionally print ring status
            printRingStatus((int)totalSubmitted.load(), (chunksPerBlock * maxBlocksInMemory * 2) + 1);

        } catch (const std::exception &e) {
            fprintf(stderr, "ERROR: exception in ioCompletionLoop: %s\n", e.what());
            std::this_thread::sleep_for(milliseconds(10));
        } catch (...) {
            fprintf(stderr, "ERROR: unknown exception in ioCompletionLoop\n");
            std::this_thread::sleep_for(milliseconds(10));
        }
    }

    //fprintf(stderr, "DEBUG: ioCompletionLoop exiting\n");
}*/