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
        // Increase ring size significantly to handle more concurrent operations
        if (io_uring_queue_init(512, &ring, 0) < 0) {
            throw std::runtime_error("Failed to initialize io_uring");
        }
        ioInitialised = true;
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

    ensureIoUringInitialised();
    const size_t chunkBytes = ampsPerChunk * sizeof(qcomp);

    io_uring_sqe* sqe = io_uring_get_sqe(&ring);
    if (!sqe) throw std::runtime_error("loadChunk: failed to get SQE");

    int fd = chunkFDs[chunkIndex];
    io_uring_prep_read(sqe, fd, alignedBuf, chunkBytes, 0);
    io_uring_sqe_set_data(sqe, nullptr);

    if (io_uring_submit_and_wait(&ring, 1) < 0)
        throw std::runtime_error("loadChunk: submit_and_wait failed");

    io_uring_cqe* cqe;
    if (io_uring_wait_cqe(&ring, &cqe) < 0)
        throw std::runtime_error("loadChunk: wait_cqe failed");

    if (cqe->res != (int)chunkBytes) {
        std::string err = "loadChunk: I/O error: expected " +
                          std::to_string(chunkBytes) + ", got " +
                          std::to_string(cqe->res);
        io_uring_cqe_seen(&ring, cqe);
        throw std::runtime_error(err);
    }

    io_uring_cqe_seen(&ring, cqe);
}

void DiskBackedState::saveChunk(size_t chunkIndex, const void* alignedBuf) const {
    if (chunkIndex >= numChunks)
        throw std::runtime_error("saveChunk: chunkIndex out of bounds");

    ensureIoUringInitialised();
    const size_t chunkBytes = ampsPerChunk * sizeof(qcomp);

    io_uring_sqe* sqe = io_uring_get_sqe(&ring);
    if (!sqe) throw std::runtime_error("saveChunk: failed to get SQE");

    int fd = chunkFDs[chunkIndex];
    io_uring_prep_write(sqe, fd, alignedBuf, chunkBytes, 0);
    io_uring_sqe_set_data(sqe, nullptr);

    if (io_uring_submit_and_wait(&ring, 1) < 0)
        throw std::runtime_error("saveChunk: submit_and_wait failed");

    io_uring_cqe* cqe;
    if (io_uring_wait_cqe(&ring, &cqe) < 0)
        throw std::runtime_error("saveChunk: wait_cqe failed");

    if (cqe->res != (int)chunkBytes) {
        std::string err = "saveChunk: I/O error: expected " +
                          std::to_string(chunkBytes) + ", got " +
                          std::to_string(cqe->res);
        io_uring_cqe_seen(&ring, cqe);
        throw std::runtime_error(err);
    }

    io_uring_cqe_seen(&ring, cqe);
}

struct IOContext {
    size_t chunkIndex;
    void* chunkBuf;
};

void DiskBackedState::loadBlock(int blockIdx, const std::vector<int>& chunkIndices, void* alignedBuf) const {
    if (chunkIndices.size() != chunksPerBlock)
        throw std::runtime_error("loadBlock: chunkIndices size mismatch");

    ensureIoUringInitialised();
    const size_t chunkBytes = ampsPerChunk * sizeof(qcomp);

    std::vector<IOContext> ioContexts(chunksPerBlock);  // stack-safe
    for (size_t i = 0; i < chunksPerBlock; ++i) {
        io_uring_sqe* sqe = io_uring_get_sqe(&ring);
        if (!sqe) throw std::runtime_error("loadBlock: failed to get SQE");

        void* chunkBuf = static_cast<char*>(alignedBuf) + i * chunkBytes;
        int fd = chunkFDs[chunkIndices[i]];

        ioContexts[i] = IOContext{chunkIndices[i], chunkBuf};
        io_uring_prep_read(sqe, fd, chunkBuf, chunkBytes, 0);
        io_uring_sqe_set_data(sqe, &ioContexts[i]);
    }

    if (io_uring_submit_and_wait(&ring, chunksPerBlock) < 0)
        throw std::runtime_error("loadBlock: submit_and_wait failed");

    for (size_t i = 0; i < chunksPerBlock; ++i) {
        io_uring_cqe* cqe;
        if (io_uring_wait_cqe(&ring, &cqe) < 0)
            throw std::runtime_error("loadBlock: wait_cqe failed");

        IOContext* ctx = static_cast<IOContext*>(io_uring_cqe_get_data(cqe));

        if (cqe->res != (int)chunkBytes) {
            std::string err = "loadBlock: I/O error on chunkIndex " +
                              std::to_string(ctx->chunkIndex) + ": expected " +
                              std::to_string(chunkBytes) + ", got " +
                              std::to_string(cqe->res);
            io_uring_cqe_seen(&ring, cqe);
            throw std::runtime_error(err);
        }

        io_uring_cqe_seen(&ring, cqe);
    }
}


void DiskBackedState::saveBlock(int blockIdx, const std::vector<int>& chunkIndices, const void* alignedBuf) const {
    if (chunkIndices.size() != chunksPerBlock)
        throw std::runtime_error("saveBlock: chunkIndices size mismatch");

    ensureIoUringInitialised();
    const size_t chunkBytes = ampsPerChunk * sizeof(qcomp);

    std::vector<IOContext> ioContexts(chunksPerBlock);  // on stack
    for (size_t i = 0; i < chunksPerBlock; ++i) {
        io_uring_sqe* sqe = io_uring_get_sqe(&ring);
        if (!sqe) throw std::runtime_error("saveBlock: failed to get SQE");

        const void* chunkBuf = static_cast<const char*>(alignedBuf) + i * chunkBytes;
        int fd = chunkFDs[chunkIndices[i]];

        ioContexts[i] = IOContext{chunkIndices[i], const_cast<void*>(chunkBuf)};
        io_uring_prep_write(sqe, fd, chunkBuf, chunkBytes, 0);
        io_uring_sqe_set_data(sqe, &ioContexts[i]);
    }

    if (io_uring_submit_and_wait(&ring, chunksPerBlock) < 0)
        throw std::runtime_error("saveBlock: submit_and_wait failed");

    for (size_t i = 0; i < chunksPerBlock; ++i) {
        io_uring_cqe* cqe;
        if (io_uring_wait_cqe(&ring, &cqe) < 0)
            throw std::runtime_error("saveBlock: wait_cqe failed");

        IOContext* ctx = static_cast<IOContext*>(io_uring_cqe_get_data(cqe));

        if (cqe->res != (int)chunkBytes) {
            std::string err = "saveBlock: I/O error on chunkIndex " +
                              std::to_string(ctx->chunkIndex) + ": expected " +
                              std::to_string(chunkBytes) + ", got " +
                              std::to_string(cqe->res);
            io_uring_cqe_seen(&ring, cqe);
            throw std::runtime_error(err);
        }

        io_uring_cqe_seen(&ring, cqe);
    }
}


DiskBackedState::~DiskBackedState() {
    for (int fd : chunkFDs) {
        if (fd >= 0) close(fd);
    }
    deleteAllChunkFiles();
    if (ioInitialised) {
        io_uring_queue_exit(&ring);
    }
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
