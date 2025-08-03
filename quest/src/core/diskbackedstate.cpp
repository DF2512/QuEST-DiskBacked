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
}

void DiskBackedState::loadChunk(size_t chunkIndex, void* alignedBuf, std::vector<qcomp>& buffer) const {
    if (chunkIndex >= numChunks)
        throw std::runtime_error("loadChunk: chunkIndex out of bounds");

    ensureIoUringInitialised();

    const std::string& path = chunkPaths[chunkIndex];
    size_t chunkBytes = ampsPerChunk * sizeof(qcomp);

    int fd = open(path.c_str(), O_RDONLY);
    if (fd < 0) {
        throw std::runtime_error("loadChunk: open failed " + path);
    }

    io_uring_sqe* sqe = io_uring_get_sqe(&ring);
    if (!sqe) {
        close(fd);
        throw std::runtime_error("loadChunk: failed to get SQE");
    }

    io_uring_prep_read(sqe, fd, alignedBuf, chunkBytes, 0);
    io_uring_sqe_set_data(sqe, nullptr);

    if (io_uring_submit_and_wait(&ring, 1) < 0) {
        close(fd);
        throw std::runtime_error("loadChunk: failed to submit_and_wait");
    }

    io_uring_cqe* cqe;
    if (io_uring_wait_cqe(&ring, &cqe) < 0) {
        close(fd);
        throw std::runtime_error("loadChunk: failed to wait for CQE");
    }

    if (cqe->res != (int)chunkBytes) {
        std::string err = "loadChunk: I/O error: expected " + std::to_string(chunkBytes) + ", got " + std::to_string(cqe->res);
        io_uring_cqe_seen(&ring, cqe);
        close(fd);
        throw std::runtime_error(err);
    }

    io_uring_cqe_seen(&ring, cqe);
    close(fd);

    buffer.resize(ampsPerChunk);
    std::memcpy(buffer.data(), alignedBuf, chunkBytes);
}


void DiskBackedState::saveChunk(size_t chunkIndex, void* alignedBuf, const std::vector<qcomp>& buffer) const {
    if (chunkIndex >= numChunks)
        throw std::runtime_error("saveChunk: chunkIndex out of bounds");

    if (buffer.size() != ampsPerChunk)
        throw std::runtime_error("saveChunk: buffer size mismatch");

    ensureIoUringInitialised();

    const std::string& path = chunkPaths[chunkIndex];
    size_t chunkBytes = ampsPerChunk * sizeof(qcomp);

    std::memcpy(alignedBuf, buffer.data(), chunkBytes);

    int fd = open(path.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
    if (fd < 0) {
        throw std::runtime_error("saveChunk: open failed " + path);
    }

    io_uring_sqe* sqe = io_uring_get_sqe(&ring);
    if (!sqe) {
        close(fd);
        throw std::runtime_error("saveChunk: failed to get SQE");
    }

    io_uring_prep_write(sqe, fd, alignedBuf, chunkBytes, 0);
    io_uring_sqe_set_data(sqe, nullptr);

    if (io_uring_submit_and_wait(&ring, 1) < 0) {
        close(fd);
        throw std::runtime_error("saveChunk: failed to submit_and_wait");
    }

    io_uring_cqe* cqe;
    if (io_uring_wait_cqe(&ring, &cqe) < 0) {
        close(fd);
        throw std::runtime_error("saveChunk: failed to wait for CQE");
    }

    if (cqe->res != (int)chunkBytes) {
        std::string err = "saveChunk: I/O error: expected " + std::to_string(chunkBytes) + ", got " + std::to_string(cqe->res);
        io_uring_cqe_seen(&ring, cqe);
        close(fd);
        throw std::runtime_error(err);
    }

    io_uring_cqe_seen(&ring, cqe);
    close(fd);
}

/*
void DiskBackedState::loadBlock(int blockIdx, const std::vector<int>& chunkIndices, std::vector<qcomp>& buffer) const {
    if (chunkIndices.size() != chunksPerBlock) {
        throw std::runtime_error("loadBlock: chunkIndices size mismatch");
    }

    buffer.resize(ampsPerChunk * chunksPerBlock);

    for (size_t i = 0; i < chunksPerBlock; ++i) {
        int chunkIdx = chunkIndices[i];
        if (chunkIdx < 0 || chunkIdx >= (int)numChunks) {
            throw std::runtime_error("loadBlock: invalid chunk index");
        }

        std::vector<qcomp> chunkBuf;
        loadChunk(chunkIdx, chunkBuf);

        std::copy(chunkBuf.begin(), chunkBuf.end(),
                  buffer.begin() + i * ampsPerChunk);
    }

    //std::cout << "[Disk] Loaded block " << blockIdx << " (chunks: ";
    //for (int c : chunkIndices) std::cout << c << " ";
    //std::cout << ")\n";
}
*/
void DiskBackedState::loadBlock(int blockIdx, const std::vector<int>& chunkIndices, void* alignedBuf, std::vector<qcomp>& buffer) const {
    if (chunkIndices.size() != chunksPerBlock) {
        throw std::runtime_error("loadBlock: chunkIndices size mismatch");
    }

    ensureIoUringInitialised();

    const size_t chunkBytes = ampsPerChunk * sizeof(qcomp);
    const size_t totalBytes = chunkBytes * chunksPerBlock;

    std::vector<int> fds(chunksPerBlock);
    for (size_t i = 0; i < chunksPerBlock; ++i) {
        const std::string& path = chunkPaths[chunkIndices[i]];
        int fd = open(path.c_str(), O_RDONLY);
        if (fd < 0) {
            throw std::runtime_error("loadBlock: failed to open file: " + path);
        }
        fds[i] = fd;
    }

    for (size_t i = 0; i < chunksPerBlock; ++i) {
        io_uring_sqe* sqe = io_uring_get_sqe(&ring);
        if (!sqe) {
            for (int fd : fds) close(fd);
            throw std::runtime_error("loadBlock: failed to get SQE");
        }

        void* chunkBuf = static_cast<char*>(alignedBuf) + i * chunkBytes;
        io_uring_prep_read(sqe, fds[i], chunkBuf, chunkBytes, 0);
        io_uring_sqe_set_data(sqe, reinterpret_cast<void*>(i));
    }

    if (io_uring_submit_and_wait(&ring, chunksPerBlock) < 0) {
        for (int fd : fds) close(fd);
        throw std::runtime_error("loadBlock: failed to submit_and_wait");
    }

    for (size_t i = 0; i < chunksPerBlock; ++i) {
        io_uring_cqe* cqe;
        if (io_uring_wait_cqe(&ring, &cqe) < 0) {
            for (int fd : fds) close(fd);
            throw std::runtime_error("loadBlock: failed to wait for CQE");
        }
        if (cqe->res != (int)chunkBytes) {
            std::string err = "loadBlock: I/O error: expected " + std::to_string(chunkBytes) + ", got " + std::to_string(cqe->res);
            io_uring_cqe_seen(&ring, cqe);
            for (int fd : fds) close(fd);
            throw std::runtime_error(err);
        }
        io_uring_cqe_seen(&ring, cqe);
    }

    buffer.resize(ampsPerChunk * chunksPerBlock);
    std::memcpy(buffer.data(), alignedBuf, totalBytes);

    for (int fd : fds) close(fd);
}

/*
// saves a block of chunks to disk files
void DiskBackedState::saveBlock(int blockIdx, const std::vector<int>& chunkIndices, const std::vector<qcomp>& buffer) const {
    if (chunkIndices.size() != chunksPerBlock) {
        throw std::runtime_error("saveBlock: chunkIndices size mismatch");
    }

    if (buffer.size() != ampsPerChunk * chunksPerBlock) {
        throw std::runtime_error("saveBlock: buffer size mismatch");
    }

    for (size_t i = 0; i < chunksPerBlock; ++i) {
        int chunkIdx = chunkIndices[i];
        if (chunkIdx < 0 || chunkIdx >= (int)numChunks) {
            throw std::runtime_error("saveBlock: invalid chunk index");
        }

        std::vector<qcomp> chunkBuf(ampsPerChunk);
        std::copy(buffer.begin() + i * ampsPerChunk,
                  buffer.begin() + (i + 1) * ampsPerChunk,
                  chunkBuf.begin());

        saveChunk(chunkIdx, chunkBuf);
    }

    //std::cout << "[Disk] Saved block " << blockIdx << " (chunks: ";
    //for (int c : chunkIndices) std::cout << c << " ";
    //std::cout << ")\n";
}
*/
void DiskBackedState::saveBlock(int blockIdx, const std::vector<int>& chunkIndices, void* alignedBuf, const std::vector<qcomp>& buffer) const {
    if (chunkIndices.size() != chunksPerBlock) {
        throw std::runtime_error("saveBlock: chunkIndices size mismatch");
    }

    if (buffer.size() != ampsPerChunk * chunksPerBlock) {
        throw std::runtime_error("saveBlock: buffer size mismatch");
    }

    ensureIoUringInitialised();

    const size_t chunkBytes = ampsPerChunk * sizeof(qcomp);
    const size_t totalBytes = chunkBytes * chunksPerBlock;

    std::memcpy(alignedBuf, buffer.data(), totalBytes);

    std::vector<int> fds(chunksPerBlock);
    for (size_t i = 0; i < chunksPerBlock; ++i) {
        const std::string& path = chunkPaths[chunkIndices[i]];
        int fd = open(path.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
        if (fd < 0) {
            throw std::runtime_error("saveBlock: failed to open file: " + path);
        }
        fds[i] = fd;
    }

    for (size_t i = 0; i < chunksPerBlock; ++i) {
        io_uring_sqe* sqe = io_uring_get_sqe(&ring);
        if (!sqe) {
            for (int fd : fds) close(fd);
            throw std::runtime_error("saveBlock: failed to get SQE");
        }

        void* chunkBuf = static_cast<char*>(alignedBuf) + i * chunkBytes;
        io_uring_prep_write(sqe, fds[i], chunkBuf, chunkBytes, 0);
        io_uring_sqe_set_data(sqe, reinterpret_cast<void*>(i));
    }

    if (io_uring_submit_and_wait(&ring, chunksPerBlock) < 0) {
        for (int fd : fds) close(fd);
        throw std::runtime_error("saveBlock: failed to submit_and_wait");
    }

    for (size_t i = 0; i < chunksPerBlock; ++i) {
        io_uring_cqe* cqe;
        if (io_uring_wait_cqe(&ring, &cqe) < 0) {
            for (int fd : fds) close(fd);
            throw std::runtime_error("saveBlock: failed to wait for CQE");
        }
        if (cqe->res != (int)chunkBytes) {
            std::string err = "saveBlock: I/O error: expected " + std::to_string(chunkBytes) + ", got " + std::to_string(cqe->res);
            io_uring_cqe_seen(&ring, cqe);
            for (int fd : fds) close(fd);
            throw std::runtime_error(err);
        }
        io_uring_cqe_seen(&ring, cqe);
    }

    for (int fd : fds) close(fd);
}


DiskBackedState::~DiskBackedState() {
    deleteAllChunkFiles();
    if (ioInitialised) {
        io_uring_queue_exit(&ring);
    }
    for (void* ptr : ioAlignedBuffers)
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

        std::vector<qcomp> buffer(chunkQureg.cpuAmps, chunkQureg.cpuAmps + ampsPerChunk);
        void* alignedBuf = getAlignedBuffer(chunk % maxBlocksInMemory);
        saveChunk(chunk, alignedBuf, buffer);
        destroyQureg(chunkQureg); 
    }
}

void DiskBackedState::diskBacked_initZeroState() {
    for (int chunk = 0; chunk < numChunks; ++chunk) {
        std::vector<qcomp> buffer(ampsPerChunk, 0.0);
        if (chunk == 0) {
            buffer[0] = 1.0;
        }
        void* alignedBuf = getAlignedBuffer(chunk % maxBlocksInMemory);
        saveChunk(chunk, alignedBuf, buffer);
    }
}

qreal DiskBackedState::diskBacked_calcTotalProbability() const {
    qreal total = 0.0;
    for (size_t i = 0; i < numChunks; ++i) {
        void* alignedBuf = getAlignedBuffer(i % maxBlocksInMemory);
        std::vector<qcomp> buffer;
        loadChunk(i, alignedBuf, buffer);
        Qureg tempQureg = createTempQureg(buffer, numQubitsPerChunk);
        total += calcTotalProb(tempQureg);
    }
    return total;
}

int DiskBackedState::diskBacked_applyQubitMeasurement(int qubit) {
    std::vector<qreal> probs(2, 0.0);

    if (qubit <= numQubitsPerChunk) {
        for (size_t i = 0; i < numChunks; ++i) {
            void* alignedBuf = getAlignedBuffer(i % maxBlocksInMemory);
            std::vector<qcomp> buffer;
            loadChunk(i, alignedBuf, buffer);
            Qureg tempQureg = createTempQureg(buffer, numQubitsPerChunk);
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
            std::vector<qcomp> buffer;
            loadChunk(physicalChunk, alignedBuf, buffer);
            Qureg tempQureg = createTempQureg(buffer, numQubitsPerChunk);

            int region = logicalChunk / regionSize;
            int outcome = region % 2;
            qreal chunkProb = calcTotalProb(tempQureg);
            probs[outcome] += chunkProb;
        }
    }

    int outcome = rand_getRandomSingleQubitOutcome(probs[0]);
    qreal correctProb = probs[outcome];
    qreal normalizationFactor = 1.0 / std::sqrt(correctProb);

    for (size_t logicalChunk = 0; logicalChunk < numChunks; ++logicalChunk) {
        size_t physicalChunk = permTracker.getCurrentChunkMap()[logicalChunk];
        void* alignedBuf = getAlignedBuffer(logicalChunk % maxBlocksInMemory);
        std::vector<qcomp> buffer;
        loadChunk(physicalChunk, alignedBuf, buffer);

        if (qubit <= numQubitsPerChunk) {
            for (size_t i = 0; i < buffer.size(); ++i) {
                bool isWrongOutcome = ((i >> qubit) & 1) != outcome;
                if (isWrongOutcome) {
                    buffer[i] = 0.0;
                } else {
                    buffer[i] *= normalizationFactor;
                }
            }
        } else {
            int qubitOffset = qubit - numQubitsPerChunk - 1;
            int regionSize = 1 << qubitOffset;
            int region = logicalChunk / regionSize;
            int chunkOutcome = region % 2;

            if (chunkOutcome != outcome) {
                std::fill(buffer.begin(), buffer.end(), 0.0);
            } else {
                for (auto& amp : buffer) {
                    amp *= normalizationFactor;
                }
            }
        }

        saveChunk(physicalChunk, alignedBuf, buffer);
    }

    return outcome;
}
