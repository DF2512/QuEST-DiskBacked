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
                                 const std::vector<std::string>& diskRoots_)
    : numQubits(numQubits_),
      numBlocks(numBlocks_),
      chunksPerBlock(chunksPerBlock_),
      diskRoots(diskRoots_),
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
size_t DiskBackedState::getNumAmplitudes() const { return numAmplitudes; }
size_t DiskBackedState::getNumChunks() const { return numChunks; }
size_t DiskBackedState::getAmpsPerChunk() const { return ampsPerChunk; }

void DiskBackedState::ensureIoUringInitialised() const {
    if (!ioInitialised) {
        if (io_uring_queue_init(32, &ring, 0) < 0) {
            throw std::runtime_error("Failed to initialize io_uring");
        }
        ioInitialised = true;
    }
}

void DiskBackedState::loadChunk(size_t chunkIndex, std::vector<qcomp>& buffer) const {
    if (chunkIndex >= numChunks)
        throw std::runtime_error("loadChunk: chunkIndex out of bounds");

    ensureIoUringInitialised(); 

    const std::string& path = chunkPaths[chunkIndex];
    size_t chunkBytes = ampsPerChunk * sizeof(qcomp);

    void* rawBuf;
    if (posix_memalign(&rawBuf, 4096, chunkBytes) != 0)
        throw std::runtime_error("loadChunk: posix_memalign failed");

    int fd = open(path.c_str(), O_RDONLY | O_DIRECT);
    if (fd < 0) {
        free(rawBuf);
        throw std::runtime_error("loadChunk: open failed " + path);
    }

    io_uring_sqe* sqe = io_uring_get_sqe(&ring);
    io_uring_prep_read(sqe, fd, rawBuf, chunkBytes, 0);
    io_uring_submit(&ring);

    io_uring_cqe* cqe;
    io_uring_wait_cqe(&ring, &cqe);
    if (cqe->res < 0) {
        std::string err = "loadChunk: I/O error: ";
        err += strerror(-cqe->res);
        io_uring_cqe_seen(&ring, cqe);
        close(fd); free(rawBuf);
        throw std::runtime_error(err);
    }
    io_uring_cqe_seen(&ring, cqe);

    buffer.resize(ampsPerChunk);
    std::memcpy(buffer.data(), rawBuf, chunkBytes);

    close(fd);
    free(rawBuf);
}


void DiskBackedState::saveChunk(size_t chunkIndex, const std::vector<qcomp>& buffer) const {
    if (chunkIndex >= numChunks)
        throw std::runtime_error("saveChunk: chunkIndex out of bounds");

    if (buffer.size() != ampsPerChunk)
        throw std::runtime_error("saveChunk: buffer size mismatch");

    ensureIoUringInitialised();

    const std::string& path = chunkPaths[chunkIndex];
    size_t chunkBytes = ampsPerChunk * sizeof(qcomp);

    void* rawBuf;
    if (posix_memalign(&rawBuf, 4096, chunkBytes) != 0)
        throw std::runtime_error("saveChunk: posix_memalign failed");

    std::memcpy(rawBuf, buffer.data(), chunkBytes);

    int fd = open(path.c_str(), O_WRONLY | O_CREAT | O_TRUNC | O_DIRECT, 0644);
    if (fd < 0) {
        free(rawBuf);
        throw std::runtime_error("saveChunk: open failed " + path);
    }

    io_uring_sqe* sqe = io_uring_get_sqe(&ring);
    io_uring_prep_write(sqe, fd, rawBuf, chunkBytes, 0);
    io_uring_submit(&ring);

    io_uring_cqe* cqe;
    io_uring_wait_cqe(&ring, &cqe);
    if (cqe->res < 0) {
        std::string err = "saveChunk: I/O error: ";
        err += strerror(-cqe->res);
        io_uring_cqe_seen(&ring, cqe);
        close(fd); free(rawBuf);
        throw std::runtime_error(err);
    }
    io_uring_cqe_seen(&ring, cqe);

    close(fd);
    free(rawBuf);
}


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

DiskBackedState::~DiskBackedState() {
    deleteAllChunkFiles();
    if (ioInitialized) {
        io_uring_queue_exit(&ring);
    }
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
        saveChunk(chunk, buffer);
        destroyQureg(chunkQureg); 
    }
}

void DiskBackedState::diskBacked_initZeroState() {
    for (int chunk = 0; chunk < numChunks; ++chunk) {
        std::vector<qcomp> buffer(ampsPerChunk, 0.0);
        if (chunk == 0) {
            buffer[0] = 1.0;
        }
        saveChunk(chunk, buffer);
    }
}

qreal DiskBackedState::diskBacked_calcTotalProbability() const {
    qreal total = 0.0;
    
    for (size_t i = 0; i < numChunks; ++i) {
        std::vector<qcomp> buffer;
        loadChunk(i, buffer);
        
        Qureg tempQureg = createTempQureg(buffer, numQubitsPerChunk);

        
        total += calcTotalProb(tempQureg);
        
    }

    return total;
}

int DiskBackedState::diskBacked_applyQubitMeasurement(int qubit) {
    std::vector<qreal> probs(2, 0.0);
    
    if (qubit <= numQubitsPerChunk) {
        // Qubit is within chunk size - simple case
        for (size_t i = 0; i < numChunks; ++i) {
            std::vector<qcomp> buffer;
            loadChunk(i, buffer);
            Qureg tempQureg = createTempQureg(buffer, numQubitsPerChunk);
            probs[0] += calcProbOfQubitOutcome(tempQureg, qubit, 0);
            probs[1] += calcProbOfQubitOutcome(tempQureg, qubit, 1);
        }
    } else {
        // Qubit is outside chunk size - need to handle chunk mapping and regions
        const auto& chunkMap = permTracker.getCurrentChunkMap();
        int qubitOffset = qubit - numQubitsPerChunk - 1;
        int regionSize = 1 << qubitOffset; // 2^(qubit - qubitsPerChunk - 1)
        
        for (size_t logicalChunk = 0; logicalChunk < numChunks; ++logicalChunk) {
            size_t physicalChunk = chunkMap[logicalChunk];
            std::vector<qcomp> buffer;
            loadChunk(physicalChunk, buffer);
            Qureg tempQureg = createTempQureg(buffer, numQubitsPerChunk);
            
            // Determine which region this chunk belongs to
            int region = logicalChunk / regionSize;
            int outcome = region % 2; // 0 or 1 based on region
            
            // Calculate probability for this chunk
            qreal chunkProb = calcTotalProb(tempQureg);
            probs[outcome] += chunkProb;
            
        }
    }
    
    // Determine measurement outcome
    int outcome = rand_getRandomSingleQubitOutcome(probs[0]);
    
    // Renormalize based on outcome
    qreal correctProb = probs[outcome];
    qreal wrongProb = probs[1 - outcome];
    qreal normalizationFactor = 1.0 / std::sqrt(correctProb);
    
    // Apply renormalization to all chunks
    
    for (size_t logicalChunk = 0; logicalChunk < numChunks; ++logicalChunk) {
        size_t physicalChunk = permTracker.getCurrentChunkMap()[logicalChunk];
        std::vector<qcomp> buffer;
        loadChunk(physicalChunk, buffer);
        
        if (qubit <= numQubitsPerChunk) {
            // For qubits within chunk size, check each amplitude
            for (size_t i = 0; i < buffer.size(); ++i) {
                // Check if this amplitude corresponds to the wrong outcome
                bool isWrongOutcome = false;
                if (outcome == 0) {
                    // Check if amplitude corresponds to qubit=1
                    if ((i >> qubit) & 1) {
                        isWrongOutcome = true;
                    }
                } else {
                    // Check if amplitude corresponds to qubit=0
                    if (!((i >> qubit) & 1)) {
                        isWrongOutcome = true;
                    }
                }
                
                if (isWrongOutcome) {
                    buffer[i] = 0.0; // Zero out wrong amplitudes
                } else {
                    buffer[i] *= normalizationFactor; // Renormalize correct amplitudes
                }
            }
        } else {
            // For qubits outside chunk size, check region
            int qubitOffset = qubit - numQubitsPerChunk - 1;
            int regionSize = 1 << qubitOffset;
            int region = logicalChunk / regionSize;
            int chunkOutcome = region % 2;
            
            if (chunkOutcome != outcome) {
                // Zero out entire chunk if it corresponds to wrong outcome
                std::fill(buffer.begin(), buffer.end(), 0.0);
            } else {
                // Renormalize entire chunk if it corresponds to correct outcome
                for (auto& amp : buffer) {
                    amp *= normalizationFactor;
                }
            }
        }
        
        saveChunk(physicalChunk, buffer);
    }
    
    return outcome;
}