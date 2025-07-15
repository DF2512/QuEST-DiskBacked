#include "diskbackedstate.h"
#include <filesystem>
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <random>
#include <omp.h>

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

/*
const PermutationTracker& DiskBackedState::getPermutationTracker() const {
    return permTracker;
}
*/

// Accessors
int DiskBackedState::getChunksPerBlock() const { return chunksPerBlock; }
int DiskBackedState::getMaxPermutableQubits() const {return maxPermutableQubits; }
int DiskBackedState::getNumQubitsPerBlock() const { return numQubitsPerBlock; }
int DiskBackedState::getNumBlocks() const { return numBlocks; }
int DiskBackedState::getNumQubits() const { return numQubits; }
size_t DiskBackedState::getNumAmplitudes() const { return numAmplitudes; }
size_t DiskBackedState::getNumChunks() const { return numChunks; }
size_t DiskBackedState::getAmpsPerChunk() const { return ampsPerChunk; }

void DiskBackedState::loadChunk(size_t chunkIndex, std::vector<qcomp>& buffer) const {
    if (chunkIndex >= numChunks) {
        throw std::runtime_error("loadChunk: chunkIndex out of bounds");
    }

    buffer.resize(ampsPerChunk);
    const std::string& path = chunkPaths[chunkIndex];

    FILE* f = fopen(path.c_str(), "rb");
    if (!f) {
        throw std::runtime_error("loadChunk: failed to open file " + path);
    }

    size_t read = fread(buffer.data(), sizeof(qcomp), ampsPerChunk, f);
    fclose(f);

    if (read != ampsPerChunk) {
        throw std::runtime_error("loadChunk: read incomplete");
    }
}

void DiskBackedState::saveChunk(size_t chunkIndex, const std::vector<qcomp>& buffer) const {
    if (chunkIndex >= numChunks) {
        throw std::runtime_error("saveChunk: chunkIndex out of bounds");
    }

    if (buffer.size() != ampsPerChunk) {
        throw std::runtime_error("saveChunk: buffer size does not match ampsPerChunk");
    }

    const std::string& path = chunkPaths[chunkIndex];

    FILE* f = fopen(path.c_str(), "wb");
    if (!f) {
        throw std::runtime_error("saveChunk: failed to open file " + path);
    }

    size_t written = fwrite(buffer.data(), sizeof(qcomp), ampsPerChunk, f);
    fflush(f);
    fclose(f);

    if (written != ampsPerChunk) {
        throw std::runtime_error("saveChunk: write incomplete");
    }
}

void DiskBackedState::initialiseRandomState() {
    std::cout << "[Init] Generating random amplitudes...\n";

    double totalNormSq = 0.0;
    std::vector<qcomp> buffer(ampsPerChunk);
    std::mt19937_64 rng(42); // fixed seed for reproducibility
    std::uniform_real_distribution<qreal> dist(-0.5, 0.5);

    // ─── Pass 1: Generate and write random complex numbers ─────────
    for (size_t chunk = 0; chunk < numChunks; ++chunk) {
        for (size_t i = 0; i < ampsPerChunk; ++i) {
            qreal re = dist(rng);
            qreal im = dist(rng);
            buffer[i] = qcomp(re, im);
            totalNormSq += std::norm(buffer[i]);
        }

        if (chunk == 0) {
            std::cout << "[Debug] First amplitude after generation: "
                      << buffer[0] << "\n";
        }

        saveChunk(chunk, buffer);
    }

    std::cout << "[Init] Total norm² = " << totalNormSq << "\n";

    // ─── Compute scaling factor ─────────────
    qreal scale = 1.0 / std::sqrt(totalNormSq);
    std::cout << "[Init] Applying scale factor: " << scale << "\n";

    // ─── Pass 2: Normalize and save ─────────
    for (int chunk = 0; chunk < static_cast<int>(numChunks); ++chunk) {
        std::vector<qcomp> buf;
        loadChunk(chunk, buf);

        #pragma omp parallel for schedule(static)
        for (int i = 0; i < static_cast<int>(buf.size()); ++i) {
            buf[i] *= scale;
        }

        if (chunk == 0) {
            std::cout << "[Debug] First amplitude after normalization: "
                      << buf[0] << "\n";
        }

        saveChunk(chunk, buf);
        buf.clear();
        buf.shrink_to_fit();
    }

    std::cout << "[Init] Normalization complete.\n";
}

double DiskBackedState::computeTotalProbability() const {
    std::cout << "[Measure] Computing total probability...\n";

    double total = 0.0;

    #pragma omp parallel for reduction(+:total)
    for (int chunk = 0; chunk < static_cast<int>(numChunks); ++chunk) {
        std::vector<qcomp> buf;
        loadChunk(chunk, buf);

        if (chunk == 0) {
            std::cout << "[Debug] First amplitude during probability check: "
                      << buf[0] << "\n";
        }

        double localSum = 0.0;
        for (const qcomp& amp : buf)
            localSum += std::norm(amp);

        total += localSum;
        buf.clear();
        buf.shrink_to_fit();
    }

    std::cout << "[Measure] Total probability = " << total << "\n";
    return total;
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

    std::cout << "[Disk] Loaded block " << blockIdx << " (chunks: ";
    for (int c : chunkIndices) std::cout << c << " ";
    std::cout << ")\n";
}

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

    std::cout << "[Disk] Saved block " << blockIdx << " (chunks: ";
    for (int c : chunkIndices) std::cout << c << " ";
    std::cout << ")\n";
}

DiskBackedState::~DiskBackedState() {
    deleteAllChunkFiles();
}

void DiskBackedState::deleteAllChunkFiles() {
    for (const auto& file : chunkPaths) {
        std::remove(file.c_str());
    }
}