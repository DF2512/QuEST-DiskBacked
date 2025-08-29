#pragma once
#include <vector>
#include <utility>
#include "quest/include/gatescheduler.h"

// Forward declaration
class DiskBackedState;

// ─── Chunk Map Generator ────────────────────────────────

// ─── Transition Object ──────────────────────────────────
struct Transition {
    std::vector<std::pair<int, int>> swap1;        // prev → interim
    std::vector<int> swapLevels;                   // area swap levels
    std::vector<std::pair<int, int>> swap2;        // interim → target
    std::vector<int> interimTarget;                // full layout after swap1
};

// ─── Permutation Tracker ────────────────────────────────
class PermutationTracker {
public:
    PermutationTracker(int numQubits, int numBlocks, int chunksPerBlock);

    void updateToPermutation(const std::vector<int>& targetPerm);
    const std::vector<int>& getCurrentPermutation() const;
    const std::vector<int>& getCurrentChunkMap() const;
    std::vector<std::vector<int>> getBlockChunkMapping() const;
    void applyTransition(const Transition& t);
    std::vector<int> chunkMap;
    std::tuple<std::vector<std::pair<int, int>>, std::vector<std::pair<int, int>>, std::vector<int>>
    alignUpperQubitsAndGenerateInterim(
        const std::vector<int>& previous_permutation,
        const std::vector<int>& target_permutation
    );
    Transition generateTransition(
    const std::vector<int>& prev_perm,
    const std::vector<int>& target_perm
    );
    std::vector<int> areaSwapShuffle(const std::vector<int>& perm, int level, int max_level);
    void print_blocks(const std::vector<int>& chunk_map);
    std::vector<Transition> generateTransitions(const std::vector<SubCircuit>& subcircuits);
    std::vector<int> currentPermutation;


private:

    int numQubits;
    int numBlocks;
    int chunksPerBlock;
    int qubitsPerBlock;
    int qubitsPerChunk;
    int numChunks;
    void computeChunkShuffle(const std::vector<int>& targetPerm);

};

// Declaration for applySwaps
std::vector<int> applySwaps(const std::vector<int>& permutation, const std::vector<std::pair<int, int>>& swaps);
