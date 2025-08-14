#pragma once

#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <set>
#include <chunkmanager.h>
#include <diskbackedstate.h>
#include <gatescheduler.h>
#include <stdexcept>
#include <cmath>
#include <numeric>
#include <algorithm>

using namespace std;

// Helper function to get swap steps to target permutation
vector<pair<int, int>> getSwapStepsToTarget(
    const vector<int>& previous_permutation,
    const vector<int>& target_permutation,
    bool verbose = false
) {
    vector<int> current = previous_permutation;
    vector<pair<int, int>> steps;

    for (size_t i = 0; i < current.size(); ++i) {
        while (current[i] != target_permutation[i]) {
            int correct_val = target_permutation[i];
            auto it = find(current.begin(), current.end(), correct_val);
            int j = distance(current.begin(), it);

            swap(current[i], current[j]);
            steps.emplace_back(i, j);
        }
    }
    return steps;
}

// This functionn is the foundation of generating transitions, it creates the memory swap steps
tuple<vector<pair<int, int>>, vector<pair<int, int>>, vector<int>>
PermutationTracker::alignUpperQubitsAndGenerateInterim(
    const vector<int>& previous_permutation,
    const vector<int>& target_permutation
) {
    // Define region boundaries
    int lower_start = 0;
    int lower_end = qubitsPerChunk;
    int middle_start = qubitsPerChunk;
    int middle_end = qubitsPerBlock;
    int upper_start = qubitsPerBlock;
    int upper_end = numQubits;

    // Make interim_prev and interim_target
    vector<int> interim_prev = previous_permutation;
    vector<int> interim_target = target_permutation;
    vector<pair<int, int>> swap1;
    vector<pair<int, int>> swap2_rev;

    // For swap1: For each qubit q in the upper region of the target permutation
    std::vector<bool> middle_used(middle_end - middle_start, false);
    std::vector<bool> upper_preserved(upper_end - upper_start, false);

    // First pass: Mark qubits that are already in their correct positions in the upper region
    for (int qidx = 0; qidx < upper_end - upper_start; ++qidx) {
        int q = target_permutation[upper_start + qidx];
        int target_pos = upper_start + qidx;

        // Check if q is already at its correct position in interim_prev
        if (interim_prev[target_pos] == q) {
            upper_preserved[qidx] = true;
           // std::cout << "[DEBUG] Preserving qubit " << q << " at position " << target_pos << " (already correct)" << std::endl;


        }
    }

    // Second pass: Process qubits that need to be moved, but skip preserved ones
    for (int qidx = 0; qidx < upper_end - upper_start; ++qidx) {
        int q = target_permutation[upper_start + qidx];

        // Skip if this qubit is already preserved in its correct position
        if (upper_preserved[qidx]) {
            continue;
        }

        // Recheck all previous qs for correctness
        for (int prev_qidx = 0; prev_qidx < qidx; ++prev_qidx) {
            int prev_q = target_permutation[upper_start + prev_qidx];
            auto it_up = std::find(interim_prev.begin() + upper_start, interim_prev.begin() + upper_end, prev_q);
           // if (it_up == interim_prev.begin() + upper_end) {
           //     std::cerr << "[alignUpperQubitsAndGenerateInterim] ERROR: prev_q " << prev_q << " not in upper region after supposed placement!" << std::endl;
           // }
        }
        // If q is present in the upper region of interim_prev, do nothing
        auto it_up = std::find(interim_prev.begin() + upper_start, interim_prev.begin() + upper_end, q);
        if (it_up != interim_prev.begin() + upper_end) {
            continue;
        }
        // Else if q is present in the middle region, mark that slot as used
        auto it_mid = std::find(interim_prev.begin() + middle_start, interim_prev.begin() + middle_end, q);
        if (it_mid != interim_prev.begin() + middle_end) {
            int mid_idx = std::distance(interim_prev.begin() + middle_start, it_mid);
            middle_used[mid_idx] = true;
            continue;
        }
        // Else, search for q in interim_prev and swap it into the first available (unused) middle region slot
        auto it = std::find(interim_prev.begin() + lower_start, interim_prev.begin() + middle_start, q);
        if (it == interim_prev.begin() + middle_start) {
            std::cerr << "[alignUpperQubitsAndGenerateInterim] ERROR: q " << q << " not found in lower or middle region!" << std::endl;
            continue;
        }
        // Find first available (unused) middle slot
        int middle_slot = -1;
        for (int m = 0; m < (middle_end - middle_start); ++m) {
            if (!middle_used[m]) {
                middle_slot = middle_start + m;
                middle_used[m] = true;
                break;
            }
        }
        if (middle_slot == -1) {
            std::cerr << "[alignUpperQubitsAndGenerateInterim] ERROR: No available middle slot for swap1!" << std::endl;
            continue;
        }
        int lower_idx = std::distance(interim_prev.begin(), it);
        std::swap(interim_prev[lower_idx], interim_prev[middle_slot]);
        swap1.emplace_back(lower_idx, middle_slot);

        // After each q, recheck all previous qs
        for (int prev_qidx = 0; prev_qidx <= qidx; ++prev_qidx) {
            int prev_q = target_permutation[upper_start + prev_qidx];
            auto it_up2 = std::find(interim_prev.begin() + upper_start, interim_prev.begin() + upper_end, prev_q);
           // if (it_up2 == interim_prev.begin() + upper_end) {
           //     std::cerr << "[alignUpperQubitsAndGenerateInterim] ERROR: prev_q " << prev_q << " not in upper region after swap!" << std::endl;
           // }
        }
    }

    // For swap2: Make interim_target's lower region match interim_prev's lower region
    while (true) {
        bool all_match = true;
        for (int i = lower_start; i < lower_end; ++i) {
            if (interim_target[i] == interim_prev[i]) continue;
            all_match = false;
            // Find where interim_prev[i] is in lower+middle region of interim_target
            auto it = std::find(interim_target.begin() + i + 1, interim_target.begin() + middle_end, interim_prev[i]);
            if (it == interim_target.begin() + middle_end) {
                std::cerr << "[alignUpperQubitsAndGenerateInterim] ERROR: value " << interim_prev[i] << " not found in lower+middle region of interim_target!" << std::endl;
                continue;
            }
            int j = std::distance(interim_target.begin(), it);
            std::swap(interim_target[i], interim_target[j]);
            swap2_rev.emplace_back(i, j);

            // After a swap, restart from the beginning
            break;
        }
        if (all_match) break;
    }

    // Check that the lower regions of interim_prev and interim_target match
    bool lower_match = true;
    for (int i = lower_start; i < lower_end; ++i) {
        if (interim_prev[i] != interim_target[i]) {
            lower_match = false;
            break;
        }
    }
    if (!lower_match) {
        std::cerr << "[alignUpperQubitsAndGenerateInterim] WARNING: Lower regions of interim_prev and interim_target do not match!\n";
        std::cerr << "  interim_prev lower: ";
        for (int i = lower_start; i < lower_end; ++i) std::cerr << interim_prev[i] << " ";
        std::cerr << "\n  interim_target lower: ";
        for (int i = lower_start; i < lower_end; ++i) std::cerr << interim_target[i] << " ";
        std::cerr << std::endl;
    }

    // Check that the sets of the combined middle+upper regions match
    std::set<int> interim_prev_middle_upper(interim_prev.begin() + middle_start, interim_prev.end());
    std::set<int> interim_target_middle_upper(interim_target.begin() + middle_start, interim_target.end());
    if (interim_prev_middle_upper != interim_target_middle_upper) {
        std::cerr << "[alignUpperQubitsAndGenerateInterim] WARNING: Combined middle+upper region sets do not match!\n";
        std::cerr << "  interim_prev middle+upper: ";
        for (auto v : interim_prev_middle_upper) std::cerr << v << " ";
        std::cerr << "\n  interim_target middle+upper: ";
        for (auto v : interim_target_middle_upper) std::cerr << v << " ";
        std::cerr << std::endl;
    }

    // Return swap1, swap2_rev, interim_target
    return {swap1, swap2_rev, interim_target};
}

// Helper function to apply swaps to a permutation
vector<int> applySwaps(const vector<int>& permutation, const vector<pair<int, int>>& swaps) {
    vector<int> perm = permutation;
    for (auto [i, j] : swaps) {
        swap(perm[i], perm[j]);
    }
    return perm;
}

// Function to generate swap steps to match a target permutation
vector<int> swapToMatch(vector<int> current, const vector<int>& target) {


    vector<int> steps;

    for (size_t i = 0; i < current.size(); ++i) {
        if (current[i] == target[i]) continue;

        auto it = find(current.begin() + i, current.end(), target[i]);
        int target_index = distance(current.begin(), it);

        //std::cout << "[DEBUG] Moving element " << target[i] << " from position " << target_index << " to " << i << "\n";

        for (int j = target_index; j > (int)i; --j) {
            swap(current[j], current[j - 1]);
            steps.push_back(j);
            //std::cout << "[DEBUG] Added swap level: " << j << "\n";
        }
    }

    //std::cout << "[DEBUG] swapToMatch returning " << steps.size() << " levels\n";
    return steps;
}

// Function containing logic for chunk index permuting
vector<int> PermutationTracker::areaSwapShuffle(const vector<int>& perm, int level, int max_level) {

    int x = 1 << (max_level - level);
    //cout << "\n=== Shuffling with x = " << x << " ===\n";
    //cout << "[Debug] Starting area swap shuffle...\n";

    if (perm.size() % x != 0)
        throw invalid_argument("Permutation length must be divisible by region size x.");

    int region_size = perm.size() / x;
    vector<int> new_perm;
    //cout << "[Debug] Processing " << x << " regions of size " << region_size << "\n";

    for (int r = 0; r < x; ++r) {
        //if (r % 10 == 0) cout << "[Debug] Processing region " << r << "/" << x << "\n";

        int start = r * region_size;
        int end = start + region_size;
        vector<int> region(perm.begin() + start, perm.begin() + end);

        // If region_size is too small, just copy the region as-is
        if (region_size < 4) {
            new_perm.insert(new_perm.end(), region.begin(), region.end());
            continue;
        }

        if (region_size % 4 != 0)
            throw invalid_argument("Region size must be divisible by 4.");

        int chunk_size = region_size / 4;

        vector<int> area0(region.begin(), region.begin() + chunk_size);
        vector<int> area1(region.begin() + chunk_size, region.begin() + 2 * chunk_size);
        vector<int> area2(region.begin() + 2 * chunk_size, region.begin() + 3 * chunk_size);
        vector<int> area3(region.begin() + 3 * chunk_size, region.end());

        new_perm.insert(new_perm.end(), area0.begin(), area0.end());
        new_perm.insert(new_perm.end(), area2.begin(), area2.end());
        new_perm.insert(new_perm.end(), area1.begin(), area1.end());
        new_perm.insert(new_perm.end(), area3.begin(), area3.end());
    }

    //cout << "[Debug] Area swap shuffle completed\n";
    return new_perm;
}

// Function to generate a transition from previous to target permutation
Transition PermutationTracker::generateTransition(
    const std::vector<int>& prev_perm,
    const std::vector<int>& target_perm
) {

    auto [swap1, swap2_rev, interim_target] = alignUpperQubitsAndGenerateInterim(prev_perm, target_perm);

    std::vector<int> interim_prev = applySwaps(prev_perm, swap1);

    std::vector<int> current(interim_prev.begin() + qubitsPerChunk, interim_prev.begin() + numQubits);
    std::vector<int> target(interim_target.begin() + qubitsPerChunk, interim_target.begin() + numQubits);

    std::vector<int> levels = swapToMatch(current, target);

    std::reverse(swap2_rev.begin(), swap2_rev.end());

    return Transition{
        swap1,
        levels,
        swap2_rev,
        interim_target
    };
}

//Function to generate all transitions for a set of subcircuits
std::vector<Transition> PermutationTracker::generateTransitions(const std::vector<SubCircuit>& subcircuits) {
    std::vector<Transition> transitions;

    // Generate transitions for all pairs, including the last one
    for (size_t i = 0; i < subcircuits.size() - 1; ++i) {
        const auto& prev = subcircuits[i];
        const auto& next = subcircuits[i + 1];

        Transition t = generateTransition(prev.permutation, next.permutation);
        transitions.push_back(t);
    }

    return transitions;
}

// Constructor for PermutationTracker
PermutationTracker::PermutationTracker(int numQubits_, int numBlocks_, int chunksPerBlock_)
    : numQubits(numQubits_), numBlocks(numBlocks_), chunksPerBlock(chunksPerBlock_)
{
    numChunks = numBlocks * chunksPerBlock;
    qubitsPerBlock = numQubits - static_cast<int>(std::log2(numBlocks));
    qubitsPerChunk = qubitsPerBlock - static_cast<int>(std::log2(chunksPerBlock));

    // identity permutation and linear chunk map
    currentPermutation.resize(numQubits);
    std::iota(currentPermutation.begin(), currentPermutation.end(), 0);

    chunkMap.resize(numChunks);
    std::iota(chunkMap.begin(), chunkMap.end(), 0);
}



// Function to determine the current block-chunk mapping
std::vector<std::vector<int>> PermutationTracker::getBlockChunkMapping() const {


    std::vector<std::vector<int>> blocks(numBlocks);

    for (int i = 0; i < numBlocks; ++i) {
        //std::cout << "  Block " << i << ":\n";
        for (int j = 0; j < chunksPerBlock; ++j) {
            int flatIndex = i * chunksPerBlock + j;

            // Check for out-of-bounds
            if (flatIndex >= static_cast<int>(chunkMap.size())) {
                std::cerr << "  [ERROR] flatIndex " << flatIndex
                          << " out of bounds for chunkMap (size: "
                          << chunkMap.size() << ")\n";
                continue;  // Or throw std::runtime_error if fatal
            }

            int chunk = chunkMap[flatIndex];
            //std::cout << "    chunkMap[" << flatIndex << "] = " << chunk << "\n";

            blocks[i].push_back(chunk);
        }
    }

    //std::cout << "[PermutationTracker] Block-chunk mapping complete.\n";
    return blocks;
}

//Accessors
const std::vector<int>& PermutationTracker::getCurrentPermutation() const {
    return currentPermutation;
}

const std::vector<int>& PermutationTracker::getCurrentChunkMap() const {
    return chunkMap;
}
