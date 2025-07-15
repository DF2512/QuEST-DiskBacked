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

            if (verbose) {
                cout << "Swap index " << i << " and " << j << " => ";
                for (int v : current) cout << v << " ";
                cout << endl;
            }
        }
    }
    return steps;
}

// This functionn is the foundation of generating transitions, it creates the memory swap steps and provides the 'levels' to provide to the area swap shuffle
tuple<vector<pair<int, int>>, vector<pair<int, int>>, vector<int>>
PermutationTracker::alignUpperQubitsAndGenerateInterim(
    const vector<int>& previous_permutation,
    const vector<int>& target_permutation
) {
    vector<int> temp_interim = target_permutation;

    int upper_start = qubitsPerBlock;
    int upper_end = numQubits;
    int middle_start = qubitsPerChunk;
    int middle_end = qubitsPerBlock;

    vector<int> upper_qubits_prev(previous_permutation.begin() + upper_start, previous_permutation.begin() + upper_end);
    set<int> upper_qubits_set(upper_qubits_prev.begin(), upper_qubits_prev.end());

    vector<pair<int, int>> secondary_swaps;
    int insertion_index = middle_start;

    for (int q : upper_qubits_prev) {
        auto it = find(temp_interim.begin(), temp_interim.end(), q);
        int current_index = distance(temp_interim.begin(), it);

        if (insertion_index >= middle_end)
            throw runtime_error("Not enough room in middle segment to relocate qubits.");

        if (current_index >= middle_start && current_index < middle_end) {
            ++insertion_index;
            continue;
        }

        if (current_index != insertion_index) {
            swap(temp_interim[insertion_index], temp_interim[current_index]);
            secondary_swaps.emplace_back(insertion_index, current_index);
        }

        ++insertion_index;
    }

    vector<int> middle_qubits_now(temp_interim.begin() + middle_start, temp_interim.begin() + middle_end);
    set<int> middle_qubits_set(middle_qubits_now.begin(), middle_qubits_now.end());

    if (middle_qubits_set != upper_qubits_set)
        throw runtime_error("Failed to align all upper qubits in the middle range.");

    vector<int> interim_target = temp_interim;

    for (int i = 0; i < upper_end - upper_start; ++i) {
        swap(temp_interim[upper_start + i], temp_interim[middle_start + i]);
    }

    vector<pair<int, int>> initial_swaps = getSwapStepsToTarget(
        vector<int>(previous_permutation.begin(), previous_permutation.begin() + qubitsPerBlock),
        vector<int>(temp_interim.begin(), temp_interim.begin() + qubitsPerBlock)
    );

    reverse(secondary_swaps.begin(), secondary_swaps.end());

    return {initial_swaps, secondary_swaps, interim_target};
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
    std::cout << "[DEBUG] swapToMatch called\n";
    std::cout << "[DEBUG] current: ";
    for (int c : current) std::cout << c << " ";
    std::cout << "\n";
    std::cout << "[DEBUG] target: ";
    for (int t : target) std::cout << t << " ";
    std::cout << "\n";

    vector<int> steps;

    for (size_t i = 0; i < current.size(); ++i) {
        if (current[i] == target[i]) continue;

        auto it = find(current.begin() + i, current.end(), target[i]);
        int target_index = distance(current.begin(), it);

        std::cout << "[DEBUG] Moving element " << target[i] << " from position " << target_index << " to " << i << "\n";

        for (int j = target_index; j > (int)i; --j) {
            swap(current[j], current[j - 1]);
            steps.push_back(j);
            std::cout << "[DEBUG] Added swap level: " << j << "\n";
        }
    }

    std::cout << "[DEBUG] swapToMatch returning " << steps.size() << " levels\n";
    return steps;
}

// Function containing logic for chunk index permuting
vector<int> PermutationTracker::areaSwapShuffle(const vector<int>& perm, int level, int max_level) {
    std::cout << "[DEBUG] areaSwapShuffle called with:\n";
    std::cout << "[DEBUG] perm size: " << perm.size() << "\n";
    std::cout << "[DEBUG] level: " << level << "\n";
    std::cout << "[DEBUG] max_level: " << max_level << "\n";
    
    int x = 1 << (max_level - level);
    cout << "\n=== Shuffling with x = " << x << " ===\n";
    cout << "[Debug] Starting area swap shuffle...\n";

    if (perm.size() % x != 0)
        throw invalid_argument("Permutation length must be divisible by region size x.");

    int region_size = perm.size() / x;
    vector<int> new_perm;
    cout << "[Debug] Processing " << x << " regions of size " << region_size << "\n";

    for (int r = 0; r < x; ++r) {
        if (r % 10 == 0) cout << "[Debug] Processing region " << r << "/" << x << "\n";
        
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
    
    cout << "[Debug] Area swap shuffle completed\n";
    return new_perm;
}

// Function to generate a transition from previous to target permutation
Transition PermutationTracker::generateTransition(
    const std::vector<int>& prev_perm,
    const std::vector<int>& target_perm
) {
    std::cout << "[DEBUG] generateTransition called\n";
    std::cout << "[DEBUG] prev_perm: ";
    for (int p : prev_perm) std::cout << p << " ";
    std::cout << "\n";
    std::cout << "[DEBUG] target_perm: ";
    for (int p : target_perm) std::cout << p << " ";
    std::cout << "\n";

    auto [swap1, swap2_rev, interim_target] = alignUpperQubitsAndGenerateInterim(prev_perm, target_perm);

    std::cout << "[DEBUG] After alignUpperQubitsAndGenerateInterim:\n";
    std::cout << "[DEBUG] swap1 size: " << swap1.size() << "\n";
    std::cout << "[DEBUG] swap2_rev size: " << swap2_rev.size() << "\n";
    std::cout << "[DEBUG] interim_target: ";
    for (int p : interim_target) std::cout << p << " ";
    std::cout << "\n";

    std::vector<int> interim_prev = applySwaps(prev_perm, swap1);

    std::cout << "[DEBUG] interim_prev: ";
    for (int p : interim_prev) std::cout << p << " ";
    std::cout << "\n";

    std::vector<int> current(interim_prev.begin() + qubitsPerChunk, interim_prev.begin() + numQubits);
    std::vector<int> target(interim_target.begin() + qubitsPerChunk, interim_target.begin() + numQubits);

    std::cout << "[DEBUG] current (middle+upper): ";
    for (int p : current) std::cout << p << " ";
    std::cout << "\n";
    std::cout << "[DEBUG] target (middle+upper): ";
    for (int p : target) std::cout << p << " ";
    std::cout << "\n";

    std::cout << "[DEBUG] generate_transition: comparing current vs target permutation\n";
    for (size_t i = 0; i < current.size(); ++i) {
        if (current[i] != target[i]) {
            std::cout << "  mismatch at index " << i << ": current = " << current[i]
                      << ", target = " << target[i] << "\n";
        }
    }

    std::vector<int> levels = swapToMatch(current, target);

    std::cout << "[DEBUG] swapToMatch returned levels: ";
    for (int level : levels) std::cout << level << " ";
    std::cout << "\n";

    std::reverse(swap2_rev.begin(), swap2_rev.end());

    std::cout << "[DEBUG] Final transition:\n";
    std::cout << "[DEBUG] swap1 size: " << swap1.size() << "\n";
    std::cout << "[DEBUG] levels size: " << levels.size() << "\n";
    std::cout << "[DEBUG] swap2_rev size: " << swap2_rev.size() << "\n";

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

    for (size_t i = 0; i < subcircuits.size() - 1; ++i) {
        const auto& prev = subcircuits[i];
        const auto& next = subcircuits[i + 1];

        Transition t = generateTransition(prev.permutation, next.permutation);
        transitions.push_back(t);
    }

    return transitions;
}

// ─── Constructor ──────────────────────────────────────
PermutationTracker::PermutationTracker(int numQubits_, int numBlocks_, int chunksPerBlock_)
    : numQubits(numQubits_), numBlocks(numBlocks_), chunksPerBlock(chunksPerBlock_)
{
    numChunks = numBlocks * chunksPerBlock;
    qubitsPerBlock = numQubits - static_cast<int>(std::log2(numBlocks));
    qubitsPerChunk = qubitsPerBlock - static_cast<int>(std::log2(chunksPerBlock));
    
    // Initial layout: identity permutation and linear chunk map
    currentPermutation.resize(numQubits);
    std::iota(currentPermutation.begin(), currentPermutation.end(), 0);

    chunkMap.resize(numChunks);
    std::iota(chunkMap.begin(), chunkMap.end(), 0);

    std::cout << "[Init] Initial chunkMap:\n";
    for (int i = 0; i < chunkMap.size(); ++i) {
        std::cout << "  chunkMap[" << i << "] = " << chunkMap[i] << "\n";
    }
    std::cout << "[Init] Initial currentPermutation:\n";
    for (int i = 0; i < currentPermutation.size(); ++i) {
        std::cout << "  currentPermutation[" << i << "] = " << currentPermutation[i] << "\n";
    }
}



// ─── Block-Chunk Mapping ──────────────────────────────
std::vector<std::vector<int>> PermutationTracker::getBlockChunkMapping() const {
    std::cout << "[PermutationTracker] Generating block-chunk mapping...\n";
    std::cout << "  numBlocks: " << numBlocks << "\n";
    std::cout << "  chunksPerBlock: " << chunksPerBlock << "\n";
    std::cout << "  chunkMap.size(): " << chunkMap.size() << "\n";

    std::vector<std::vector<int>> blocks(numBlocks);

    for (int i = 0; i < numBlocks; ++i) {
        std::cout << "  Block " << i << ":\n";
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
            std::cout << "    chunkMap[" << flatIndex << "] = " << chunk << "\n";

            blocks[i].push_back(chunk);
        }
    }

    std::cout << "[PermutationTracker] Block-chunk mapping complete.\n";
    return blocks;
}

