#pragma once

#include "diskbackedstate.h"
#include "gatescheduler.h"

#include <thread>
#include <queue>
#include <mutex>
#include <condition_variable>

std::vector<GateOp> scheduleSwaps(const std::vector<std::pair<int, int>>& swaps);

void runCircuit(GateScheduler& scheduler, DiskBackedState& state, bool verbose = true);

void applySubCircuitToBlock(const SubCircuit& sub, std::vector<qcomp>& buffer, int qubitsPerBlock);


template <typename T>
class ThreadSafeQueue {
    std::queue<T> queue;
    std::mutex mtx;
    std::condition_variable cv_not_empty;
    std::condition_variable cv_not_full;
    bool finished = false;
    size_t maxSize;

public:
    ThreadSafeQueue(size_t maxSize_) : maxSize(maxSize_) {}

    void push(T item) {
        std::unique_lock<std::mutex> lock(mtx);
        cv_not_full.wait(lock, [&]() { return queue.size() < maxSize; });
        queue.push(std::move(item));
        cv_not_empty.notify_one();
    }

    bool pop(T& item) {
        std::unique_lock<std::mutex> lock(mtx);
        cv_not_empty.wait(lock, [&]() { return !queue.empty() || finished; });
        if (queue.empty()) return false;
        item = std::move(queue.front());
        queue.pop();
        cv_not_full.notify_one();
        return true;
    }

    void setFinished() {
        {
            std::lock_guard<std::mutex> lock(mtx);
            finished = true;
        }
        cv_not_empty.notify_all();
    }
};
struct BlockData {
    int blockIdx;
    std::vector<int> chunkIndices;
    int bufferIndex = -1;
};