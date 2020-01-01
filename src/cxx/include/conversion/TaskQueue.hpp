/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file TaskQueue.hpp
 *
 * \brief Declaration of task queue.
 *
 * \author Aaron Day
 */


#ifndef BCL2FASTQ_CONVERSION_TASKQUEUE_HPP
#define BCL2FASTQ_CONVERSION_TASKQUEUE_HPP


#include "conversion/TaskQueue.hh"
#include <numeric>

namespace bcl2fastq {
namespace conversion {

template<typename T, uint8_t NumPriorityLevels>
size_t PriorityQueue<T, NumPriorityLevels>::size()
{
    return std::accumulate(queue_.begin(), queue_.end(), 0, [](size_t total, std::deque<T>& queue) { return total + queue.size(); });
}

template<typename T, uint8_t NumPriorityLevels>
bool PriorityQueue<T, NumPriorityLevels>::empty()
{
    for (auto& queue : queue_)
    {
        if (!queue.empty()) return false;
    }

    return true;
}

template<typename T, uint8_t NumPriorityLevels>
T& PriorityQueue<T, NumPriorityLevels>::front()
{
    for (auto& queue : queue_)
    {
        if (!queue.empty()) return queue.front();
    }
}

template<typename T, uint8_t NumPriorityLevels>
void PriorityQueue<T, NumPriorityLevels>::push_back(T& data)
{
    queue_.at(data->getPriority()).push_back(data);
}

template<typename T, uint8_t NumPriorityLevels>
void PriorityQueue<T, NumPriorityLevels>::pop_front()
{
    for (auto& queue : queue_)
    {
        if (!queue.empty())
        {
            queue.pop_front();
            return;
        }
    }
}

template<typename T, uint8_t NumPriorityLevels>
void PriorityQueue<T, NumPriorityLevels>::clear()
{
    for (auto& queue : queue_)
    {
        queue.clear();
    }
}

typedef ThreadSafeQueue<TaskPtr, PriorityQueue<TaskPtr, Task::PriorityLevel::NUM_PRIORITY_LEVELS>> TaskQueue;

} // namespace conversion
} // namespace bcl2fastq


#endif // BCL2FASTQ_CONVERSION_TASKQUEUE_HPP


