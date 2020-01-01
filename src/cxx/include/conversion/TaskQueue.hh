/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file TaskQueue.hh
 *
 * \brief Declaration of task queue.
 *
 * \author Aaron Day
 */


#ifndef BCL2FASTQ_CONVERSION_TASKQUEUE_HH
#define BCL2FASTQ_CONVERSION_TASKQUEUE_HH


#include "conversion/Task.hh"
#include "conversion/ThreadSafeQueue.hh"
#include <boost/noncopyable.hpp>
#include <array>
#include <deque>

namespace bcl2fastq {
namespace conversion {

template<typename T, uint8_t NumPriorityLevels>
class PriorityQueue : private boost::noncopyable
{
public:
    PriorityQueue() : queue_() { }

    size_t size();
    bool empty();
    T& front();

    void push_back(T& data);
    void pop_front();

    void clear();

private:
    std::array<typename std::deque<T>, NumPriorityLevels> queue_;
};

} // namespace conversion
} // namespace bcl2fastq

#include "conversion/TaskQueue.hpp"

#endif // BCL2FASTQ_CONVERSION_TASKQUEUE_HH


