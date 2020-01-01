/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file ThreadSafeQueue.hh
 *
 * \brief Declaration of task queue.
 *
 * \author Aaron Day
 */


#ifndef BCL2FASTQ_CONVERSION_THREAD_SAFE_QUEUE_HH
#define BCL2FASTQ_CONVERSION_THREAD_SAFE_QUEUE_HH


#include <boost/noncopyable.hpp>
#include <deque>
#include <atomic>
#include <numeric>
#include <map>

namespace bcl2fastq {
namespace conversion {


/// \brief Thread safe queue.
template<typename T, class Queue=std::deque<T>>
class ThreadSafeQueue : boost::noncopyable
{
public:

    /// \brief Default constructor.
    ThreadSafeQueue(uint32_t maxQueueDepth = std::numeric_limits<uint32_t>::max());

    /// \brief Add task.
    /// \param task Task to be added.
    /// \note Task queue assumes ownership of the task instance.
    void addData(T data);
    void addData(T data, const std::string &description);

    /// \brief Get task.
    /// \return Task at the top of the queue.
    /// \note Ownership of the task is retained by the task queue throughout
    /// whole task's lifetime.
    bool tryGetData(T& data);
    bool tryGetData(T& data, const std::string &description);

    bool tryGetDataCheckEmpty(T& data);
    bool tryGetDataCheckEmpty(T& data, const std::string &description);

    bool isTerminated() const { return terminated_; }

    void setFinished();
    void terminate();

private:
    const uint32_t maxQueueDepth_;
    std::atomic<bool> done_;
    std::atomic<bool> terminated_;
    mutable std::mutex mut_;
    std::condition_variable cvEmptyQueue_;
    std::condition_variable cvFullQueue_;
    Queue queue_;
};

template<typename KeyType, typename DataType>
class ThreadSafeMap : boost::noncopyable
{
public:
    ThreadSafeMap();

    void addData(const KeyType& key, DataType& data);

    bool tryGetData(const KeyType& key, DataType& data);

    void setFinished(uint32_t lastGroupNumber);
    void terminate();

private:
    std::atomic<uint32_t> lastGroupNumber_;
    std::atomic<bool> done_;
    std::atomic<bool> terminated_;
    mutable std::mutex mut_;
    std::condition_variable cvEmpty_;
    std::map<KeyType, DataType> map_;    
};

} // namespace conversion
} // namespace bcl2fastq


#include "conversion/ThreadSafeQueue.hpp"

#endif // BCL2FASTQ_CONVERSION_THREAD_SAFE_QUEUE_HH


