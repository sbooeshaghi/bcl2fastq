/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file ThreadSafeQueue.hpp
 *
 * \brief Implementation of thread safe queue queue.
 *
 * \author Aaron Day
 */

#ifndef BCL2FASTQ_CONVERSION_THREAD_SAFE_QUEUE_HPP
#define BCL2FASTQ_CONVERSION_THREAD_SAFE_QUEUE_HPP

#include "common/Debug.hh"
#include "conversion/ThreadSafeQueue.hh"


namespace bcl2fastq {
namespace conversion {


template<typename T, class Queue>
ThreadSafeQueue<T, Queue>::ThreadSafeQueue(uint32_t maxQueueDepth /*= std::numceric_limits<uint32_t>::max()*/)
: maxQueueDepth_(maxQueueDepth)
, done_(false)
, terminated_(false)
, mut_()
, cvEmptyQueue_()
, cvFullQueue_()
, queue_()
{
}

template<typename T, class Queue>
void ThreadSafeQueue<T, Queue>::addData(T data, const std::string &description)
{
    BCL2FASTQ_LOG(common::LogLevel::TRACE) << description << std::endl;
    return addData(data);
}

template<typename T, class Queue>
void ThreadSafeQueue<T, Queue>::addData(T data)
{
    {
        std::unique_lock<std::mutex> lock(mut_);
        if (terminated_)
        {
            return;
        }

        if (queue_.size() >= maxQueueDepth_)
        {
            cvFullQueue_.wait(lock, [this] { return queue_.size() < maxQueueDepth_ || terminated_; });
        }

        if (terminated_)
        {
            return;
        }

        queue_.push_back(data);
    }

    cvEmptyQueue_.notify_one();
}

template<typename T, class Queue>
bool ThreadSafeQueue<T, Queue>::tryGetData(T& data, const std::string &description)
{
    BCL2FASTQ_LOG(common::LogLevel::TRACE) << description << std::endl;
    return tryGetData(data);
}

template<typename T, class Queue>
bool ThreadSafeQueue<T, Queue>::tryGetData(T& data)
{
    {
        std::unique_lock<std::mutex> lock(mut_);
        cvEmptyQueue_.wait(lock, [this] { return !queue_.empty() || done_ || terminated_; });

        if (queue_.empty() || terminated_) return false;

        data = queue_.front();
        queue_.pop_front();
    }

    cvFullQueue_.notify_one();

    return true;
}

template<typename T, class Queue>
bool ThreadSafeQueue<T, Queue>::tryGetDataCheckEmpty(T& data, const std::string &description)
{
    BCL2FASTQ_LOG(common::LogLevel::TRACE) << description << std::endl;
    return tryGetDataCheckEmpty(data);
}

template<typename T, class Queue>
bool ThreadSafeQueue<T, Queue>::tryGetDataCheckEmpty(T& data)
{
    {
        std::unique_lock<std::mutex> lock(mut_);

        if (queue_.empty()) return false;

        data = queue_.front();
        queue_.pop_front();
    }

    cvFullQueue_.notify_one();

    return true;
}

template<typename T, class Queue>
void ThreadSafeQueue<T, Queue>::setFinished()
{
    done_ = true;
    cvEmptyQueue_.notify_all();
}

template<typename T, class Queue>
void ThreadSafeQueue<T, Queue>::terminate()
{
    {
        std::unique_lock<std::mutex> lock(mut_);

        queue_.clear();
        terminated_ = true;
    }

    cvEmptyQueue_.notify_all();
    cvFullQueue_.notify_all();
}

template<typename KeyType, typename DataType>
ThreadSafeMap<KeyType, DataType>::ThreadSafeMap()
: lastGroupNumber_(std::numeric_limits<uint32_t>::max())
, done_(false)
, terminated_(false)
, mut_()
, cvEmpty_()
, map_()
{
}

template<typename KeyType, typename DataType>
void ThreadSafeMap<KeyType, DataType>::addData(const KeyType& key, DataType& data)
{
    {
        std::lock_guard<std::mutex> lock(mut_);
        if (!terminated_)
        {
            map_.insert(std::make_pair(key, data));
        }
    }

    cvEmpty_.notify_all();
}

template<typename KeyType, typename DataType>
bool ThreadSafeMap<KeyType, DataType>::tryGetData(const KeyType& key, DataType& data)
{
    {
        std::unique_lock<std::mutex> lock(mut_);

        typename std::map<KeyType, DataType>::iterator found = map_.end();

        cvEmpty_.wait(lock, [this, key, &found] { found = map_.find(key); return done_ || terminated_ || (found != map_.end() || key > lastGroupNumber_); });

        if (found == map_.end() || terminated_)
        {
            return false;
        }

        data = found->second;
        map_.erase(found);
    }

    cvEmpty_.notify_all();

    return true;
}

template<typename KeyType, typename DataType>
void ThreadSafeMap<KeyType, DataType>::setFinished(uint32_t lastGroupNumber)
{
    lastGroupNumber_ = lastGroupNumber;
    done_ = true;
    cvEmpty_.notify_all();
}

template<typename KeyType, typename DataType>
void ThreadSafeMap<KeyType, DataType>::terminate()
{
    {
        std::unique_lock<std::mutex> lock(mut_);

        map_.clear();
        terminated_ = true;
    }

    cvEmpty_.notify_all();
}

} // namespace conversion
} // namespace bcl2fastq

#endif // BCL2FASTQ_CONVERSION_THREAD_SAFE_QUEUE_HPP

