/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file Threads.hpp
 *
 * \brief Implementation of helpers for thread management.
 *
 * \author Roman Petorvski
 * \author Marek Balint
 */


#ifndef BCL2FASTQ_COMMON_THREADS_HPP
#define BCL2FASTQ_COMMON_THREADS_HPP


#include <algorithm>

#include <boost/bind.hpp>
#include <boost/thread/locks.hpp>

#include "common/SystemCompatibility.hh"
#include "common/Logger.hh"
#include "common/Debug.hh"


namespace bcl2fastq {


namespace common {


template<typename Mutex>
unlock_guard<Mutex>::unlock_guard(Mutex& m)
: m_(m)
{
    m_.unlock();
}

template<typename Mutex>
unlock_guard<Mutex>::~unlock_guard()
{
    m_.lock();
}


template<bool CrashOnException>
BasicThreadVector<CrashOnException>::BasicThreadVector(typename BasicThreadVector<CrashOnException>::size_type size)
: pool_()
, mutex_()
, stateChangedCondition_()
, executor_(NULL)
, busyThreads_(size)
, neededThreads_(0)
, lowestBlockedThreadNumber_(0)
, terminateRequested_(false)
, currentRequest_(0)
, firstThreadException_()
{
    BCL2FASTQ_ASSERT_MSG(size > 0, "Inadequate pool size: " << size);
    pool_.reserve(size);

    boost::unique_lock<boost::mutex> lock(mutex_);
    while (size--)
    {
        pool_.push_back(new std::thread(
            &BasicThreadVector<CrashOnException>::threadFunction, this, size)
        );
    }
    // Initial wait for all threads to initialize their 'processedRequest'
    this->waitAll(lock);
}

template<bool CrashOnException>
BasicThreadVector<CrashOnException>::~BasicThreadVector()
{
    {
        /// \todo Refactoring: Locking is not really needed here, but it helps the assert to be a bit more definite.
        boost::unique_lock<boost::mutex> lock(mutex_);
        BCL2FASTQ_ASSERT_MSG(!busyThreads_, "Workers must not be running at this point");

        terminateRequested_ = true;
        stateChangedCondition_.notify_all();
    }

    std::for_each(
        pool_.begin(),
        pool_.end(),
        boost::bind(&std::thread::join, _1)
    );
}

template<bool CrashOnException>
typename BasicThreadVector<CrashOnException>::size_type BasicThreadVector<CrashOnException>::size() const
{
    return pool_.size();
}

template<bool CrashOnException>
template<typename FunctionT>
void BasicThreadVector<CrashOnException>::execute(FunctionT function, size_type threads)
{
    BCL2FASTQ_ASSERT_MSG(threads <= pool_.size(), "Request must not exceed the amount of threads available");
    BCL2FASTQ_ASSERT_MSG(executor_ == NULL, "Queueing is not supported");

    class FunctionExecutor : public Executor
    {
    public:

        FunctionExecutor(FunctionT &function)
        : function_(function)
        {
        }

    public:

        virtual void execute(BasicThreadVector<CrashOnException>::size_type threadNum)
        {
            function_(threadNum);
        }

    private:

        FunctionT &function_;

    } executor(function);

    executor_ = &executor;
    this->cycle(threads);
    executor_ = 0;
}

template<bool CrashOnException>
template<typename FunctionT>
void BasicThreadVector<CrashOnException>::execute(FunctionT function)
{
    this->execute(function, pool_.size());
}

template<bool CrashOnException>
void BasicThreadVector<CrashOnException>::cycle(size_type threads)
{
    boost::unique_lock<boost::mutex> lock(mutex_);
    BCL2FASTQ_ASSERT_MSG(busyThreads_ == 0, "Only one at a time outstanding request is allowed");

    firstThreadException_ = boost::exception_ptr();
    if (threads == 1)
    {
        // Special case for one to simplify debugging. Just do it on the calling thread.
        executor_->execute(0);
    }
    else
    {
        lowestBlockedThreadNumber_ = threads;
        neededThreads_ = threads;
        ++currentRequest_;
        stateChangedCondition_.notify_all();
        this->waitAll(lock);
    }
    if (firstThreadException_)
    {
        BCL2FASTQ_LOG(LogLevel::WARNING) << "Rethrowing a thread exception " << std::endl;
        boost::rethrow_exception(firstThreadException_);
    }
}

template<bool CrashOnException>
void BasicThreadVector<CrashOnException>::waitAll(boost::unique_lock<boost::mutex> &lock)
{
    while (busyThreads_ || neededThreads_)
    {
        stateChangedCondition_.wait(lock);
    }
}

template<bool CrashOnException>
void BasicThreadVector<CrashOnException>::threadFunction(size_type threadNum)
{
    BCL2FASTQ_LOG(LogLevel::DEBUG) << "Thread #" << threadNum << " created" << std::endl;

    boost::unique_lock<boost::mutex> lock(mutex_);
    while (!terminateRequested_)
    {
        BCL2FASTQ_ASSERT_MSG(busyThreads_ > 0, "Thread is not accounted for!!!");
        --busyThreads_;

        const std::size_t processedRequest = currentRequest_;

        //BCL2FASTQ_LOG(LogLevel::DEBUG) << "Thread #" << threadNum << " waiting for new request" << std::endl;
        stateChangedCondition_.notify_all();
        while (!terminateRequested_ && processedRequest == currentRequest_ )
        {
            stateChangedCondition_.wait(lock);
        }
        ++busyThreads_;

        if (!terminateRequested_)
        {
            //BCL2FASTQ_LOG(LogLevel::DEBUG) << "Thread #" << threadNum << " unblocked" << std::endl;
            if (lowestBlockedThreadNumber_ > threadNum)
            {
                BCL2FASTQ_ASSERT_MSG(neededThreads_, "If thread is allowed to run, there must be a need for it!");
                --neededThreads_;

                if (CrashOnException)
                {
                    this->unsafeExecute(threadNum, lock);
                }
                else
                {
                    this->safeExecute(threadNum, lock);
                }
                // we're back under lock
            }
        }
    }

    BCL2FASTQ_LOG(LogLevel::DEBUG) << "Thread #" << threadNum << " terminated" << std::endl;
}

template<bool CrashOnException>
void BasicThreadVector<CrashOnException>::unsafeExecute(size_type threadNum, boost::unique_lock<boost::mutex> &lock)
{
    unlock_guard<boost::unique_lock<boost::mutex> > unlock(lock);
    executor_->execute(threadNum);
}

template<bool CrashOnException>
void BasicThreadVector<CrashOnException>::safeExecute(size_type threadNum, boost::unique_lock<boost::mutex> &lock)
{
    try
    {
        this->unsafeExecute(threadNum, lock);
    }
    catch (...)
    {
        if (!firstThreadException_)
        {
            firstThreadException_ = boost::current_exception();
            BCL2FASTQ_LOG(LogLevel::ERROR_TYPE)
                << "Thread: " << threadNum
                << " caught an exception first: " << boost::current_exception_diagnostic_information()
                << std::endl
            ;
        }
        else
        {
            BCL2FASTQ_LOG(LogLevel::ERROR_TYPE)
                << "Thread: " << threadNum
                << " also caught an exception: " << boost::current_exception_diagnostic_information()
                << std::endl
            ;
        }
    }
}


} // namespace common


} // namespace bcl2fastq


#endif // BCL2FASTQ_COMMON_THREADS_HPP


