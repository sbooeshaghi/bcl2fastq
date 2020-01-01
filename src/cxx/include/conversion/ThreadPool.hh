/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file ThreadPool.hh
 *
 * \brief Definition of helpers for thread pool.
 *
 * \author Aaron Day
 */

#ifndef BCL2FASTQ_CONVERSION_THREAD_POOL_HH
#define BCL2FASTQ_CONVERSION_THREAD_POOL_HH

#include "conversion/Stage.hh"
#include "conversion/TaskQueue.hh"
#include "common/SystemCompatibility.hh"
#include <boost/noncopyable.hpp>
#include <vector>
#include <thread>
#include <mutex>
#include <boost/exception_ptr.hpp>

namespace bcl2fastq
{
namespace conversion
{

/*
 * Generic thread pool. Specifiy number of threads on construction. The worker
 * threads will pull tasks off of the TaskQueue and exectue them.
 */

class ThreadPool : boost::noncopyable
{
public:
    ThreadPool(uint32_t   numThreads,
               TaskQueue& dataQueue,
               Terminator& terminator);

    ~ThreadPool() { }

private:
class JoinThreads : boost::noncopyable
{
public:
    JoinThreads(std::vector<std::thread>& threads,
                boost::exception_ptr& firstThreadException);

    ~JoinThreads();

private:
    boost::exception_ptr& firstThreadException_;
    std::vector<std::thread>& threads_;
};

    void workerThread(uint32_t threadNum);

    TaskQueue& dataQueue_;

    Terminator& terminator_;

    /// \brief Thread synchronisation mutex.
    std::mutex mut_;

    /// \brief Exception thrown by worker thread.
    boost::exception_ptr firstThreadException_;

    std::vector<std::thread> threads_;
    JoinThreads              joiner_;
};

} // namespace conversion
} // namespace bcl2fastq

#endif // BCL2FASTQ_CONVERSION_THREAD_POOL_HH
