/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file Task.cpp
 *
 * \brief Implementation of a thread pool.
 *
 * \author Aaron Day
 */

#include "conversion/ThreadPool.hh"
#include "conversion/TaskQueue.hh"
#include "conversion/Task.hh"
#include "common/Debug.hh"

namespace bcl2fastq
{
namespace conversion
{

ThreadPool::JoinThreads::JoinThreads(std::vector<std::thread>& threads,
                                     boost::exception_ptr& firstThreadException)
 : firstThreadException_(firstThreadException)
 , threads_(threads)
{
}

ThreadPool::JoinThreads::~JoinThreads()
{
    for (std::thread& thread : threads_)
    {
        if (thread.joinable()) thread.join();
    }

    if (firstThreadException_)
    {
        boost::rethrow_exception(firstThreadException_);
    }
}

ThreadPool::ThreadPool(uint32_t   numThreads,
                       TaskQueue& dataQueue,
                       Terminator& terminator)
 : dataQueue_(dataQueue)
 , terminator_(terminator)
 , mut_()
 , firstThreadException_()
 , threads_()
 , joiner_(threads_,
           firstThreadException_)
{
    for (uint32_t i = 0; i < numThreads; ++i)
    {
        threads_.push_back(std::thread(&ThreadPool::workerThread, this, i));
    }
}

void ThreadPool::workerThread(uint32_t threadNum)
{
    try
    {
        TaskPtr task;
        while (dataQueue_.tryGetData(task))
        {
            task->execute(threadNum);
            task.reset();
        }
    }
    catch(...)
    {
        std::lock_guard<std::mutex> lock(mut_);

        if (!firstThreadException_)
        {
            firstThreadException_ = boost::current_exception();

            terminator_.terminate();
            dataQueue_.setFinished();

            BCL2FASTQ_LOG(common::LogLevel::ERROR_TYPE)
                << "Thread: " << threadNum
                << " caught an exception first: " << boost::current_exception_diagnostic_information()
                << std::endl;
        }
        else
        {
            BCL2FASTQ_LOG(common::LogLevel::ERROR_TYPE)
                << "Thread: " << threadNum
                << " also caught an exception: " << boost::current_exception_diagnostic_information()
                << std::endl;
        }
    }
}

} // namespace conversion
} // namespace bcl2fastq

