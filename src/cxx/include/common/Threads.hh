/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file Threads.hh
 *
 * \brief Definition of helpers for thread management.
 *
 * \author Roman Petorvski
 * \author Marek Balint
 */


#ifndef BCL2FASTQ_COMMON_THREADS_HH
#define BCL2FASTQ_COMMON_THREADS_HH


#include <cstddef>

#include <thread>
#include <boost/noncopyable.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/condition_variable.hpp>
#include <boost/exception_ptr.hpp>


namespace bcl2fastq {


namespace common {


/// \brief Inversion of the boost::lock_guard.
/// \todo Refactoring: Use Boost::ASIO: http://mostlycoding.blogspot.de/2009/05/asio-library-has-been-immensely-helpful.html
template<typename Mutex>
class unlock_guard : private boost::noncopyable
{
public:

    /// \brief Constructor.
    /// \param m Mutex.
    explicit unlock_guard(Mutex& m);

    /// \brief Destructor.
    ~unlock_guard();

private:

    /// \brief Mutex.
    Mutex& m_;
};


/// \brief Thread pool.
template<bool CrashOnException>
class BasicThreadVector : boost::noncopyable
{
public:

    /// \brief Size type definition.
    typedef std::size_t size_type;

public:

    /// \brief Constructor.
    /// \param size Number of threads (a.k.a pool size).
    explicit BasicThreadVector(size_type size);

    /// \brief Destructor.
    ~BasicThreadVector();

public:

    /// \brief Get pool size.
    /// \return Pool size.
    size_type size() const;

public:

    /// \brief Execute given function on given number of threads.
    /// \param function The function to be executed.
    /// \param threads Number of threads to be used.
    /// \pre <tt>threads <= this->size()</tt>
    template<typename FunctionT>
    void execute(FunctionT function, size_type threads);

    /// \brief Execute given function on all threads.
    /// \param function The function to be executed.
    template<typename FunctionT>
    void execute(FunctionT function);

private:

    /// \brief Signal threads to execute workload.
    /// \param threads Number of threads to use for workload execution.
    void cycle(size_type threads);

    /// \brief Wait for all threads to finish their current task.
    /// \param lock Synchronisation lock.
    void waitAll(boost::unique_lock<boost::mutex> &lock);

    /// \brief Thread function for individual threads.
    /// \param threadNum Thread number (0-based).
    void threadFunction(size_type threadNum);

    /// \brief Execute workload by the thread.
    /// \param threadNum Thread number (0-based).
    /// \param lock Synchronisation lock.
    /// \note No exception handling.
    void unsafeExecute(size_type threadNum, boost::unique_lock<boost::mutex> &lock);

    /// \brief Execute workload by the thread.
    /// \param threadNum Thread number (0-based).
    /// \param lock Synchronisation lock.
    /// \note With exception handling (exception thrown by the worker thread
    /// will be rethrown on calling thread).
    void safeExecute(size_type threadNum, boost::unique_lock<boost::mutex> &lock);

private:

    /// \brief Thread workload executor.
    class Executor
    {
    public:

        /// \brief Virtual destructor.
        virtual ~Executor() {}

    public:

        /// \brief Execute thread workload.
        /// \param threadNum Thread number (0-based).
        virtual void execute(size_type threadNum) = 0;
    };

private:

    /// \brief Container of threads.
    boost::ptr_vector<std::thread> pool_;

    /// \brief Thread synchronisation mutex.
    boost::mutex mutex_;

    /// \brief Thread synchronisation condition variable.
    boost::condition_variable stateChangedCondition_;

    /// \brief Current thread workload executor.
    Executor *executor_;

    /// \brief Number of threads currently processing the request.
    size_type busyThreads_;

    /// \brief Number of threads still required to process the request
    size_type neededThreads_;

    /// \brief Lowest idle thread
    /// \note When executing with less threads than available,
    /// this prevents the higher number threads to carry out the request.
    size_type lowestBlockedThreadNumber_;

    /// \brief True when the whole thing goes down.
    bool terminateRequested_;

    /// \brief Number of the current request.
    /// \note Constantly-incrementing number to make sure
    /// each thread processes one master call only once.
    std::size_t currentRequest_;

    /// \brief Exception thrown by worker thread.
    boost::exception_ptr firstThreadException_;
};

/// \brief Thread pool with exception handling type definition.
typedef BasicThreadVector<false> SafeThreadVector;

/// \brief Thread pool without exception handling type definition.
typedef BasicThreadVector<true> UnsafeThreadVector;

/// \brief Thread pool type definition.
typedef SafeThreadVector ThreadVector;


} // namespace common


} // namespace bcl2fastq


#include "common/Threads.hpp"


#endif // BCL2FASTQ_COMMON_THREADS_HH


