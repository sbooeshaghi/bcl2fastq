/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file Task.hh
 *
 * \brief Declaration of a workload item a.k.a. task.
 *
 * \author Marek Balint
 */


#ifndef BCL2FASTQ_CONVERSION_TASK_HH
#define BCL2FASTQ_CONVERSION_TASK_HH


#include <memory>
#include <atomic>
#include <mutex>
#include <condition_variable>
#include <limits>
#include <boost/noncopyable.hpp>

namespace bcl2fastq {
namespace conversion {

class Stage;

class TaskManager
{
public:

    class GroupSerializer
    {
    public:

        /// \brief Default constructor.
        GroupSerializer();

        /// \brief Mutex lock waiting for nextGroupNumber_ to match groupNumber.
        void waitForNextGroupNumber(uint32_t groupNumber);

        /// \brief signal lock is released and group number incremented.
        void signalProceed();

        /// \brief Value of next group number.
        uint32_t GetNextGroupNumber() const;

    private:

        uint32_t nextGroupNumber_;
        std::condition_variable cvGroupNumber_;
        std::mutex groupMutex_;

    }; // end class GroupSerializer


    virtual ~TaskManager() {}
};

/// \brief Task.
class Task : private boost::noncopyable
{
public:

    enum PriorityLevel
    {
        Zero = 0, // Highest priority
        One,
        Two,
        Three,
        Four,
        LOWEST_PRIORITY = Four,
        NUM_PRIORITY_LEVELS
    };

    /// \brief (Default) constructor.
    Task(std::shared_ptr<TaskManager> taskManager);

    /// \brief Virtual destructor.
    virtual ~Task();

    /// \brief Execute the task.
    /// \retval true Continue with execution on current thread.
    /// \retval false Finish current execution iteration on current thread.
    /// \note This function should be overriden in order for task to actually
    /// do something.
    virtual bool execute(int32_t threadNum) = 0;

    virtual PriorityLevel getPriority() const { return PriorityLevel::LOWEST_PRIORITY; }

private:
    std::shared_ptr<TaskManager> taskManager_;
};

typedef std::shared_ptr<Task> TaskPtr;

} // namespace conversion
} // namespace bcl2fastq


#endif // BCL2FASTQ_CONVERSION_TASK_HH


