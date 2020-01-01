/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file Stage.hh
 *
 * \brief Declaration of a processing stage.
 *
 * \author Marek Balint
 * \author Aaron Day
 */


#ifndef BCL2FASTQ_CONVERSION_STAGE_HH
#define BCL2FASTQ_CONVERSION_STAGE_HH


#include <boost/noncopyable.hpp>
#include <boost/thread/mutex.hpp>
#include <future>
#include <mutex>

#include "conversion/Task.hh"
#include "conversion/TaskQueue.hh"


namespace bcl2fastq {
namespace conversion {

class TaskManager;

/// \brief Processing stage.
class Stage : private boost::noncopyable
{
public:

    /// \brief Constructor.
    /// \param stageName Name of the stage.
    explicit Stage(const std::string& stageName,
                   TaskQueue& taskQueue);

    /// \brief Virtual destructor.
    virtual ~Stage() = 0;

    /// \brief Run the processing stage.
    void run();

protected:

    virtual bool startNewTasks() = 0;

    /// \brief Get task queue.
    /// \return Task queue.
    TaskQueue & getTaskQueue();

public:

    /// \brief Announce that an error has occurred and thread is going down.
    /// \note This function should be overridden
    /// by processing stages to let other threads that this thread is going
    /// down and no more data will be either consumed nor produced by it.
    virtual void terminate() = 0;

private:
    /// \brief Task queue.
    TaskQueue& taskQueue_;

protected:
    /// \brief Name of stage
    std::string stageName_;
};

class Terminator : private boost::noncopyable
{
public:
    Terminator() : stages_() { }

    void add(Stage* stage) { stages_.push_back(stage); }
    void terminate() { for (auto& stage : stages_) { stage->terminate(); } }

private:
    std::vector<Stage*> stages_;
};

} // namespace conversion
} // namespace bcl2fastq


#endif // BCL2FASTQ_CONVERSION_STAGE_HH

