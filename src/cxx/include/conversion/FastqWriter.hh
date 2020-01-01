/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file FastqWriter.hh
 *
 * \brief Declaration of FASTQ writer.
 *
 * \author Marek Balint
 */


#ifndef BCL2FASTQ_CONVERSION_FASTQWRITER_HH
#define BCL2FASTQ_CONVERSION_FASTQWRITER_HH


#include "layout/Layout.hh"
#include "data/FastqFile.hh"
#include "conversion/FastqBuffer.hh"
#include "conversion/Stage.hh"
#include "conversion/Task.hh"

#include <condition_variable>
#include <mutex>

namespace bcl2fastq {


namespace conversion {

class FastqWriteTaskManager : public TaskManager
{
public:
    FastqWriteTaskManager(std::shared_ptr<FastqBuffer>& inputBuffer,
                          ThreadSafeQueue<std::shared_ptr<FastqBuffer>>& inputBuffersToRecycle);

    virtual ~FastqWriteTaskManager();

    static void waitForAllTasksToFinish();

    FastqBuffer& getInputBuffer() { return *inputBuffer_; }

private:
    std::shared_ptr<FastqBuffer> inputBuffer_;
    ThreadSafeQueue<std::shared_ptr<FastqBuffer>>& inputBuffersToRecycle_;

    static std::atomic<uint32_t> numTaskManagers_;
    static std::condition_variable cvAllTaskManagersDone_;
    static std::mutex mut_;
};

/// \brief Task: Write FASTQs from buffer to file.
class FastqWriteTask : public Task
{
public:

    /// \brief FASTQ file external state with metadata type definition.
    struct FastqFileStateInfo
    {
    public:

        /// \brief Default constructor.
        FastqFileStateInfo();

        /// \brief Constructor.
        explicit FastqFileStateInfo(data::FastqFile::FastqFileStatePtr fastqFileStatePtr);

        FastqFileStateInfo(FastqFileStateInfo&& tmp);

        void waitForGroupNumber(uint32_t groupNumber);

        void signalWriteComplete();

    public:

        /// \brief FASTQ file external state.
        data::FastqFile::FastqFileStatePtr fastqFileStatePtr_;

        uint32_t nextGroupNumber_;
        std::condition_variable cvGroupNumber_;
        std::mutex mtx_;
    };

public:

    /// \brief Constructor.
    /// \param fastqFiles FASTQ files (one per thread).
    /// \param fastqFileStateInfo External state and metadata of FASTQ file to be written to.
    FastqWriteTask(std::shared_ptr<FastqWriteTaskManager>& taskManager,
                   FastqFileStateInfo &fastqFileStateInfo,
                   const FastqBuffer::FastqsContainer::value_type::value_type &inputBuffer,
                   uint32_t groupNumber);

public:

    virtual bool execute(int32_t threadNum);

    virtual PriorityLevel getPriority() const { return PriorityLevel::Zero; }

private:

    /// \brief File state meta data of the FASTQ file to be written to.
    FastqFileStateInfo& fastqFileStateInfo_;

    /// \brief Buffer to read data from.
    const FastqBuffer::FastqsContainer::value_type::value_type& inputBuffer_;

    uint32_t groupNumber_;
};


/// \brief FASTQ writer.
class FastqWriter : public Stage // public SinkStage<FastqBuffer>
{
public:

    /// \brief Constructor.
    /// \param threadsCount Number of threads.
    /// \param layout Flowcell layout.
    /// \param laneInfo Lane meata data.
    /// \param createFastqsForIndexReads Create FASTQ files also for index reads flag.
    /// \param outputDir Output directory.
    FastqWriter(
        TaskQueue& taskQueue,
        ThreadSafeMap<uint32_t, std::shared_ptr<FastqBuffer>>& inputBuffersToUse,
        ThreadSafeQueue<std::shared_ptr<FastqBuffer>>& inputBuffersToRecycle,
        const layout::Layout &layout,
        const layout::LaneInfo &laneInfo,
        bool createFastqsForIndexReads,
        const boost::filesystem::path& outputDir,
        bool noLaneSplitting
    );

public:

    virtual std::shared_ptr<FastqWriteTaskManager> getNewTaskManager();

    virtual bool startNewTasks();

private:

    virtual void terminate();

    /// \brief FASTQ file external states container type definition.
    typedef std::vector<std::vector<std::shared_ptr<FastqWriteTask::FastqFileStateInfo>>> FastqFileStatesContainer;

private:
// Class to delete an empty file on destruction
class EmptyFileDeleter
{
public:
    EmptyFileDeleter() : filePaths_() { }
    ~EmptyFileDeleter();

    void addFile(boost::filesystem::path& filePath) { filePaths_.push_back(filePath); }

private:
    std::vector<boost::filesystem::path> filePaths_;
};

    ThreadSafeMap<uint32_t, std::shared_ptr<FastqBuffer>>& inputBuffersToUse_;
    ThreadSafeQueue<std::shared_ptr<FastqBuffer>>& inputBuffersToRecycle_;

    /// \brief Layout.
    const layout::Layout &layout_;

    /// \brief Delete empty files on destruction
    EmptyFileDeleter emptyFileDeleter_;

    /// \brief FASTQ file states.
    FastqFileStatesContainer fastqFileStates_;

    uint32_t groupNumber_;
};


} // namespace conversion


} // namespace bcl2fastq


#endif // BCL2FASTQ_CONVERSION_FASTQWRITER_HH


