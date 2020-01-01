/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file FastqWriter.cpp
 *
 * \brief Implementation of FASTQ writer.
 *
 * \author Marek Balint
 */


#include <utility>

#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>

#include "common/Debug.hh"
#include "conversion/FastqWriter.hh"


namespace bcl2fastq {
namespace conversion {

std::atomic<uint32_t> FastqWriteTaskManager::numTaskManagers_(0);
std::condition_variable FastqWriteTaskManager::cvAllTaskManagersDone_;
std::mutex FastqWriteTaskManager::mut_;

FastqWriteTaskManager::FastqWriteTaskManager(std::shared_ptr<FastqBuffer>& inputBuffer,
                                             ThreadSafeQueue<std::shared_ptr<FastqBuffer>>& inputBuffersToRecycle)
: TaskManager()
, inputBuffer_(inputBuffer)
, inputBuffersToRecycle_(inputBuffersToRecycle)
{
    ++numTaskManagers_;
}

FastqWriteTaskManager::~FastqWriteTaskManager()
{
    static int sequentialId = 0; ++sequentialId;

    // We're done with the inputBuffer; Recycle it.
    inputBuffersToRecycle_.addData(inputBuffer_);
    --numTaskManagers_;

    if (numTaskManagers_ == 0)
    {
        // This could be called more than once on different threads. I think
        // this is ok.
        cvAllTaskManagersDone_.notify_all();
    }
}

void FastqWriteTaskManager::waitForAllTasksToFinish()
{
    std::unique_lock<std::mutex> lock(mut_);
    cvAllTaskManagersDone_.wait(lock, [] { return FastqWriteTaskManager::numTaskManagers_ == 0; });
}

FastqWriteTask::FastqFileStateInfo::FastqFileStateInfo()
: fastqFileStatePtr_()
, nextGroupNumber_(0)
, cvGroupNumber_()
, mtx_()
{
}

FastqWriteTask::FastqFileStateInfo::FastqFileStateInfo(FastqWriteTask::FastqFileStateInfo&& tmp)
: fastqFileStatePtr_(tmp.fastqFileStatePtr_)
, nextGroupNumber_(tmp.nextGroupNumber_)
, cvGroupNumber_()
, mtx_()
{
}

FastqWriteTask::FastqFileStateInfo::FastqFileStateInfo(data::FastqFile::FastqFileStatePtr fastqFileStatePtr)
: fastqFileStatePtr_(fastqFileStatePtr)
, nextGroupNumber_(0)
, cvGroupNumber_()
, mtx_()
{
}

void FastqWriteTask::FastqFileStateInfo::waitForGroupNumber(uint32_t groupNumber)
{
    
    std::unique_lock<std::mutex> lock(mtx_);
    cvGroupNumber_.wait(lock, [this, groupNumber] { return nextGroupNumber_ == groupNumber; });
}

void FastqWriteTask::FastqFileStateInfo::signalWriteComplete()
{
    {
        std::unique_lock<std::mutex> lock(mtx_);
        ++nextGroupNumber_;
    }

    cvGroupNumber_.notify_all();
}

FastqWriteTask::FastqWriteTask(std::shared_ptr<FastqWriteTaskManager>& taskManager,
                               FastqFileStateInfo &fastqFileStateInfo,
                               const FastqBuffer::FastqsContainer::value_type::value_type &inputBuffer,
                               uint32_t groupNumber)
: Task(taskManager)
, fastqFileStateInfo_(fastqFileStateInfo)
, inputBuffer_(inputBuffer)
, groupNumber_(groupNumber)
{
}

bool FastqWriteTask::execute(int32_t)
{
    fastqFileStateInfo_.waitForGroupNumber(groupNumber_);

    for (const auto& read : inputBuffer_)
    {
        data::FastqFile::write(*fastqFileStateInfo_.fastqFileStatePtr_, &*read.begin(), read.size());
    }

    fastqFileStateInfo_.signalWriteComplete();

    return true;
}

FastqWriter::FastqWriter(
    TaskQueue& taskQueue,
    ThreadSafeMap<uint32_t, std::shared_ptr<FastqBuffer>>& inputBuffersToUse,
    ThreadSafeQueue<std::shared_ptr<FastqBuffer>>& inputBuffersToRecycle,
    const layout::Layout &layout,
    const layout::LaneInfo &laneInfo,
    bool createFastqsForIndexReads,
    const boost::filesystem::path& outputDir,
    bool noLaneSplitting
)
: Stage("Fastq writing", taskQueue)
, inputBuffersToUse_(inputBuffersToUse)
, inputBuffersToRecycle_(inputBuffersToRecycle)
, layout_(layout)
, emptyFileDeleter_()
, fastqFileStates_()
, groupNumber_(0)
{
    layout::LaneInfo::SampleInfosContainer::size_type sampleInfosCounts = laneInfo.getSampleInfos().size();

    BOOST_FOREACH (const layout::SampleInfo &sampleInfo, std::make_pair(laneInfo.sampleInfosBegin(), laneInfo.sampleInfosEnd()))
    {
        fastqFileStates_.push_back(FastqFileStatesContainer::value_type());
        for (const auto& readInfo : laneInfo.readInfos())
        {
            if (!readInfo.shouldCreateFastqForRead(createFastqsForIndexReads))
            {
                continue;
            }
            std::ostringstream fileName;
            fileName << (noLaneSplitting ?
                boost::format("%s_S%d_%c%d_001.fastq.gz")
                    % sampleInfo.getName()
                    % sampleInfo.getNumber()
                    % (readInfo.isIndexRead() ? 'I' : 'R')
                    % readInfo.getNumber() :
                boost::format("%s_S%d_%s_%c%d_001.fastq.gz")
                    % sampleInfo.getName()
                    % sampleInfo.getNumber()
                    % laneInfo.getDirName().string()
                    % (readInfo.isIndexRead() ? 'I' : 'R')
                    % readInfo.getNumber())
            ;

            boost::filesystem::path fastqPath(outputDir);
            if (sampleInfo.hasProject())
            {
                fastqPath /= sampleInfo.getProject();
                if( sampleInfo.getId() != sampleInfo.getName() )
                {
                    fastqPath /= sampleInfo.getId();
                }
            }

            boost::filesystem::path filePath(fastqPath / boost::filesystem::path(fileName.str()));
            emptyFileDeleter_.addFile(filePath);

            boost::filesystem::path dirPath = filePath.parent_path();
            if (!dirPath.empty())
            {
                boost::filesystem::create_directories(dirPath);
            }

            fastqFileStates_.back().push_back(std::make_shared<FastqWriteTask::FastqFileStateInfo>( FastqWriteTask::FastqFileStateInfo(
                    data::FastqFile::openFile(filePath,
                                       std::ios_base::out | std::ios_base::binary | (noLaneSplitting && laneInfo.getNumber() > 1 ?
                                           std::ios_base::app :
                                           std::ios_base::trunc )))));
        }
    }
}

void FastqWriter::terminate()
{
    this->getTaskQueue().terminate();
    inputBuffersToUse_.terminate();
}

std::shared_ptr<FastqWriteTaskManager> FastqWriter::getNewTaskManager()
{
    std::shared_ptr<FastqBuffer> inputBuffer;
    if (!inputBuffersToUse_.tryGetData(groupNumber_++, inputBuffer))
    {
        FastqWriteTaskManager::waitForAllTasksToFinish();

        return std::shared_ptr<FastqWriteTaskManager>();
    }

    return std::make_shared<FastqWriteTaskManager>(inputBuffer,
                                                   inputBuffersToRecycle_);
}

bool FastqWriter::startNewTasks()
{
    static int sequentialId = 0; ++sequentialId;
    BCL2FASTQ_LOG(common::LogLevel::DEBUG) << "FastqWriter::startNewTask " << sequentialId  << std::endl;

    TaskQueue &taskQueue = this->getTaskQueue();
    if (taskQueue.isTerminated())
    {
        // Only occurs if there was an error
        return false;
    }

    auto taskManager = getNewTaskManager();
    if (taskManager.get() == NULL)
    {
        return false;
    }

    FastqBuffer &inputBuffer = taskManager->getInputBuffer();

    BCL2FASTQ_LOG(common::LogLevel::DEBUG) << " Start FASTQ-Writer " << inputBuffer.getGroupNumber() << std::endl;

    FastqBuffer::FastqsContainer::value_type::size_type sampleIndex = 0;
    for (const auto& sample : inputBuffer.fastqs_)
    {
        FastqBuffer::FastqsContainer::value_type::value_type::size_type readIndex = 0;
        for (const auto& read : sample)
        {
            taskQueue.addData(std::make_shared<FastqWriteTask>(taskManager,
                                                               *fastqFileStates_.at(sampleIndex).at(readIndex),
                                                               read,
                                                               inputBuffer.getGroupNumber()));

            ++readIndex;
        }

        ++sampleIndex;
    }
    return true;
}

FastqWriter::EmptyFileDeleter::~EmptyFileDeleter()
{
    for (const boost::filesystem::path& filePath : filePaths_)
    {
        if(boost::filesystem::exists(filePath) && boost::filesystem::file_size(filePath) == 0)
        {
            boost::filesystem::remove(filePath);
        }
    }
}


} // namespace conversion
} // namespace bcl2fastq

