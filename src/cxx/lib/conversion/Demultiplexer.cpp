/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file Demultiplexer.cpp
 *
 * \brief Implementation of demultiplexer.
 *
 * \author Marek Balint
 * \author Mauricio Varea
 */


#include <algorithm>
#include <numeric>
#include <utility>

#include <boost/foreach.hpp>
#include <boost/assign/list_of.hpp>


#include "common/FileSystem.hh"
#include "layout/Barcode.hh"
#include "conversion/BclBaseConversion.hh"
#include "conversion/Converter.hh"
#include "conversion/Demultiplexer.hh"


namespace bcl2fastq
{
namespace conversion
{

DemultiplexTask::DemuxStatsSubset::DemuxStatsSubset(DemultiplexTask::DemuxStats& summaryDemuxStats)
: summaryDemuxStats_(summaryDemuxStats)
, stats_(summaryDemuxStats.size())
{
    size_t numSamples = summaryDemuxStats.size();
    for (size_t sampleIdx = 0; sampleIdx < numSamples; ++sampleIdx)
    {
        size_t numBarcodes = summaryDemuxStats[sampleIdx].size();
        stats_[sampleIdx].resize(numBarcodes);
        for (size_t barcodeIdx = 0; barcodeIdx < numBarcodes; ++barcodeIdx)
        {
            stats_[sampleIdx][barcodeIdx].resize(summaryDemuxStats[sampleIdx][barcodeIdx].size());
        }
    }
}

DemultiplexTask::DemuxStatsSubset::~DemuxStatsSubset()
{
    // Note: It is important that these are all created/destroyed on the same thread

    size_t numSamples = stats_.size();
    for (size_t sampleIdx = 0; sampleIdx < numSamples; ++sampleIdx)
    {
        size_t numBarcodes = stats_[sampleIdx].size();
        for (size_t barcodeIdx = 0; barcodeIdx < numBarcodes; ++barcodeIdx)
        {
            size_t numTiles = stats_[sampleIdx][barcodeIdx].size();
            for (size_t tileIdx = 0; tileIdx < numTiles; ++tileIdx)
            {
                summaryDemuxStats_[sampleIdx][barcodeIdx][tileIdx] += stats_[sampleIdx][barcodeIdx][tileIdx];
            }
        }
    }
}

std::atomic<uint32_t> DemuxTaskManager::numTaskManagers_(0);
std::condition_variable DemuxTaskManager::cvAllTaskManagersDone_;
std::mutex DemuxTaskManager::mut_;

DemuxTaskManager::DemuxTaskManager(std::shared_ptr<DemuxBuffer>& buffer,
                                   ThreadSafeQueue<std::shared_ptr<DemuxBuffer>>& outputBufferQueue,
                                   TaskManager::GroupSerializer &groupSerializer)
: TaskManager()
, buffer_(buffer)
, outputBufferQueue_(outputBufferQueue)
, groupSerializer_(groupSerializer)
{
    ++numTaskManagers_;
}

DemuxTaskManager::~DemuxTaskManager()
{
    static int sequentialId = 0; ++sequentialId;

    buffer_->calculateIndex();

    BCL2FASTQ_LOG(common::LogLevel::DEBUG) << "dQ Demux Tile: " << buffer_->bclBuffers_.begin()->tileInfo_->getNumber() << std::endl;
    BCL2FASTQ_LOG(common::LogLevel::DEBUG) << " Q FASTQ-Create: " << buffer_->getGroupNumber() << std::endl;

    // Ensure data is passed between queues in order.
    groupSerializer_.waitForNextGroupNumber(buffer_->getGroupNumber());
    outputBufferQueue_.addData(buffer_);
    groupSerializer_.signalProceed();
    --numTaskManagers_;

    if (numTaskManagers_ == 0)
    {
        // This could be called more than once on different threads. I think
        // this is ok.
        cvAllTaskManagersDone_.notify_all();
    }

}

void DemuxTaskManager::waitForAllTasksToFinish()
{
    std::unique_lock<std::mutex> lock(mut_);
    cvAllTaskManagersDone_.wait(lock, [] { return DemuxTaskManager::numTaskManagers_ == 0; });
}


DemultiplexTask::DemultiplexTask(
    std::shared_ptr<DemuxTaskManager>& taskManager,
    const data::BclBuffer &inputBuffer,
    const layout::ReadInfosContainer& readInfos,
    data::BclBuffer::BclCycleContainer::size_type offsetBegin,
    data::BclBuffer::BclCycleContainer::size_type offsetEnd,
    const layout::BarcodeTranslationTable &barcodeTranslationTable,
    layout::LaneInfo::TileInfosContainer::const_iterator tileInfosBegin,
    size_t numIndexReads,
    bool includeNonPfClusters,
    common::NumBasesPerByte numBasesPerByte,
    std::vector<std::shared_ptr<DemultiplexTask::DemuxStatsSubset>>& demuxStats,
    data::BclBuffer::SamplesContainer &outputBuffer
)
: Task(taskManager)
, inputBuffer_(inputBuffer)
, readInfos_(readInfos)
, offsetBegin_(offsetBegin)
, offsetEnd_(offsetEnd)
, barcodeTranslationTable_(barcodeTranslationTable)
, tileInfosBegin_(tileInfosBegin)
, numIndexReads_(numIndexReads)
, includeNonPfClusters_(includeNonPfClusters)
, numBasesPerByte_(numBasesPerByte)
, demuxStats_(demuxStats)
, outputBuffer_(outputBuffer)
{
}

bool DemultiplexTask::execute(int32_t threadNum)
{
    switch(numBasesPerByte_)
    {
    case common::NumBasesPerByte::ONE:
        return execute<common::NumBasesPerByte::ONE>(threadNum);
        break;
    case common::NumBasesPerByte::TWO:
        return execute<common::NumBasesPerByte::TWO>(threadNum);
        break;
    case common::NumBasesPerByte::FOUR:
        return execute<common::NumBasesPerByte::FOUR>(threadNum);
        break;
    default:
        BCL2FASTQ_ASSERT_MSG(false, "Unsuported number of bases per byte. Supported values include: One, Two, or Four");
        break;
    };

    return false;
}

Demultiplexer::Demultiplexer(
    TaskQueue& taskQueue,
    ThreadSafeQueue<std::shared_ptr<DemuxBuffer>>& inputBuffersToUse,
    ThreadSafeQueue<std::shared_ptr<DemuxBuffer>>& outputBuffersToSubmit,
    const layout::Layout &layout,
    const layout::LaneInfo &laneInfo,
    const layout::BarcodeTranslationTable &barcodeTranslationTable,
    bool includeNonPfClusters,
    int32_t numThreads,
    DemultiplexTask::DemuxStats &summaryDemuxStats
)
: Stage("Demultiplexing", taskQueue)
, inputBuffersToUse_(inputBuffersToUse)
, outputBuffersToSubmit_(outputBuffersToSubmit)
, layout_(layout)
, laneInfo_(laneInfo)
, barcodeTranslationTable_(barcodeTranslationTable)
, includeNonPfClusters_(includeNonPfClusters)
, summaryDemuxStats_(summaryDemuxStats)
, doDemultiplex_(false)
, numIndexReads_(0)
, stats_(numThreads)
, groupSerializer_()
{
    layout::LaneInfo::SampleInfosContainer::size_type sampleIndex = 0;
    for (const layout::SampleInfo &sampleInfo : laneInfo_.getSampleInfos())
    {
        summaryDemuxStats_.push_back(DemultiplexTask::DemuxStats::value_type());

        layout::BarcodesContainer::size_type barcodesCount = sampleInfo.getBarcodes().size();
        if( 0==barcodesCount )
        {
            summaryDemuxStats_.back().push_back(DemultiplexTask::DemuxStats::value_type::value_type());
            createTileStats( layout::BarcodeTranslationTable::SampleMetadata(sampleIndex) );
        } else {
            for( layout::BarcodesContainer::size_type barcodeIndex = 0;
                 barcodeIndex < barcodesCount;
                 ++barcodeIndex )
            {
                summaryDemuxStats_.back().push_back(DemultiplexTask::DemuxStats::value_type::value_type());
                createTileStats( layout::BarcodeTranslationTable::SampleMetadata(sampleIndex,barcodeIndex) );
            }
        }
        ++sampleIndex;
    }

    if (laneInfo.getSampleInfos().size() > 1)
    {
        for (const auto& readInfo : laneInfo.readInfos())
        {
            if ((readInfo.isIndexRead()) && !readInfo.cycleInfos().empty())
            {
                ++numIndexReads_;
                doDemultiplex_ = true;
            }
        }
    }
    else
    {
        BCL2FASTQ_ASSERT_MSG( laneInfo.getSampleInfos().size() == 1,
                              "There must be at least default sample" );
        // in case there is sample with barcode defined in sample sheet, there will be also default sample (i.e. at least 2 samples in total) and this branch will never get executed
        BCL2FASTQ_ASSERT_MSG( laneInfo.sampleInfosBegin()->getBarcodes().empty(),
                              "There should be either default sample with no barcode or the only sample defined in sample sheet without barcode");
    }

    for (auto& stats : stats_)
    {
        stats = std::make_shared<DemultiplexTask::DemuxStatsSubset>(summaryDemuxStats_);
    }
}

Demultiplexer::~Demultiplexer()
{
}

void Demultiplexer::createTileStats(const layout::BarcodeTranslationTable::SampleMetadata &sampleMetadata)
{
    layout::LaneInfo::TileInfosContainer::size_type tilesCount = laneInfo_.getTileInfos().size();
    for( layout::LaneInfo::TileInfosContainer::size_type tileIndex = 0;
         tileIndex < tilesCount;
         ++tileIndex )
    {
        summaryDemuxStats_.back().back().push_back(DemultiplexTask::DemuxStats::value_type::value_type::value_type(sampleMetadata,tileIndex));
    }
}

bool Demultiplexer::startNewTasks()
{
    static int sequentialId = 0; ++sequentialId;
    BCL2FASTQ_LOG(common::LogLevel::DEBUG) << "Demultiplexer::startNewTask " << sequentialId  << std::endl;

    TaskQueue &taskQueue = this->getTaskQueue();
    if (taskQueue.isTerminated())
    {
        outputBuffersToSubmit_.terminate();

        // Only occurs if there was an error
        return false;
    }

    auto taskManager = getNewTaskManager();
    if (taskManager.get() == NULL)
    {
        DemuxTaskManager::waitForAllTasksToFinish();
        outputBuffersToSubmit_.setFinished();
        return false;
    }

    auto &outputBufferVec = taskManager->getBuffer().bclBuffers_;

    // If no input, clear the output and return true...
    bool hasCycleData = false;
    for (data::BclBuffer &bclBuffer : outputBufferVec) {
        hasCycleData = bclBuffer.hasCycleData();
        if (hasCycleData) break;
    }
    if (!hasCycleData) {
        return true;
    }

    calcFilterOffsets(outputBufferVec);

    uint32_t taskNum = 0;

    for (auto& outputBuffer : outputBufferVec)
    {
        BCL2FASTQ_ASSERT_MSG(outputBuffer.bcls_.size() > 0, "No BCL cycles to process");

        if (doDemultiplex_)
        {
            const data::BclBuffer::BclCycleContainer::size_type dataSize = outputBuffer.bcls_.front().getNumClusters();
            BCL2FASTQ_LOG(common::LogLevel::DEBUG) << "Queue new demux tasks " << sequentialId  << std::endl;
            for (data::BclBuffer::BclCycleContainer::size_type offset = 0; offset < dataSize; offset += data::ClustersPerTask)
            {
                taskQueue.addData(std::make_shared<DemultiplexTask>(
                    taskManager,
                    outputBuffer,
                    laneInfo_.readInfos(),
                    offset,
                    std::min(offset+data::ClustersPerTask, dataSize),
                    barcodeTranslationTable_,
                    laneInfo_.getTileInfos().begin(),
                    numIndexReads_,
                    includeNonPfClusters_,
                    layout_.getNumBasesPerByte(),
                    stats_,
                    outputBuffer.samples_
                ));

                ++taskNum;
            }
        }
    }

    return true;
}

// TODO: Make this multithreaded
void Demultiplexer::calcFilterOffsets(data::BclBufferVec& buffers)
{
    size_t bufferIndex = 0;
    for (auto& buffer : buffers)
    {
        if (buffer.bcls_.back().includeNonPf_)
        {
            continue;
        }

        // We reserve a bit more than the necessary size in general to avoid
        // reallocating when we find a bigger buffer.
        if (buffer.cbclFilterOffsets_.empty())
        {
            buffer.cbclFilterOffsets_.reserve(buffer.filters_.capacity());
        }

        buffer.cbclFilterOffsets_.resize(buffer.filters_.size(), 0);

        size_t passedFilterCount = 0;
        size_t filterIndex = 0;
        for (auto& filterRecord : buffer.filters_)
        {
            if (filterRecord.data_)
            {
                // Passed filter
                buffer.cbclFilterOffsets_[filterIndex] = passedFilterCount;
                ++passedFilterCount;
            }

            ++filterIndex;
        }

        ++bufferIndex;
    }
}

void Demultiplexer::terminate()
{
    this->getTaskQueue().terminate();
    inputBuffersToUse_.terminate();
}

std::shared_ptr<DemuxTaskManager> Demultiplexer::getNewTaskManager()
{
    std::shared_ptr<DemuxBuffer> inputBuffer;
    if (!inputBuffersToUse_.tryGetData(inputBuffer))
    {
        // Only get here if we're done. Otherwise, tryGetData would wait.
        return std::shared_ptr<DemuxTaskManager>();
    }

    BCL2FASTQ_LOG(common::LogLevel::TRACE) << "Group: " << inputBuffer->getGroupNumber() << " New DemuxBuffer" << std::endl;
    return std::make_shared<DemuxTaskManager>(inputBuffer,
                                              outputBuffersToSubmit_, groupSerializer_);
}

} // namespace conversion
} // namespace bcl2fastq


