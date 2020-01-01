/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file BclReader.cpp
 *
 * \brief Implementation of BCL reader.
 *
 * \author Marek Balint
 * \author Mauricio Varea
 */

#include "conversion/BclReader.hh"
#include "data/CbclFile.hh"
#include "data/TileBclFileReader.hh"
#include "data/AggregatedBclFileReader.hh"
#include "data/BclFile.hh"
#include "data/Instrument.hh"

#include "common/Debug.hh"
#include "common/Exceptions.hh"
#include "common/Types.hh"
#include "io/GzipCompressor.hh"

#include <errno.h>
#include <algorithm>
#include <utility>

#include <boost/format.hpp>
#include <boost/filesystem.hpp>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

namespace bcl2fastq {
namespace conversion {

namespace 
{
    std::string getTileSpec(
        const layout::LaneInfo &laneInfo,
        const layout::TileInfo &tileInfo)
    {
        std::stringstream ss;
        ss << "Lane: " << laneInfo.getNumber() << " Tile: " << tileInfo.getNumber();

        return ss.str();
    }
}


BclReadTask::BclReadTask(
    std::shared_ptr<BclReaderTaskManager>& taskManager,
    data::BclFileReader& bclFileReader,
    data::RawBclBufferGroup& outputBuffer
)
: Task(taskManager)
, bclFileReader_(bclFileReader)
, outputBuffer_(outputBuffer)
{
}

bool BclReadTask::execute(int32_t)
{
    return bclFileReader_.read(outputBuffer_);
}

PositionsReadTask::PositionsReadTask(
    std::shared_ptr<BclReaderTaskManager>& taskManager,
    std::shared_ptr<io::SyncFile> &positionsFile,
    const boost::filesystem::path &intensitiesDir,
    bool aggregateTilesFlag,
    bool isPatternedFlowcell,
    const layout::LaneInfo &laneInfo,
    layout::LaneInfo::TileInfosContainer::const_iterator &tileInfoIter,
    bool ignoreMissingPositions,
    data::RawBclBufferGroup& outputBuffer,
    std::shared_ptr<PatternedPositionsContainer>& patternedFlowcellPositions
)
: Task(taskManager)
, positionsFile_(positionsFile)
, intensitiesDir_(intensitiesDir)
, aggregateTilesFlag_(aggregateTilesFlag)
, isPatternedFlowcell_(isPatternedFlowcell)
, laneInfo_(laneInfo)
, tileInfoIter_(tileInfoIter)
, ignoreMissingPositions_(ignoreMissingPositions)
, outputBuffer_(outputBuffer)
, patternedFlowcellPositions_(patternedFlowcellPositions)
{
}

bool PositionsReadTask::execute(int32_t)
{
    static int sequentialId = 0; ++sequentialId;
    BCL2FASTQ_LOG(common::LogLevel::DEBUG) << "PositionsReadTask::execute " << sequentialId  << std::endl;
    if (isPatternedFlowcell_ && patternedFlowcellPositions_)
    {
        if (readEntirePosFile(*patternedFlowcellPositions_))
        {
            patternedFlowcellPositions_->setReady();
            return true;
        }
        else
        {
            patternedFlowcellPositions_->setReady();
            return false;
        }
    }
    else
    {
        for (auto& outputBuffer : outputBuffer_)
        {
            if (!outputBuffer.positions_)
            {
                outputBuffer.positions_ = std::make_shared<common::RawDataBuffer>();
            }

            if (!(aggregateTilesFlag_ ? readPartOfAggregatedPosFile(*outputBuffer.positions_) : readEntirePosFile(*outputBuffer.positions_)))
            {
                return false;
            }

            ++tileInfoIter_;
        }
    }

    BCL2FASTQ_LOG(common::LogLevel::DEBUG) << "PositionsReadTask::execute done" << sequentialId  << std::endl;
    return true;
}

bool PositionsReadTask::readEntirePosFile(common::RawDataBuffer& outputBuffer)
{
    static int sequentialId = 0; ++sequentialId;

    size_t headerSize;
    boost::filesystem::path posFilePath;
    if (!data::PositionsFileFactory::doesFileExist(intensitiesDir_,
                                                   aggregateTilesFlag_,
                                                   isPatternedFlowcell_,
                                                   laneInfo_.getNumber(),
                                                   tileInfoIter_->getNumber(),
                                                   headerSize,
                                                   posFilePath))
    {
        BCL2FASTQ_LOG(common::LogLevel::DEBUG) << "FAIL PositionsReadTask readEntireFile " << sequentialId  << " " << posFilePath.string() << std::endl;
        if (ignoreMissingPositions_)
        {
            return true;
        }
        else
        {
            BOOST_THROW_EXCEPTION(common::IoError(ENOENT, "Unable to find positions file: '" + posFilePath.string() + "'"));
        }
    }

    BCL2FASTQ_LOG(common::LogLevel::DEBUG) << "PositionsReadTask readEntireFile " << sequentialId  << " " << posFilePath.string() << std::endl;

    io::UnprocessedFile positionsFile(posFilePath,
                                      ignoreMissingPositions_);

    outputBuffer.path_ = positionsFile.getPath();
    bool success = positionsFile.readEntireFile(outputBuffer);
    return success;
}

bool PositionsReadTask::readPartOfAggregatedPosFile(common::RawDataBuffer& outputBuffer)
{
    static int sequentialId = 0; ++sequentialId;

    io::SyncFile::SyncFileReader fileReader(*positionsFile_,
                                            std::distance(laneInfo_.getTileInfos().begin(), tileInfoIter_));

    size_t headerSize = 0;
    if (!positionsFile_->isOpen())
    {
        boost::filesystem::path filePath;
        if (!(data::PositionsFileFactory::doesFileExist(intensitiesDir_,
                                                        aggregateTilesFlag_,
                                                        isPatternedFlowcell_,
                                                        laneInfo_.getNumber(),
                                                        tileInfoIter_->getNumber(),
                                                        headerSize,
                                                        filePath)))
        {
            if (ignoreMissingPositions_)
            {
    BCL2FASTQ_LOG(common::LogLevel::DEBUG) << "SKIP PositionsReadTask readPart...File " << sequentialId  << " " << filePath.string() << std::endl;
                return true;
            }
            else
            {
                BOOST_THROW_EXCEPTION(common::IoError(ENOENT, "Unable to find position file: '" + filePath.string() + "'"));
            }
        }

        fileReader.openFile(filePath,
                            ignoreMissingPositions_);
    }

    BCL2FASTQ_LOG(common::LogLevel::DEBUG) << "PositionsReadTask readPart...File " << sequentialId << std::endl;
    std::size_t skippedClustersCount = tileInfoIter_->getSkippedClustersCount();
    bool skipHeader = (tileInfoIter_ == laneInfo_.getTileInfos().begin()) && (tileInfoIter_->getIndex() != 0);
    if (skippedClustersCount || skipHeader)
    {
        fileReader.seek(skippedClustersCount*sizeof(data::PositionsFile::Record) + (skipHeader ? headerSize : 0));
    }

    outputBuffer.path_ = positionsFile_->getPath();
    fileReader.read(outputBuffer,
                    (skipHeader ? 0 : headerSize) + sizeof(data::PositionsFile::Record)*tileInfoIter_->getClustersCount());

    return true;
}

FilterReadTask::FilterReadTask(
    std::shared_ptr<BclReaderTaskManager>& taskManager,
    std::shared_ptr<io::SyncFile> &filterFile,
    const boost::filesystem::path &inputDir,
    common::TileAggregationMode tileAggregationMode,
    const layout::LaneInfo &laneInfo,
    layout::LaneInfo::TileInfosContainer::const_iterator &tileInfoIter,
    bool ignoreMissingFilters,
    data::RawBclBufferGroup& outputBuffer
)
: Task(taskManager)
, filterFile_(filterFile)
, inputDir_(inputDir)
, tileAggregationMode_(tileAggregationMode)
, aggregateTilesFlag_(data::Instrument::isFilterFileAggregated(tileAggregationMode))
, laneInfo_(laneInfo)
, tileInfoIter_(tileInfoIter)
, ignoreMissingFilters_(ignoreMissingFilters)
, outputBuffer_(outputBuffer)
{
}

bool FilterReadTask::execute(int32_t)
{
    static int sequentialId = 0; ++sequentialId;
    BCL2FASTQ_LOG(common::LogLevel::DEBUG) << "FilterReadTask::execute " << sequentialId  << std::endl;
    for (data::RawBclBuffer& rawBclBuffer : outputBuffer_)
    {
        if (!(aggregateTilesFlag_ ? readPartOfAggregatedFilterFile(rawBclBuffer.filters_) : readEntireFilterFile(rawBclBuffer)))
        {
            return false;
        }

        ++tileInfoIter_;
    }

    BCL2FASTQ_LOG(common::LogLevel::DEBUG) << "FilterReadTask::execute done" << sequentialId  << std::endl;
    return true;
}

bool FilterReadTask::readEntireFilterFile(data::RawBclBuffer& rawBclBuffer)
{
    static int sequentialId = 0; ++sequentialId;

    common::RawDataBuffer& rawFiltersBuffer = rawBclBuffer.filters_;

    // Ensure there isn't any data from previous tiles
    rawFiltersBuffer.clear();

    size_t headerSize;
    boost::filesystem::path filePath;
    if (!(data::FilterFile::doesFileExist(inputDir_,
                                          aggregateTilesFlag_,
                                          laneInfo_.getNumber(),
                                          tileInfoIter_->getNumber(),
                                          headerSize,
                                          filePath)))
    {
        if (ignoreMissingFilters_)
        {
            if (tileAggregationMode_ == common::TileAggregationMode::CBCL)
            {
                // If cbcl filter file is missing then we can't map the tile indices so ignore this tile. 
                BCL2FASTQ_LOG(common::LogLevel::WARNING) << "Unable to find CBCL filter file: '" 
                        << filePath.string()  << "'. This tile will be ignored." << std::endl;
            }
            return true;
        }
        else
        {
            BOOST_THROW_EXCEPTION(common::IoError(ENOENT, "Unable to find filter file: '" + filePath.string() + "'"));
        }
    }

    BCL2FASTQ_LOG(common::LogLevel::DEBUG) << "FilterReadTask::execute " << sequentialId << " " << filePath.string() << std::endl;

    bool status =  false;
    if (!aggregateTilesFlag_) 
    {
        // A.  NO WAIT APPROACH 
        io::UnprocessedFile filterFile(filePath,
                                       ignoreMissingFilters_);

        status = filterFile.readEntireFile(rawFiltersBuffer);
    }
    else 
    {
        io::SyncFile::SyncFileReader fileReader(*filterFile_,
                                                std::distance(laneInfo_.getTileInfos().begin(), tileInfoIter_));

        // B.  Sync read.
        fileReader.openFile(filePath,
                            ignoreMissingFilters_);
        status = fileReader.readEntireFile(rawFiltersBuffer);
    }

    rawFiltersBuffer.path_ = filePath;

    BCL2FASTQ_LOG(common::LogLevel::DEBUG) << "...done FilterReadTask::execute " << sequentialId  << std::endl;

    return status;
}

bool FilterReadTask::readPartOfAggregatedFilterFile(common::RawDataBuffer& outputBuffer)
{
    io::SyncFile::SyncFileReader fileReader(*filterFile_,
                                            std::distance(laneInfo_.getTileInfos().begin(), tileInfoIter_));

    size_t headerSize = 0;
    if (!filterFile_->isOpen())
    {
        boost::filesystem::path filePath;
        if (!(data::FilterFile::doesFileExist(inputDir_,
                                              aggregateTilesFlag_,
                                              laneInfo_.getNumber(),
                                              tileInfoIter_->getNumber(),
                                              headerSize,
                                              filePath)))
        {
            if (ignoreMissingFilters_)
            {
                return true;
            }
            else
            {
                BOOST_THROW_EXCEPTION(common::IoError(ENOENT, "Unable to find filter file: '" + filePath.string() + "'"));
            }
        }

        fileReader.openFile(filePath,
                            ignoreMissingFilters_);
    }


    bool skipHeader = (tileInfoIter_ == laneInfo_.getTileInfos().begin()) && (tileInfoIter_->getIndex() != 0);
    std::size_t skippedClustersCount = tileInfoIter_->getSkippedClustersCount();
    if (skippedClustersCount || skipHeader)
    {
        fileReader.seek(skippedClustersCount + (skipHeader ? headerSize : 0));
    }

    outputBuffer.path_ = filterFile_->getPath();
    fileReader.read(outputBuffer,
                    (skipHeader ? 0 : headerSize) + tileInfoIter_->getClustersCount());

    return true;
}

ControlReadTask::ControlReadTask(
    std::shared_ptr<BclReaderTaskManager>& taskManager,
    const boost::filesystem::path &inputDir,
    const layout::LaneInfo &laneInfo,
    layout::LaneInfo::TileInfosContainer::const_iterator &tileInfoIter,
    bool ignoreMissingControls,
    data::RawBclBufferGroup& outputBuffer
)
: Task(taskManager)
, inputDir_(inputDir)
, laneInfo_(laneInfo)
, tileInfoIter_(tileInfoIter)
, ignoreMissingControls_(ignoreMissingControls)
, outputBuffer_(outputBuffer)
{
}

bool ControlReadTask::execute(int32_t)
{
    for (auto& outputBuffer : outputBuffer_)
    {
        outputBuffer.controls_.clear();

        if (!readEntireControlFile(outputBuffer.controls_))
        {
            return false;
        }

        ++tileInfoIter_;
    }

    return true;
}

bool ControlReadTask::readEntireControlFile(common::RawDataBuffer& outputBuffer)
{
    boost::filesystem::path filePath;
    if (!(data::ControlFile::doesFileExist(inputDir_,
                                           laneInfo_.getNumber(),
                                           tileInfoIter_->getNumber(),
                                           filePath)))
    {
        return true;
    }

    io::UnprocessedFile controlFile(filePath,
                                    ignoreMissingControls_);

    return controlFile.readEntireFile(outputBuffer);
}
 
std::atomic<uint32_t> BclReaderTaskManager::numTaskManagers_(0);
std::condition_variable BclReaderTaskManager::cvAllTaskManagersDone_;
std::mutex BclReaderTaskManager::mut_;

BclReaderTaskManager::BclReaderTaskManager(bool ignoreMissingBcls,
                                           bool ignoreMissingFilters,
                                           bool ignoreMissingPositions,
                                           bool ignoreMissingControls,
                                           std::shared_ptr<data::RawBclBufferGroup>& outputBuffer,
                                           ThreadSafeQueue<std::shared_ptr<data::RawBclBufferGroup>>& outputBufferQueue,
                                           std::shared_ptr<PatternedPositionsContainer>& patternedFlowcellPositions,
                                           TaskManager::GroupSerializer &groupSerializer)
: TaskManager()
, ignoreMissingBcls_(ignoreMissingBcls)
, ignoreMissingFilters_(ignoreMissingFilters)
, ignoreMissingPositions_(ignoreMissingPositions)
, ignoreMissingControls_(ignoreMissingControls)
, outputBuffer_(outputBuffer)
, outputBufferQueue_(outputBufferQueue)
, patternedFlowcellPositions_(patternedFlowcellPositions)
, groupSerializer_(groupSerializer)
{
    ++numTaskManagers_;
}

BclReaderTaskManager::~BclReaderTaskManager()
{
    // All tasks are complete.
    postExecute();

    BCL2FASTQ_LOG(common::LogLevel::DEBUG) << " Q Raw Tile: " << outputBuffer_->begin()->tileInfo_->getNumber() << std::endl;

    // Inconsistent latency in reading files can result in a dead-lock writing fastq files.
    // Wait to ensure tiles are loaded in order and minimize latency issues.
    groupSerializer_.waitForNextGroupNumber(outputBuffer_->getGroupNumber());
    outputBufferQueue_.addData(outputBuffer_);
    groupSerializer_.signalProceed();
    --numTaskManagers_;

    if (numTaskManagers_ == 0)
    {
        // This could be called more than once on different threads. I think
        // this is ok.
        cvAllTaskManagersDone_.notify_all();
    }
}

void BclReaderTaskManager::waitForAllTasksToFinish()
{
    std::unique_lock<std::mutex> lock(mut_);
    cvAllTaskManagersDone_.wait(lock, [] { return BclReaderTaskManager::numTaskManagers_ == 0; });
}

void BclReaderTaskManager::postExecute()
{
    if (!patternedFlowcellPositions_)
    {
        // nothing to do
        return;
    }

    patternedFlowcellPositions_->waitForReady();
    for (auto& outputBufferForTile : *outputBuffer_)
    {
        outputBufferForTile.positions_ = patternedFlowcellPositions_;
    }
}

BclReader::BclReader(
    TaskQueue& taskQueue,
    ThreadSafeQueue<std::shared_ptr<data::RawBclBufferGroup>>& outputBuffersToSubmit,
    ThreadSafeQueue<std::shared_ptr<data::RawBclBufferGroup>>& outputBuffersToUse,
    const layout::Layout &layout,
    const layout::LaneInfo &laneInfo,
    bool ignoreMissingBcls,
    bool ignoreMissingFilters,
    bool ignoreMissingPositions,
    bool ignoreMissingControls,
    const boost::filesystem::path &inputDir,
    const boost::filesystem::path &intensitiesDir
)
: Stage("Bcl reading", taskQueue)
, outputBuffersToSubmit_(outputBuffersToSubmit)
, outputBuffersToUse_(outputBuffersToUse)
, layout_(layout)
, laneInfo_(laneInfo)
, ignoreMissingBcls_(ignoreMissingBcls)
, ignoreMissingFilters_(ignoreMissingFilters)
, ignoreMissingPositions_(ignoreMissingPositions)
, ignoreMissingControls_(ignoreMissingControls)
, currentTileInfo_(laneInfo_.getTileInfos().begin())
, inputDir_(inputDir)
, intensitiesDir_(intensitiesDir)
, bclFileReaders_()
, positionsFile_(std::make_shared<io::SyncFile>(std::ios_base::in | std::ios_base::binary))
, filterFile_(std::make_shared<io::SyncFile>(std::ios_base::in | std::ios_base::binary))
, patternedFlowcellPositions_(layout.getFlowcellInfo().isPatternedFlowcell() ?
                                  std::make_shared<PatternedPositionsContainer>() :
                                  std::shared_ptr<PatternedPositionsContainer>())
, numTilesPerBuffer_(1)
, numBuffersCreated_(0)
, groupNumber_(0)
, groupSerializer_()
{
    common::CycleNumber cyclesCount = laneInfo.getNumCyclesToLoad();
    BCL2FASTQ_ASSERT_MSG(cyclesCount != 0, "There are no cycles to be processed");

    auto aggregateTilesMode = layout_.getFlowcellInfo().getAggregateTilesMode();

    size_t cycleIndex = 0;
    for (const auto& readInfo : laneInfo.readInfos())
    {
        for (const auto& cycleInfo : readInfo.cyclesToLoad())
        {
            std::shared_ptr<io::SyncFile> bclFile = std::make_shared<io::SyncFile>(std::ios_base::in | std::ios_base::binary);

            switch (aggregateTilesMode)
            {
            case common::TileAggregationMode::AGGREGATED:
                bclFileReaders_.push_back(new data::AggregatedBclFileReader(inputDir,
                                                                            laneInfo,
                                                                            cycleInfo.getNumber(),
                                                                            cycleIndex,
                                                                            ignoreMissingBcls_,
                                                                            bclFile,
                                                                            std::make_shared<data::CycleBCIFile>(ignoreMissingBcls_)));
                break;
            case common::TileAggregationMode::CBCL:
                bclFileReaders_.push_back(new data::CbclFileReader(inputDir,
                                                                   laneInfo,
                                                                   cycleInfo.getNumber(),
                                                                   cycleIndex,
                                                                   ignoreMissingBcls_,
                                                                   layout.getTileFileMap(),
                                                                   bclFile));
                break;
            case common::TileAggregationMode::NON_AGGREGATED:
                bclFileReaders_.push_back(new data::TileBclFileReader(inputDir,
                                                                      laneInfo,
                                                                      cycleInfo.getNumber(),
                                                                      cycleIndex,
                                                                      ignoreMissingBcls_,
                                                                      bclFile));
                break;
            default:
                BCL2FASTQ_ASSERT_MSG(false, "Unexpected aggregation mode type");
            }

            ++cycleIndex;
        }
    }

    calculateNumTilesPerBuffer();
}

void BclReader::calculateNumTilesPerBuffer()
{
    // We want a minimum number of clusters per sample
    size_t desiredClustersCount = data::ClustersPerReadBuffer;
    size_t totalClusters = 0;
    for (auto& tile : laneInfo_.getTileInfos())
    {
        if (!tile.haveClustersCount() || tile.getClustersCount() == 0)
        {
            // Don't have the clusters counts
            return;
        }

        totalClusters += tile.getClustersCount();
    }

    if (totalClusters <= desiredClustersCount)
    {
        numTilesPerBuffer_ = laneInfo_.getTileInfos().size();
        return;
    }

    numTilesPerBuffer_ = laneInfo_.getTileInfos().size() / (totalClusters / desiredClustersCount);
    if (numTilesPerBuffer_ == 0) numTilesPerBuffer_ = 1;
}

std::shared_ptr<BclReaderTaskManager> BclReader::getNewTaskManager()
{
    const std::string longTileId = getTileSpec(laneInfo_, *currentTileInfo_);
    BCL2FASTQ_LOG(common::LogLevel::TRACE) << longTileId << " BclReader::getNewTaskManager 1 " << std::endl;

    std::shared_ptr<data::RawBclBufferGroup> outputBuffer;

    if (outputBuffersToSubmit_.isTerminated())
    {
        // This should only occur if an error was encountered in the BclLoader
        return std::shared_ptr<BclReaderTaskManager>();
    }

    if (numBuffersCreated_ < 3)
    {
        if (!outputBuffersToUse_.tryGetDataCheckEmpty(outputBuffer, 
            longTileId + " Group: " + std::to_string(numBuffersCreated_) + " Push RawBclBufferGroup "))
        {
            // There aren't any buffers available. Make a new one.
            // The memory allocation isn't too bad here, since we're recycling the buffers.
            outputBuffer = std::make_shared<data::RawBclBufferGroup>(numTilesPerBuffer_, laneInfo_.getNumCyclesToLoad());
            ++numBuffersCreated_;
        }
        else
        {
            outputBuffer->waitForFinished();
        }
    }
    else
    {
        if (!outputBuffersToUse_.tryGetData(outputBuffer,
            longTileId + " Group: " + std::to_string(numBuffersCreated_) + " Push RawBclBufferGroup "))
        {
            return std::shared_ptr<BclReaderTaskManager>();
        }
        else
        {
            outputBuffer->waitForFinished();
        }
    }

    outputBuffer->setGroupNumber(groupNumber_++);

    return std::make_shared<BclReaderTaskManager>(ignoreMissingBcls_,
                                                  ignoreMissingFilters_,
                                                  ignoreMissingPositions_,
                                                  ignoreMissingControls_,
                                                  outputBuffer,
                                                  outputBuffersToSubmit_,
                                                  patternedFlowcellPositions_,
                                                  groupSerializer_);
}

void BclReader::terminate()
{
    this->getTaskQueue().terminate();
    outputBuffersToUse_.terminate();    
}

bool BclReader::startNewTasks()
{
    const std::string longTileId = getTileSpec(laneInfo_, *currentTileInfo_);
    BCL2FASTQ_LOG(common::LogLevel::TRACE) << longTileId << " BclReader::startNewTask " << std::endl;

    TaskQueue &taskQueue = this->getTaskQueue();
    if (taskQueue.isTerminated())
    {
        outputBuffersToSubmit_.terminate();

        // Only occurs if there was an error
        return false;
    }

    if (currentTileInfo_ == laneInfo_.getTileInfos().end())
    {
        BclReaderTaskManager::waitForAllTasksToFinish();
        outputBuffersToSubmit_.setFinished();

        return false;
    }

    auto taskManager = getNewTaskManager();

    if (!taskManager)
    {
        return false;
    }

    auto& outputBuffer = taskManager->getOutputBuffer();

    size_t remainingTiles = std::distance(currentTileInfo_, laneInfo_.getTileInfos().end());
    if (remainingTiles < numTilesPerBuffer_)
    {
        outputBuffer.resize(remainingTiles);
    }

    const auto aggregateTilesMode = layout_.getFlowcellInfo().getAggregateTilesMode();
    const bool isPatternedFlowcell = layout_.getFlowcellInfo().isPatternedFlowcell();

    auto tileInfoForBuffers = currentTileInfo_;
    for (auto& buffer : outputBuffer)
    {
        buffer.cycleData_.resize(laneInfo_.getNumCyclesToLoad());
        buffer.gzipDecompressors_.resize(buffer.cycleData_.size());
        buffer.tileInfo_ = tileInfoForBuffers;
        ++tileInfoForBuffers;
    }

    // Add the positions task first, since if it's a patterned flowcell, we
    // only need to load it once, but all tiles will need the data.
    if (!isPatternedFlowcell || numBuffersCreated_ == 1)
    {
        taskQueue.addData(std::make_shared<PositionsReadTask>(
            taskManager,
            positionsFile_,
            intensitiesDir_,
            data::Instrument::isPositionFileAggregated(aggregateTilesMode),
            //aggregateTilesMode == common::TileAggregationMode::AGGREGATED ||
            //    aggregateTilesMode == common::TileAggregationMode::CBCL,
            isPatternedFlowcell,
            laneInfo_,
            currentTileInfo_,
            ignoreMissingPositions_,
            outputBuffer,
            patternedFlowcellPositions_
        ), longTileId + " BclReader Task: Position read");
    }

        BclFilesContainer::size_type bclIdx = 0;
        for (const auto& readInfo : laneInfo_.readInfos())
        {
            // We want all the cycles for the index reads, even if they were masked.
            for (const auto& cycleInfo : readInfo.cyclesToLoad())
            {
                taskQueue.addData(std::make_shared<BclReadTask>(
                    taskManager,
                    bclFileReaders_.at(bclIdx),
                    outputBuffer
                ));
                //), longTileId + " Cycle: " + std::to_string(cycleInfo.getNumber()) + "  Read: " + std::to_string(readInfo.getNumber()) + " Task: bcl-read");

                ++bclIdx;
            }
        }

        taskQueue.addData(std::make_shared<FilterReadTask>(
            taskManager,
            filterFile_,
            inputDir_,
            aggregateTilesMode,
            laneInfo_,
            currentTileInfo_,
            ignoreMissingFilters_,
            outputBuffer
        ), longTileId + " Task: filter-read");

        taskQueue.addData(std::make_shared<ControlReadTask>(
            taskManager,
            inputDir_,
            laneInfo_,
            currentTileInfo_,
            ignoreMissingControls_,
            outputBuffer
        ));

    currentTileInfo_ += (remainingTiles < numTilesPerBuffer_ ? remainingTiles : numTilesPerBuffer_);

    return true;
}

} // namespace conversion
} // namespace bcl2fastq

