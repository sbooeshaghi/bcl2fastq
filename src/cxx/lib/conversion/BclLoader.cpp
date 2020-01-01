/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file BclLoader.cpp
 *
 * \brief Implementation of BCL loader.
 *
 * \author Marek Balint
 * \author Mauricio Varea
 */


#include <errno.h>
#include <algorithm>
#include <utility>
#include <chrono>
#include <thread>

#include <boost/format.hpp>
#include <boost/filesystem.hpp>

#include "common/Debug.hh"
#include "common/Exceptions.hh"
#include "common/Types.hh"
#include "conversion/BclLoader.hh"
#include "data/Instrument.hh"


namespace bcl2fastq {


namespace conversion {

typedef std::shared_ptr<BclLoaderTaskManager> BclLoaderTaskManagerPtr;

namespace
{
    std::string getTileSpec(
        const layout::LaneInfo &laneInfo)
    {
        std::stringstream ss;
        ss << "Lane: " << laneInfo.getNumber() << " Tile: xxxx";

        return ss.str();
    }

    std::string getTileSpec(
        const layout::LaneInfo &laneInfo,
        const layout::TileInfo &tileInfo)
    {
        std::stringstream ss;
        ss << "Lane: " << laneInfo.getNumber() << " Tile: " << tileInfo.getNumber();

        return ss.str();
    }
}

BclLoadTask::BclLoadTask(
    std::shared_ptr<BclLoaderTaskManager>& taskManager,
    data::RawBclBufferGroup& bclData,
    std::mutex* aggregatedBclFileMutex,
    common::TileAggregationMode aggregateTilesMode,
    const layout::LaneInfo &laneInfo,
    const layout::CycleInfo &cycleInfo,
    layout::LaneInfo::TileInfosContainer::const_iterator tileInfoIter,
    bool ignoreMissingBcls,
    size_t bclIdx,
    DemuxBuffer& outputBuffer
)
: Task(taskManager)
, bclData_(bclData)
, aggregatedBclFileMutex_(aggregatedBclFileMutex)
, aggregateTilesMode_(aggregateTilesMode)
, laneInfo_(laneInfo)
, cycleInfo_(cycleInfo)
, tileInfoIter_(tileInfoIter)
, ignoreMissingBcls_(ignoreMissingBcls)
, bclIdx_(bclIdx)
, outputBuffer_(outputBuffer)
{
}

bool BclLoadTask::execute(int32_t)
{
    switch (aggregateTilesMode_)
    {
    case common::TileAggregationMode::AGGREGATED:
        return readAggregatedTiles();
        break;
    case common::TileAggregationMode::NON_AGGREGATED:
        return readSeparateTiles();
        break;
    case common::TileAggregationMode::CBCL:
        return readCbclTiles();
        break;
    default:
        BCL2FASTQ_ASSERT_MSG(false, "Unexpected aggregation mode type");
        break;
    }

    return false;
}

bool BclLoadTask::readSeparateTiles()
{
    static int sequentialId = 0; ++sequentialId;

    size_t tileIdx = 0;
    for (auto& rawDataForTile : bclData_)
    {
        data::BclFile bclFile(rawDataForTile.cycleData_[bclIdx_].bcls_,
                              rawDataForTile.gzipDecompressors_[bclIdx_],
                              ignoreMissingBcls_,
                              true);

        auto clustersCount = bclFile.getClustersCount();

        common::resizeAndRealloc(outputBuffer_.bclBuffers_[tileIdx].bcls_[bclIdx_].bcls_, clustersCount);
        auto bytesRead = bclFile.read(outputBuffer_.bclBuffers_[tileIdx].bcls_[bclIdx_].bcls_.data(), clustersCount);
        if (bytesRead != clustersCount)
        {
            int errnum = errno;
            if (ignoreMissingBcls_)
            {
                // Assume the truncated file is corrupt.
                std::stringstream ss;
                ss << "BCL file '" << bclFile.getPath() << "' ";
                if (bytesRead > 0) {
                    ss << "corrupt: " << std::endl;
                }
                else {
                    ss << "truncated: " << std::endl;
                }
                ss << "bytes_read=" << bytesRead << " bytes_expected=" << clustersCount;

                BCL2FASTQ_LOG(common::LogLevel::WARNING) << ss.str() << std::endl;

                bytesRead = 0;
                outputBuffer_.bclBuffers_[tileIdx].bcls_[bclIdx_].bcls_.resize(bytesRead > 0 ? bytesRead : 0);
            }
            else
            {
                BOOST_THROW_EXCEPTION(common::InputDataError(errnum, (boost::format("BCL file '%s' truncated: bytes_read=%d bytes_expected=%d") % bclFile.getPath().string() % bytesRead % clustersCount).str()));
            }
        }

        ++tileIdx;
    }

    return true;
}

bool BclLoadTask::readAggregatedTiles()
{
    std::lock_guard<std::mutex> lock(*aggregatedBclFileMutex_);

    int32_t tileIdx = 0;
    auto rawDataIter = bclData_.begin();
    while (rawDataIter != bclData_.end())
    {
        common::ClustersCount clustersCount = tileInfoIter_->getClustersCount();

        if (clustersCount != 0)
        {
            std::vector<char> &bclVector = outputBuffer_.bclBuffers_[tileIdx].bcls_[bclIdx_].bcls_;
            common::resizeAndRealloc(bclVector, clustersCount);

            std::shared_ptr<data::BclFile> bclFile(std::make_shared<data::BclFile>(
                rawDataIter->cycleData_[bclIdx_].bcls_,
                rawDataIter->gzipDecompressors_[bclIdx_],
                ignoreMissingBcls_,
                tileInfoIter_->getIndex() == 0));

            std::streamsize bytesRead = 0;
            if (bclFile->seek(0, rawDataIter->uncompressedBclOffset_))
            {
                bytesRead = bclFile->read(bclVector.data(), clustersCount);
            }
            else
            {
                // Explicitly clear the empty tile.
                outputBuffer_.bclBuffers_[tileIdx].bcls_[bclIdx_].bcls_.resize(0);
            }

            checkBytesRead(bytesRead,
                           clustersCount,
                           tileIdx,
                           rawDataIter->cycleData_[bclIdx_].bcls_.path_);
        }
        else
        {
            // Explicitly clear the empty tile.
            outputBuffer_.bclBuffers_[tileIdx].bcls_[bclIdx_].bcls_.resize(0);
        }

        ++rawDataIter;
        ++tileInfoIter_;
        ++tileIdx;
    }

    return true;
}

std::stringstream &BclLoadTask::getTileDescriptor(std::stringstream &msg, data::PerCycleData perCycleData)
{
    if (perCycleData.bcls_.path_.string().empty()) {
        msg << "BCL data for lane " << laneInfo_.getNumber()
            << " cycle " << cycleInfo_.getNumber()
            << " tile " << tileInfoIter_->getNumber();
    } else {
        msg << "Tile: " << tileInfoIter_->getNumber() << ", BCL file '" << perCycleData.bcls_.path_ << "'";
    }

    return msg;
}

bool BclLoadTask::readCbclTiles()
{
    size_t tileIdx = 0;
    for (data::RawBclBuffer& rawDataForTile : bclData_)
    {
        auto& gzipDecompressor = rawDataForTile.gzipDecompressors_[bclIdx_];

        gzipDecompressor.flush();
        gzipDecompressor.reset();

        auto& cycleData = rawDataForTile.cycleData_[bclIdx_];

        boost::iostreams::basic_array_source<char> inputSrc(cycleData.bcls_.data(),
                                                            cycleData.bcls_.size());
        boost::iostreams::stream<boost::iostreams::basic_array_source<char>> istr(inputSrc);

        outputBuffer_.bclBuffers_[tileIdx].bcls_[bclIdx_].bcls_.resize(cycleData.uncompressedBlockSize_);
        outputBuffer_.bclBuffers_[tileIdx].bcls_[bclIdx_].includeNonPf_ = cycleData.includeNonPf_;
        outputBuffer_.bclBuffers_[tileIdx].bcls_[bclIdx_].numBitsPerQscore_ = cycleData.numBitsPerQscore_;
        outputBuffer_.bclBuffers_[tileIdx].bcls_[bclIdx_].remappedQscores_ = cycleData.remappedQscores_;

        std::streamsize bytesRead = 0;
        bool decompressed = true;
        try
        {
            bytesRead = gzipDecompressor.read(istr,
                                              outputBuffer_.bclBuffers_[tileIdx].bcls_[bclIdx_].bcls_.data(),
                                              outputBuffer_.bclBuffers_[tileIdx].bcls_[bclIdx_].bcls_.size());
        }
        catch (...)
        {
            bytesRead = 0;
            decompressed = false;

            std::stringstream msg;
            getTileDescriptor(msg, cycleData);
            msg << " corrupt.";

            if (ignoreMissingBcls_)
            {
                BCL2FASTQ_LOG(common::LogLevel::WARNING) << msg.str() << std::endl;
                outputBuffer_.bclBuffers_[tileIdx].bcls_[bclIdx_].bcls_.resize(0);
            }
            else
            {
                int errnum = errno;
                BOOST_THROW_EXCEPTION(common::InputDataError(errnum, msg.str()));
            }
        }

        if (decompressed && bytesRead != cycleData.uncompressedBlockSize_)
        {
            std::stringstream msg;
            getTileDescriptor(msg, cycleData);
            msg << " truncated: bytes_read=" << bytesRead << " bytes_expected=" << cycleData.uncompressedBlockSize_;

            int errnum = errno;
            if (ignoreMissingBcls_)
            {
                BCL2FASTQ_LOG(common::LogLevel::WARNING) << msg.str() << std::endl;
                outputBuffer_.bclBuffers_[tileIdx].bcls_[bclIdx_].bcls_.resize(bytesRead > 0 ? bytesRead : 0);
            }
            else
            {
                BOOST_THROW_EXCEPTION(common::InputDataError(errnum, msg.str()));
            }
        }

        ++tileIdx;
    }

    return true;
}

void BclLoadTask::checkBytesRead(int32_t bytesRead,
                                 int32_t expectedBytes,
                                 int32_t tileIdx,
                                 const boost::filesystem::path& bclPath)
{
        if (bytesRead < expectedBytes)
        {
            int errnum = errno;
            if (ignoreMissingBcls_) {
                if (!bclPath.string().empty()) {
                    BCL2FASTQ_LOG(common::LogLevel::WARNING) << "BCL file '" << bclPath << "' truncated: bytes_read=" << bytesRead << " bytes_expected=" << expectedBytes << ":" << std::strerror(errnum) << " (" << errnum << ")" << std::endl;
                }
                outputBuffer_.bclBuffers_[tileIdx].bcls_[bclIdx_].bcls_.resize(bytesRead > 0 ? bytesRead : 0);
            }
            else
            {
                BOOST_THROW_EXCEPTION(common::InputDataError(errnum, (boost::format("BCL file '%s' truncated: bytes_read=%d bytes_expected=%d") % bclPath.string() % bytesRead % expectedBytes).str()));
            }
        }
}

PositionsLoadTask::PositionsLoadTask(
    std::shared_ptr<BclLoaderTaskManager>& taskManager,
    data::RawBclBufferGroup& bclData,
    bool aggregateTilesFlag,
    const layout::LaneInfo &laneInfo,
    layout::LaneInfo::TileInfosContainer::const_iterator &tileInfoIter,
    const std::vector<common::ClustersCount>& clustersCounts,
    bool ignoreMissingPositions,
    data::BclBufferVec& outputBuffer,
    std::shared_ptr<data::BclBuffer::PatternedPositionsContainer>& patternedFlowcellPositions
)
: Task(taskManager)
, bclData_(bclData)
, aggregateTilesFlag_(aggregateTilesFlag)
, laneInfo_(laneInfo)
, tileInfoIter_(tileInfoIter)
, clustersCounts_(clustersCounts)
, ignoreMissingPositions_(ignoreMissingPositions)
, outputBuffer_(outputBuffer)
, patternedFlowcellPositions_(patternedFlowcellPositions)
{
}

bool PositionsLoadTask::execute(int32_t)
{
    if (patternedFlowcellPositions_)
    {
        if (execute(*bclData_.begin()->positions_,
                    *patternedFlowcellPositions_))
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
        size_t bufferIndex = 0;
        for (auto& inputBuffer : bclData_)
        {
            auto& outputBuffer = outputBuffer_[bufferIndex].positions_;
            if (!outputBuffer)
            {
                outputBuffer = std::make_shared<data::BclBuffer::PositionsContainer>();
            }

            execute(*inputBuffer.positions_,
                    *outputBuffer);

            ++bufferIndex;
            ++tileInfoIter_;
        }
    }

    return true;
}

bool PositionsLoadTask::execute(const common::RawDataBuffer& inputBuffer,
                                data::BclBuffer::PositionsContainer& outputBuffer)
{
    bool skipParsingHeader = !patternedFlowcellPositions_ && (aggregateTilesFlag_ && (tileInfoIter_->getIndex() != 0));
    auto bytesRead = data::PositionsFileFactory::read(inputBuffer,
                                                      ignoreMissingPositions_,
                                                      skipParsingHeader,
                                                      outputBuffer,
                                                      tileInfoIter_->getClustersCount());

    return true;
}


FilterLoadTask::FilterLoadTask(
    std::shared_ptr<BclLoaderTaskManager>& taskManager,
    data::RawBclBufferGroup& bclData,
    common::TileAggregationMode tileAggregationMode,
    const layout::LaneInfo &laneInfo,
    layout::LaneInfo::TileInfosContainer::const_iterator &tileInfoIter,
    const std::vector<common::ClustersCount>& clustersCounts,
    bool ignoreMissingFilters,
    data::BclBufferVec &outputBuffer
)
: Task(taskManager)
, bclData_(bclData)
, tileAggregationMode_(tileAggregationMode)
, aggregateTilesFlag_(data::Instrument::isFilterFileAggregated(tileAggregationMode))
, laneInfo_(laneInfo)
, tileInfoIter_(tileInfoIter)
, clustersCounts_(clustersCounts)
, ignoreMissingFilters_(ignoreMissingFilters)
, outputBuffer_(outputBuffer)
{
}

bool FilterLoadTask::execute(int32_t)
{
    size_t bufferIndex = 0;
    for (auto& rawDataForTile : bclData_)
    {
        readFilterFile(rawDataForTile.filters_,
                       *tileInfoIter_,
                       outputBuffer_[bufferIndex],
                       aggregateTilesFlag_,
                       aggregateTilesFlag_ && (tileInfoIter_->getIndex() != 0));

        ++bufferIndex;
        ++tileInfoIter_;
    }

    return true;
}

void FilterLoadTask::readFilterFile(common::RawDataBuffer& filterData,
                                    const layout::TileInfo& tileInfo,
                                    data::BclBuffer& outputBuffer,
                                    bool isAggregatedTiles,
                                    bool skipHeader)
{
    data::FilterFile filterFile(filterData,
                                ignoreMissingFilters_,
                                skipHeader);

    auto clustersCount = isAggregatedTiles ? tileInfo.getClustersCount() : filterFile.getClustersCount();

    common::resizeAndRealloc(outputBuffer.filters_, clustersCount);
    auto recordsRead = filterFile.read(outputBuffer.filters_.data(), clustersCount);

    validateFilterRecords(recordsRead,
                          clustersCount,
                          filterData.path_,
                          outputBuffer);
}

void FilterLoadTask::validateFilterRecords(std::size_t                        recordsRead,
                                           common::ClustersCount              clustersCount,
                                           const boost::filesystem::path&     filePath,
                                           data::BclBuffer& outputBuffer)
{
    int errnum = errno;
    if (recordsRead != clustersCount)
    {
        if (ignoreMissingFilters_)
        {
            if (tileAggregationMode_ == common::TileAggregationMode::CBCL)
            {
                // If cbcl filter file is corrupt then we can't map the tile indices so ignore this tile.
                BCL2FASTQ_LOG(common::LogLevel::WARNING) << "CBCL filter file '"
                        << filePath << "' corrupt: this tile will be ignored." << std::endl;

                outputBuffer.filters_.resize(0);

                for (data::PerCycleData &perCycleData : outputBuffer.bcls_)
                {
                    perCycleData.bcls_.resize(0);
                    perCycleData.uncompressedBlockSize_ = 0;
                }
            }
            else
            {
                BCL2FASTQ_LOG(common::LogLevel::WARNING) << "Filter file '" << filePath << "' truncated: records_read=" << recordsRead << " records_expected=" << clustersCount << ":" << std::strerror(errnum) << " (" << errnum << ")" << std::endl;
                outputBuffer.filters_.resize(recordsRead);
                static const data::FilterFile::Record defaultRecord(0x01);
                outputBuffer.filters_.resize(clustersCount, defaultRecord);
            }

        }
        else
        {
            BOOST_THROW_EXCEPTION(common::InputDataError(errnum, (boost::format("Filter file '%s' truncated: records_real=%d records_expected=%d") % filePath.string() % recordsRead % clustersCount).str()));
        }
    }
}

ControlLoadTask::ControlLoadTask(
    std::shared_ptr<BclLoaderTaskManager>& taskManager,
    data::RawBclBufferGroup& bclData,
    const layout::LaneInfo &laneInfo,
    bool ignoreMissingControls,
    data::BclBufferVec &outputBuffer
)
: Task(taskManager)
, bclData_(bclData)
, laneInfo_(laneInfo)
, ignoreMissingControls_(ignoreMissingControls)
, outputBuffer_(outputBuffer)
{
}

bool ControlLoadTask::execute(int32_t)
{

    int32_t bufferIndex = -1;
    for (auto& rawDataForTile : bclData_)
    {
        ++bufferIndex;

        if (rawDataForTile.controls_.empty())
        {
            outputBuffer_[bufferIndex].controls_.clear();

            // There might not be a control file
            continue;
        }

        data::ControlFile controlFile(rawDataForTile.controls_,
                                      ignoreMissingControls_);

        common::ClustersCount clustersCount = controlFile.getClustersCount();

        common::resizeAndRealloc(outputBuffer_[bufferIndex].controls_, clustersCount);
        std::size_t recordsRead = controlFile.read(&*outputBuffer_[bufferIndex].controls_.begin(), clustersCount);

        int errnum = errno;
        if (recordsRead != clustersCount)
        {
            if (ignoreMissingControls_)
            {
                BCL2FASTQ_LOG(common::LogLevel::WARNING) << "Control file '" << rawDataForTile.controls_.path_ << "' truncated: records_read=" << recordsRead << " records_expected=" << clustersCount << ":" << std::strerror(errnum) << " (" << errnum << ")" << std::endl;
                outputBuffer_[bufferIndex].controls_.resize(recordsRead);
            }
            else
            {
                BOOST_THROW_EXCEPTION(common::InputDataError(errnum, (boost::format("Control file '%s' truncated: records_real=%d records_expected=%d") % rawDataForTile.controls_.path_.string() % recordsRead % clustersCount).str()));
            }
        }
    }

    return true;
}

std::atomic<uint32_t> BclLoaderTaskManager::numTaskManagers_(0);
std::condition_variable BclLoaderTaskManager::cvAllTaskManagersDone_;
std::mutex BclLoaderTaskManager::mut_;

BclLoaderTaskManager::BclLoaderTaskManager(bool ignoreMissingBcls,
                                           bool ignoreMissingFilters,
                                           bool ignoreMissingPositions,
                                           bool ignoreMissingControls,
                                           common::TileAggregationMode aggregateTilesFlag,
                                           std::shared_ptr<data::RawBclBufferGroup>& inputBuffer,
                                           ThreadSafeQueue<std::shared_ptr<data::RawBclBufferGroup>>& inputBuffersToRecycle,
                                           std::shared_ptr<DemuxBuffer>& outputBuffer,
                                           ThreadSafeQueue<std::shared_ptr<DemuxBuffer>>& outputBufferQueue,
                                           std::shared_ptr<data::BclBuffer::PatternedPositionsContainer>& patternedFlowcellPositions,
                                           TaskManager::GroupSerializer &groupSerializer)
: TaskManager()
, ignoreMissingBcls_(ignoreMissingBcls)
, ignoreMissingFilters_(ignoreMissingFilters)
, ignoreMissingPositions_(ignoreMissingPositions)
, ignoreMissingControls_(ignoreMissingControls)
, aggregateTilesFlag_(aggregateTilesFlag)
, uniqueFailedReadIndex_(0)
, inputBuffer_(inputBuffer)
, inputBuffersToRecycle_(inputBuffersToRecycle)
, outputBuffer_(outputBuffer)
, outputBufferQueue_(outputBufferQueue)
, patternedFlowcellPositions_(patternedFlowcellPositions)
, groupSerializer_(groupSerializer)
{
    ++numTaskManagers_;
}

BclLoaderTaskManager::~BclLoaderTaskManager()
{
    auto prevInputBuffer = inputBuffer_->getPrevBuffer();
    if (prevInputBuffer)
    {
        prevInputBuffer->waitForFinished();
        prevInputBuffer->reset();
        inputBuffersToRecycle_.addData(prevInputBuffer);
    }

    // All tasks are complete.
    // Never touch this buffer again. It could be recycled after this point.
    inputBuffer_->setFinished();

    BCL2FASTQ_LOG(common::LogLevel::DEBUG) << "dQ Raw Tile: " << inputBuffer_->begin()->tileInfo_->getNumber() << std::endl;
    BCL2FASTQ_LOG(common::LogLevel::DEBUG) << " Q Demux Tile: " << outputBuffer_->bclBuffers_.begin()->tileInfo_->getNumber() << std::endl;

    postExecute(*outputBuffer_);

    // Get ready to set the demux info
    for (auto& buffer : outputBuffer_->bclBuffers_)
    {
        buffer.samples_.resize(buffer.bcls_.front().getNumClusters());
    }

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

void BclLoaderTaskManager::waitForAllTasksToFinish()
{
    std::unique_lock<std::mutex> lock(mut_);
    cvAllTaskManagersDone_.wait(lock, [] { return BclLoaderTaskManager::numTaskManagers_ == 0; });
}

void BclLoaderTaskManager::validateCycleMetadata(data::PerCycleData& cycleData,
                                                 data::PerCycleData& prevCycleData)
{
    if (prevCycleData.remappedQscores_.empty() || cycleData.remappedQscores_.empty())
    {
        return;
    }

    if (cycleData.numBitsPerQscore_ != prevCycleData.numBitsPerQscore_ ||
        cycleData.remappedQscores_ != prevCycleData.remappedQscores_ ||
        (cycleData.includeNonPf_ && !prevCycleData.includeNonPf_))
    {
        BOOST_THROW_EXCEPTION(common::InputDataError("Corrupt cbcl header detected."));
    }
}

void BclLoaderTaskManager::postExecute(DemuxBuffer& postExBuffer)
{

    for (data::BclBuffer& outputBuffer : postExBuffer.bclBuffers_)
    {
        if (patternedFlowcellPositions_)
        {
            patternedFlowcellPositions_->waitForReady();
            if (outputBufferQueue_.isTerminated())
            {
                return;
            }

            outputBuffer.positions_ = patternedFlowcellPositions_;
        }

        const auto largestBuffer = std::max_element(
            outputBuffer.bcls_.begin(),
            outputBuffer.bcls_.end(),
            [] (const data::PerCycleData& a, const data::PerCycleData& b)
                { return a.bcls_.size() < b.bcls_.size(); }
        );

        const auto bufferSizeNonFiltered = largestBuffer->bcls_.size();

        // Calculate the number of clusters (could be more than 1 per byte)
        const auto clustersCount = bufferSizeNonFiltered * largestBuffer->getNumClustersPerByte();

        const size_t largestFilteredSize = std::max_element(
            outputBuffer.bcls_.begin(),
            outputBuffer.bcls_.end(),
            [] (const data::PerCycleData& a, const data::PerCycleData& b)
                {
                    if (!a.includeNonPf_ && !b.includeNonPf_) return a.bcls_.size() < b.bcls_.size();
                    if (!a.includeNonPf_) return false;
                    if (!b.includeNonPf_) return true;
                    return false;
                }
        )->bcls_.size();

        common::CycleNumber cycleNumber = 1;
        for (auto &perCycleData : outputBuffer.bcls_)
        {
            if (cycleNumber > 1)
            {
                validateCycleMetadata(perCycleData, outputBuffer.bcls_[cycleNumber-2]);
            }

            auto& bclBuffer = perCycleData.bcls_;

            const auto currentSize = bclBuffer.size();
            const auto bufferSize = perCycleData.includeNonPf_ ? bufferSizeNonFiltered : largestFilteredSize;
            if (bufferSize != currentSize)
            {
                bclMismatchCount("BCL",cycleNumber,outputBuffer.tileInfo_->getNumber(),currentSize,bufferSize);
                common::resizeAndRealloc(bclBuffer, bufferSize, 0);
            }
            ++cycleNumber;
        }
        if (outputBuffer.positions_ && outputBuffer.positions_->size() != clustersCount)
        {
            const std::string type = (!ignoreMissingPositions_ && ignoreMissingBcls_) ? "BCL" : "positions";
            bclMismatchCount(type,1,outputBuffer.tileInfo_->getNumber(), outputBuffer.positions_->size(), clustersCount);
            createUniqueFakePositions(outputBuffer, clustersCount);
        }
        if (outputBuffer.filters_.size() != clustersCount)
        {
            // For cbcl files, mark reads as failed filter if the filter file is missing.
            // We can't use this data.
            const std::string type = (!ignoreMissingFilters_ && ignoreMissingBcls_) ? "BCL" : "filter";
            bclMismatchCount(type,1,outputBuffer.tileInfo_->getNumber(),outputBuffer.filters_.size(), clustersCount);
            common::resizeAndRealloc(outputBuffer.filters_,
                                     clustersCount,
                                     data::FilterFile::Record(aggregateTilesFlag_ == common::TileAggregationMode::CBCL ? 0 : 1));
        }
        if (outputBuffer.controls_.size() != clustersCount)
        {
            const std::string type = (!ignoreMissingControls_ && ignoreMissingBcls_) ? "BCL" : "control";
            if (!outputBuffer.controls_.empty())
            {
                bclMismatchCount(type,1,outputBuffer.tileInfo_->getNumber(),outputBuffer.controls_.size(),clustersCount);
            }
            common::resizeAndRealloc(outputBuffer.controls_, clustersCount, 0);
        }
    }
}

void BclLoaderTaskManager::bclMismatchCount(const std::string& fileType,
                                            common::CycleNumber cycleNumber,
                                            common::TileNumber tileNumber,
                                            common::RawDataBuffer::size_type realSize,
                                            common::RawDataBuffer::size_type expectedSize)
{
    if ((ignoreMissingBcls_ && "BCL" == fileType)
    ||  (ignoreMissingPositions_ && "positions" == fileType)
    ||  (ignoreMissingFilters_ && "filter" == fileType)
    ||  (ignoreMissingControls_ && "control" == fileType))
    {
        BCL2FASTQ_LOG(common::LogLevel::WARNING) << "Mismatching cluster count in " << fileType << " file: Cycle #" << cycleNumber << ", Tile#" << tileNumber << ": bytes_read=" << realSize << " bytes_expected=" << expectedSize << std::endl;
    } else {
        BOOST_THROW_EXCEPTION(common::InputDataError((boost::format("Mismatching cluster count in %s file: Cycle #%d, Tile#%d: bytes_read=%d bytes_expected=%d") % fileType % cycleNumber % tileNumber % realSize % expectedSize).str()));
    }
}

void BclLoaderTaskManager::createUniqueFakePositions(data::BclBuffer& outputBuffer,
                                                     data::PerCycleData::BclsContainer::size_type bufferSize)
{
    common::resizeAndRealloc(*outputBuffer.positions_, bufferSize);

    for (data::PerCycleData::BclsContainer::size_type i = 0; i < bufferSize; ++i)
    {
        (*outputBuffer.positions_)[i].y_ = uniqueFailedReadIndex_++;
    }

}


BclLoader::BclLoader(
    TaskQueue& taskQueue,
    ThreadSafeQueue<std::shared_ptr<data::RawBclBufferGroup>>& inputBuffersToUse,
    ThreadSafeQueue<std::shared_ptr<data::RawBclBufferGroup>>& inputBuffersToRecycle,
    ThreadSafeQueue<std::shared_ptr<DemuxBuffer>>& outputBuffersToSubmit,
    ThreadSafeQueue<std::shared_ptr<DemuxBuffer>>& outputBuffersToUse,
    const layout::Layout &layout,
    const layout::LaneInfo &laneInfo,
    bool ignoreMissingBcls,
    bool ignoreMissingFilters,
    bool ignoreMissingPositions,
    bool ignoreMissingControls,
    uint32_t numTilesPerBuffer
)
: Stage("Bcl loading", taskQueue)
, inputBuffersToUse_(inputBuffersToUse)
, inputBuffersToRecycle_(inputBuffersToRecycle)
, outputBuffersToSubmit_(outputBuffersToSubmit)
, outputBuffersToUse_(outputBuffersToUse)
, layout_(layout)
, laneInfo_(laneInfo)
, ignoreMissingBcls_(ignoreMissingBcls)
, ignoreMissingFilters_(ignoreMissingFilters)
, ignoreMissingPositions_(ignoreMissingPositions)
, ignoreMissingControls_(ignoreMissingControls)
, currentTileInfo_(laneInfo_.getTileInfos().begin())
, aggregatedBclMutexes_()
, cycleBciFiles_()
, patternedFlowcellPositions_(layout.getFlowcellInfo().isPatternedFlowcell() ?
                                  std::make_shared<data::BclBuffer::PatternedPositionsContainer>() :
                                  std::shared_ptr<data::BclBuffer::PatternedPositionsContainer>())
, clustersCounts_()
, numBuffersProcessed_(0)
, numBuffersCreated_(0)
, numTilesPerBuffer_(numTilesPerBuffer)
, prevInputBuffer_()
, terminated_(false)
, groupSerializer_()
{
    common::CycleNumber cyclesCount = laneInfo.getNumCyclesToLoad();
    BCL2FASTQ_ASSERT_MSG(cyclesCount != 0, "There are no cycles to be processed");

    for (common::CycleNumber i = 0; i < cyclesCount; ++i)
    {
        aggregatedBclMutexes_.push_back(new std::mutex);
        cycleBciFiles_.push_back(new data::CycleBCIFile(ignoreMissingBcls_));
    }
}

void BclLoader::terminate()
{
    this->getTaskQueue().terminate();
    inputBuffersToUse_.terminate();

    outputBuffersToSubmit_.terminate();

    if (patternedFlowcellPositions_)
    {
        patternedFlowcellPositions_->setReady();
    }

    terminated_ = true;
}

std::shared_ptr<BclLoaderTaskManager> BclLoader::getNewTaskManager()
{
    std::string longTileId = getTileSpec(laneInfo_);
    BCL2FASTQ_LOG(common::LogLevel::TRACE) << longTileId << " BclLoader::getNewTaskManager " << std::endl;

    std::shared_ptr<data::RawBclBufferGroup> inputBuffer;
    if (!inputBuffersToUse_.tryGetData(inputBuffer))
    {
        if (!terminated_)
        {
            BclLoaderTaskManager::waitForAllTasksToFinish();
        }

        outputBuffersToSubmit_.setFinished();

        // Only get here if we're done. Otherwise, tryGetData would wait.
        return std::shared_ptr<BclLoaderTaskManager>();
    }

    layout::LaneInfo::TileInfosContainer::const_iterator tileInfoIter = inputBuffer->begin()->tileInfo_;
    longTileId = getTileSpec(laneInfo_, *tileInfoIter);

    std::shared_ptr<DemuxBuffer> outputBuffer;
    if (numBuffersCreated_ < 3)
    {
        if (!outputBuffersToUse_.tryGetDataCheckEmpty(outputBuffer,
            longTileId + " Group: " + std::to_string(numBuffersCreated_) + " Push DemuxBuffer"))
        {
            // There aren't any buffers available. Make a new one.
            // The memory allocation isn't too bad here, since we're recycling the buffers.

            outputBuffer = std::make_shared<DemuxBuffer>(laneInfo_.getSampleInfos().size(),
                                                         data::ClustersPerTask,
                                                         laneInfo_.getNumCyclesToLoad());
            ++numBuffersCreated_;
        }
    }
    else
    {
        if (!outputBuffersToUse_.tryGetData(outputBuffer,
            longTileId + " Group: " + std::to_string(numBuffersProcessed_) + " Push DemuxBuffer"))
        {
            return std::shared_ptr<BclLoaderTaskManager>();
        }
    }

    outputBuffer->setGroupNumber(inputBuffer->getGroupNumber());

    // Set previous buffer
    inputBuffer->setPrevBuffer(prevInputBuffer_);
    prevInputBuffer_ = inputBuffer;

    ++numBuffersProcessed_;
    return std::make_shared<BclLoaderTaskManager>(ignoreMissingBcls_,
                                                  ignoreMissingFilters_,
                                                  ignoreMissingPositions_,
                                                  ignoreMissingControls_,
                                                  layout_.getFlowcellInfo().getAggregateTilesMode(),
                                                  inputBuffer,
                                                  inputBuffersToRecycle_,
                                                  outputBuffer,
                                                  outputBuffersToSubmit_,
                                                  patternedFlowcellPositions_,
                                                  groupSerializer_);
}

bool BclLoader::startNewTasks()
{
    std::string longTileId = getTileSpec(laneInfo_);
    BCL2FASTQ_LOG(common::LogLevel::TRACE) << longTileId << " BclLoader::startNewTask " << std::endl;

    TaskQueue &taskQueue = this->getTaskQueue();
    if (taskQueue.isTerminated() || outputBuffersToSubmit_.isTerminated())
    {
        // Only occurs if there was an error
        terminate();

        return false;
    }

    // TODO: Tidy. currentTileInfo_ is not guaranteed to define the correct tile spec. It is used as a counter here.
    if (currentTileInfo_ == laneInfo_.getTileInfos().end())
    {
        BclLoaderTaskManager::waitForAllTasksToFinish();
        outputBuffersToSubmit_.setFinished();

        return false;
    }

    // Does clustersCounts_ always work here?
    BclLoaderTaskManagerPtr taskManager = getNewTaskManager();

    if (!taskManager)
    {
        return false;
    }

    auto allOutputBuffers = taskManager->getOutputBuffer();
    auto& outputBuffer = allOutputBuffers->bclBuffers_;
    data::RawBclBufferGroup& inputBuffer = *taskManager->getInputBuffer();

    layout::LaneInfo::TileInfosContainer::const_iterator myCurrentTileInfo = inputBuffer.begin()->tileInfo_;
    uint32_t numTiles = std::min(numTilesPerBuffer_, (uint32_t)std::distance(inputBuffer.begin(), inputBuffer.end()));
    outputBuffer.resize(numTiles);

    longTileId = getTileSpec(laneInfo_, *myCurrentTileInfo);

    layout::LaneInfo::TileInfosContainer::const_iterator tileInfoIter = myCurrentTileInfo;

    size_t bufferNum = 0;
    for (auto& buffer : outputBuffer)
    {
        buffer.bcls_.resize(laneInfo_.getNumCyclesToLoad());
        buffer.tileInfo_ = myCurrentTileInfo + bufferNum;
        ++bufferNum;
    }

    const auto aggregateTilesMode = layout_.getFlowcellInfo().getAggregateTilesMode();
    const bool isPatternedFlowcell = layout_.getFlowcellInfo().isPatternedFlowcell();

    clustersCounts_.clear();

    // Add the positions task first, since if it's a patterned flowcell, we
    // only need to load it once, but all tiles will need the data.
    const bool addPositionsTask = 
        (!isPatternedFlowcell && inputBuffer.hasCycleData()) || // only if data is present for non-pattern flow cells. 
        (numBuffersProcessed_ == 1); // only once for patternflow cell
    if (addPositionsTask)
    {
        taskQueue.addData(std::make_shared<PositionsLoadTask>(
            taskManager,
            inputBuffer,
            data::Instrument::isPositionFileAggregated(aggregateTilesMode),
            laneInfo_,
            myCurrentTileInfo,
            clustersCounts_,
            ignoreMissingPositions_,
            outputBuffer,
            patternedFlowcellPositions_
            ), longTileId + " Task: Position load");
    }

    // If there's no cycle data (missing bcl tile file),
    // ensure the output buffer is empty and return true - i.e., don't submit the empty task.
    if (!inputBuffer.hasCycleData())
    {
        BCL2FASTQ_LOG(common::LogLevel::DEBUG) << " Empty tile: " << tileInfoIter->getNumber() << std::endl;
        for (data::BclBuffer& bclBuffer : outputBuffer)
        {
            for (auto& bufferForSingleCycle : bclBuffer.bcls_)
            {
                bufferForSingleCycle.bcls_.resize(0);
            }
        }

        return true;
    }

    BclFilesContainer::size_type bclIdx = 0;
    for (const auto& readInfo : laneInfo_.readInfos())
    {
        // We want all the cycles for the index reads, even if they were masked.
        for (const auto& cycleInfo : readInfo.cyclesToLoad())
        {
            taskQueue.addData(std::make_shared<BclLoadTask>(
                taskManager,
                inputBuffer,
                (aggregateTilesMode == common::TileAggregationMode::AGGREGATED ||
                    aggregateTilesMode == common::TileAggregationMode::CBCL) ? &aggregatedBclMutexes_[bclIdx] : NULL,
                aggregateTilesMode,
                laneInfo_,
                cycleInfo,
                tileInfoIter,
                ignoreMissingBcls_,
                bclIdx,
                *allOutputBuffers
            ));

            ++bclIdx;
        }
    }

    taskQueue.addData(std::make_shared<FilterLoadTask>(
        taskManager,
        inputBuffer,
        aggregateTilesMode,
        laneInfo_,
        myCurrentTileInfo,
        clustersCounts_,
        ignoreMissingFilters_,
        outputBuffer
    ), longTileId + " Task: filter-load");

    taskQueue.addData(std::make_shared<ControlLoadTask>(
        taskManager,
        inputBuffer,
        laneInfo_,
        ignoreMissingControls_,
        outputBuffer
    ));

    currentTileInfo_ += outputBuffer.size();

    return true;
}

} // namespace conversion
} // namespace bcl2fastq


