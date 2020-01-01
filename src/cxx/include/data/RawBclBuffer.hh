/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file BclReader.hh
 *
 * \brief Declaration of raw BCL buffer.
 *
 * \author Aaron Day
 */

#ifndef BCL2FASTQ_DATA_RAW_BCL_BUFFER_HH
#define BCL2FASTQ_DATA_RAW_BCL_BUFFER_HH


#include "common/Types.hh"
#include "layout/LaneInfo.hh"
#include "io/GzipDecompressor.hh"
#include "data/BclFile.hh"

#include <utility>
#include <memory>
#include <vector>
#include <atomic>
#include <mutex>
#include <condition_variable>

namespace bcl2fastq {
namespace data {

struct PerCycleData
{
    /// \brief Return the number of clusters per byte
    common::ClustersCount getNumClustersPerByte() const { return (8 / (numBitsPerQscore_ + 2)); }

    /// \brief Return the number of clusters
    common::ClustersCount getNumClusters() const { return bcls_.size() * getNumClustersPerByte(); }

    /// \brief BCL cycles data container type definition.
    typedef std::vector<char> BclsContainer;

    PerCycleData() : bcls_(), includeNonPf_(true), numBitsPerQscore_(6), remappedQscores_() { }

    /// \brief Buffers for cycles data.
    common::RawDataBuffer bcls_;

    /// \brief Indicates which cycles include non-pf clusters.
    bool includeNonPf_;

    /// \brief Number of bits per qscore for each cycle.
    uint8_t numBitsPerQscore_;

    /// \brief Remapped qscores for each cycle.
    std::vector<uint32_t> remappedQscores_;

    /// \brief Uncompressed block size.
    uint32_t uncompressedBlockSize_;
};

struct RawBclBuffer
{
    RawBclBuffer()
     : cycleData_()
     , gzipDecompressors_()
     , controls_()
     , filters_()
     , positions_()
     , tileOffsets_()
     , tileInfo_()
     , groupNumber_(0)
     , uncompressedBclOffset_(0)
    {
    }

    uint32_t getGroupNumber() const { return groupNumber_; }
    void setGroupNumber(uint32_t groupNumber) { groupNumber_ = groupNumber; }

    bool hasCycleData() const {
        for (PerCycleData perCycleData : cycleData_) {
            if (!(perCycleData.bcls_.empty())) return true;
        }
        return false;
    }

    // One for each cycle
    std::vector<PerCycleData> cycleData_;
    std::vector<io::GzipDecompressor> gzipDecompressors_;

    common::RawDataBuffer controls_;
    common::RawDataBuffer filters_;

    // Could be shared by all tiles if it is a patterned flowcell
    std::shared_ptr<common::RawDataBuffer> positions_;

    std::vector<std::pair<size_t, common::TileNumber>> tileOffsets_;

    /// \brief Tile meta data.
    layout::LaneInfo::TileInfosContainer::const_iterator tileInfo_;

    uint32_t groupNumber_;

    uint32_t uncompressedBclOffset_;
};

typedef std::vector<RawBclBuffer> RawBclBufferVec;

class RawBclBufferGroup : private boost::noncopyable
{
public:
    RawBclBufferGroup(size_t numTilesPerBuffer, size_t numCycles) : outputBuffer_(numTilesPerBuffer), groupNumber_(0), prevBuffer_(), prevBclFiles_(numCycles), isFinished_(false), mut_(), cv_() { }

    uint32_t getGroupNumber() const { return groupNumber_; }
    void setGroupNumber(uint32_t groupNumber) { isFinished_ = false; groupNumber_ = groupNumber; }

    size_t size() const { return outputBuffer_.size(); }

    RawBclBufferVec::iterator begin() { return outputBuffer_.begin(); }
    RawBclBufferVec::iterator end() { return outputBuffer_.end(); }

    RawBclBufferVec::const_iterator begin() const { return outputBuffer_.begin(); }
    RawBclBufferVec::const_iterator end() const { return outputBuffer_.end(); }

    bool hasCycleData() {
        for (RawBclBuffer &rawBclBuffer : outputBuffer_) {
            if (rawBclBuffer.hasCycleData()) return true;
        }
        return false;
    }

    bool hasCycleData(common::CycleNumber cycleIndex) {
        for (RawBclBuffer &rawBclBuffer : outputBuffer_) {
            PerCycleData &perCycleData = rawBclBuffer.cycleData_[cycleIndex];
            if (!(perCycleData.bcls_.empty())) return true;
        }
        return false;
    }

    // This will clear BCL cluster data for the specified cycle and current tile.
    template<class FileReaderType>
    void clearClusters(common::CycleNumber cycleIndex,
                       std::shared_ptr<FileReaderType> bclFile)
    {
        BCL2FASTQ_LOG(common::LogLevel::TRACE)
            << "Clearing BCL data for cycle " << (1+cycleIndex) << std::endl;

        for (data::RawBclBuffer& rawBclBuffer : outputBuffer_)
        {
            rawBclBuffer.cycleData_[cycleIndex].uncompressedBlockSize_ = 0;
            rawBclBuffer.cycleData_[cycleIndex].bcls_.resize(0);
            if (bclFile)
            {
                rawBclBuffer.cycleData_[cycleIndex].bcls_.path_ = bclFile->getPath();
            }
        }
    }

    void resize(size_t newSize) { outputBuffer_.resize(newSize); }

    void reset()
    {
        for (auto& buffer : outputBuffer_)
        {
            for (auto& decompressor : buffer.gzipDecompressors_)
            {
                decompressor.clear();
            }
        }

        prevBuffer_ = std::shared_ptr<RawBclBufferGroup>();

        for (auto& bclFile : prevBclFiles_)
        {
            bclFile = std::shared_ptr<data::BclFile>();
        }
    }

    void setPrevBuffer(std::shared_ptr<RawBclBufferGroup>& prevBuffer) { prevBuffer_ = prevBuffer; }
    std::shared_ptr<RawBclBufferGroup> getPrevBuffer() { return prevBuffer_; }

    std::shared_ptr<data::BclFile> getPrevBclFileForCycle(common::CycleNumber cycleNumber) { return prevBclFiles_[cycleNumber]; }
    void setPrevBclFileForCycle(common::CycleNumber cycleNumber, std::shared_ptr<data::BclFile>& bclFile) { prevBclFiles_[cycleNumber] = bclFile; }

    void setFinished() { isFinished_ = true; cv_.notify_one(); }
    void waitForFinished() { std::unique_lock<std::mutex> lock(mut_); cv_.wait(lock, [this] { return isFinished_ == true; }); }

private:
    RawBclBufferVec outputBuffer_;
    uint32_t groupNumber_;

    std::shared_ptr<RawBclBufferGroup> prevBuffer_;
    std::vector<std::shared_ptr<data::BclFile>> prevBclFiles_;

    std::atomic<bool> isFinished_;
    std::mutex mut_;
    std::condition_variable cv_;
};

typedef std::shared_ptr<RawBclBufferGroup> RawBclBufferGroupPtr;

} // namespace conversion
} // namespace bcl2fastq

#endif // BCL2FASTQ_DATA_RAW_BCL_BUFFER_HH

