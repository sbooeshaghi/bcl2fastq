/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA^M
 * and certain third party copyright/licenses, and any user of this^M
 * source file is bound by the terms therein.^M
 *
 * \file AggregatedBclFileReader.cpp
 *
 * \brief Implementation of Bcl file reader for aggregated tiles.
 *
 * \author Aaron Day
 */

#include "data/AggregatedBclFileReader.hh"

namespace bcl2fastq {
namespace data {

AggregatedBclFileReader::AggregatedBclFileReader(const::boost::filesystem::path&     inputDir,
                                                 const layout::LaneInfo&             laneInfo,
                                                 common::CycleNumber                 cycleNumber,
                                                 size_t                              cycleIndex,
                                                 bool                                ignoreMissingBcls,
                                                 std::shared_ptr<io::SyncFile>       bclFile,
                                                 std::shared_ptr<data::CycleBCIFile> cycleBciFile)
: BclFileReaderT<io::SyncFile>(inputDir,
                               laneInfo,
                               cycleNumber,
                               cycleIndex,
                               ignoreMissingBcls,
                               bclFile)
, cycleBciFile_(cycleBciFile)
{
}

bool AggregatedBclFileReader::read(RawBclBufferGroup& outputBuffer)
{
    if (outputBuffer.size() == 0) return false;

    if (!openAggFileIfNeeded())
    {
        outputBuffer.clearClusters(cycleIndex_, bclFile_);
        return false;
    }

    bool success = true;
    for (auto& outputBufferForTile : outputBuffer)
    {
        auto tileInfoIter = outputBufferForTile.tileInfo_;

        io::SyncFile::SyncFileReader bclFileReader(*bclFile_,
                                                   std::distance(laneInfo_.getTileInfos().begin(),
                                                                 tileInfoIter));

        outputBufferForTile.cycleData_[cycleIndex_].bcls_.path_ = bclFile_->getPath();

        // 1st gzip block with data for tile
        data::CycleBCIFile::Record firstBlockBciRecord = cycleBciFile_->getRecord(tileInfoIter->getIndex());
        bclFileReader.seek(firstBlockBciRecord.compressedOffset, std::ios_base::beg);

        outputBufferForTile.uncompressedBclOffset_ = firstBlockBciRecord.uncompressedOffset;

        // Look for the next tile that is inside of the last gzip block used by the current tile.
        auto nextTileInsideBlock = tileInfoIter + 1;

        // Determine where to stop reading.
        //
        // nextStartPos      is the start of a gzip-block containing the next tile.
        //                   Unless the uncompressed offset is 0, the current tile will
        //                   end in this block so we need to read this too.
        //                   A value of zero is interpretted as EOF.
        // --------------------------------
        int64_t nextStartPos = 0;

        const bool isLastTileInBciFile = (tileInfoIter->getIndex()+1 == cycleBciFile_->getTilesCount());
        if (!isLastTileInBciFile) // && nextTileInsideBlock != laneInfo_.getTileInfos().end())
        {
            // Read enough gzp blocks to load tile.
            if (nextTileInsideBlock != laneInfo_.getTileInfos().end()) {
                data::CycleBCIFile::Record endRecord = cycleBciFile_->getRecord(nextTileInsideBlock->getIndex());

                // This is start of the next tile. Some clusters will be in the next gzip block. readGzipBlocks will read this too.
                nextStartPos = endRecord.compressedOffset;
            }
        }

        int64_t startPos = firstBlockBciRecord.compressedOffset;
//BCL2FASTQ_LOG(common::LogLevel::NONE) << "Tile: " << tileInfoIter->getNumber() << " startPos: " << startPos << " next: " << nextStartPos << std::endl;

        if (!bclFileReader.readGzipBlocks(outputBufferForTile.cycleData_[cycleIndex_].bcls_,
                                          startPos, nextStartPos))
        {
            success = false;
            break;
        }
    }

    if (!success) {

        std::stringstream msg;
        msg << "Found corrupt aggregated BCL data in cycle " << (1 + cycleIndex_)
            << " starting at tile " << outputBuffer.begin()->tileInfo_->getNumber() << ". ";

        if (ignoreMissingBcls_) {

            msg << "Data for this cycle will be will be ignored. ";
            msg << "File: '" << bclFile_->getPath() << "'" << std::endl;
            BCL2FASTQ_LOG(common::LogLevel::WARNING) << msg.str();

            // This will clear cluster bcl_ data for the current cycle and tile.
            outputBuffer.clearClusters(cycleIndex_, bclFile_);

            // We also want to ignore ALL remaining tiles for this cycle.
            bclFile_.reset();
        }
        else
        {
            msg << "Aborting. ";
            msg << "File: '" << bclFile_->getPath() << "'" << std::endl;
            BCL2FASTQ_LOG(common::LogLevel::WARNING) << msg.str();

            BOOST_THROW_EXCEPTION(common::InputDataError(msg.str()));
        }
    }

    return success;
}

bool AggregatedBclFileReader::openAggFileIfNeeded()
{
    if (!bclFile_)
    {
        return false;
    }

    if (!bclFile_->isOpen())
    {
        if (!cycleBciFile_->openFile(inputDir_, laneInfo_.getNumber(), cycleNumber_))
        {
            return false;
        }

        boost::filesystem::path bclPath;
        if (!data::BclFile::getAndVerifyFileName(inputDir_,
                                            laneInfo_.getNumber(),
                                            cycleNumber_,
                                            ignoreMissingBcls_,
                                            bclPath))
        {
            return false;
        }

        bclFile_->openFile(bclPath,
                           ignoreMissingBcls_,
                           0); // There is only 1 file. It is definitely index 0.
    }

    return true;
}

} // namespace data
} // namespace bcl2fastq

