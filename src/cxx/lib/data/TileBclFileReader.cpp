/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA^M
 * and certain third party copyright/licenses, and any user of this^M
 * source file is bound by the terms therein.^M
 *
 * \file TileBclFileReader.cpp
 *
 * \brief Implementation of Bcl file reader for a single tile.
 *
 * \author Aaron Day
 */

#include "data/TileBclFileReader.hh"

namespace bcl2fastq {
namespace data {

TileBclFileReader::TileBclFileReader(const::boost::filesystem::path&     inputDir,
                                     const layout::LaneInfo&             laneInfo,
                                     common::CycleNumber                 cycleNumber,
                                     size_t                              cycleIndex,
                                     bool                                ignoreMissingBcls,
                                     std::shared_ptr<io::FileReaderBase> bclFile)
: BclFileReaderT<io::FileReaderBase>(inputDir,
                                     laneInfo,
                                     cycleNumber,
                                     cycleIndex,
                                     ignoreMissingBcls,
                                     bclFile)
{
}

bool TileBclFileReader::read(RawBclBufferGroup& outputBuffer)
{
    static int sequentialId = 0; ++sequentialId;
    for (auto& outputBufferForTile : outputBuffer)
    {
        common::RawDataBuffer &outputBuffer = outputBufferForTile.cycleData_.at(cycleIndex_).bcls_;

        boost::filesystem::path bclPath;
        if (!data::BclFile::getAndVerifyFileName(inputDir_,
                                                 laneInfo_.getNumber(),
                                                 outputBufferForTile.tileInfo_->getNumber(),
                                                 cycleNumber_,
                                                 ignoreMissingBcls_,
                                                 bclPath))
        {
            if (ignoreMissingBcls_) {
                BCL2FASTQ_LOG(common::LogLevel::WARNING) << "Unable to read tile file (" << 
                      outputBufferForTile.tileInfo_->getNumber()  << "): " << bclPath.string() << std::endl;
                outputBuffer.resize(0);
                continue;
            }
            else 
            {
                BOOST_THROW_EXCEPTION(common::IoError(ENOENT, (boost::format("Unable to open tile file '%s' for reading") % bclPath.string()).str()));
            }
        }

        BCL2FASTQ_LOG(common::LogLevel::DEBUG) << "TileBclFileReader::read " << sequentialId << " " << bclPath.string() << std::endl;

        io::UnprocessedFile bclFile(bclPath,
                                    ignoreMissingBcls_);

        bool status = bclFile.readEntireFile(outputBuffer);
        BCL2FASTQ_LOG(common::LogLevel::DEBUG) << "... done TileBclFileReader::read " << sequentialId << std::endl;
        return status;
    }

    return true;
}

} // namespace data
} // namespace bcl2fastq

