/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file TileBclFileReader.hh
 *
 * \brief Declaration of BCL reader for a single tile.
 *
 * \author Aaron Day
 */

#ifndef BCL2FASTQ_DATA_TILE_BCL_FILE_READER_HH
#define BCL2FASTQ_DATA_TILE_BCL_FILE_READER_HH

#include "data/BclFileReader.hh"
#include "data/RawBclBuffer.hh"
#include "io/SyncFile.hh"

#include <boost/filesystem/path.hpp>
#include <boost/noncopyable.hpp>

#include <memory>

namespace bcl2fastq {
namespace data {


class TileBclFileReader : public BclFileReaderT<io::FileReaderBase>
{
public:
    TileBclFileReader(const::boost::filesystem::path&     inputDir,
                      const layout::LaneInfo&             laneInfo,
                      common::CycleNumber                 cycleNumber,
                      size_t                              cycleIndex,
                      bool                                ignoreMissingBcls,
                      std::shared_ptr<io::FileReaderBase> bclFile);

    virtual bool read(RawBclBufferGroup& outputBuffer);
};

} // namespace data
} // namespace bcl2fastq

#endif // BCL2FASTQ_DATA_TILE_BCL_FILE_READER_HH

