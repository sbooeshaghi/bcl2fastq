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
 * \brief Declaration of BCL reader.
 *
 * \author Aaron Day
 */

#ifndef BCL2FASTQ_DATA_BCL_FILE_READER_HH
#define BCL2FASTQ_DATA_BCL_FILE_READER_HH

#include "data/RawBclBuffer.hh"
#include "io/SyncFile.hh"
#include "layout/LaneInfo.hh"
#include "common/Types.hh"

#include <boost/filesystem/path.hpp>
#include <boost/noncopyable.hpp>

#include <memory>

namespace bcl2fastq {
namespace data {

class BclFileReader : private boost::noncopyable
{
public:
    BclFileReader(const::boost::filesystem::path& inputDir,
                  const layout::LaneInfo&         laneInfo,
                  common::CycleNumber             cycleNumber,
                  size_t                          cycleIndex,
                  bool                            ignoreMissingBcls)
        : inputDir_(inputDir)
        , laneInfo_(laneInfo)
        , cycleNumber_(cycleNumber)
        , cycleIndex_(cycleIndex)
        , ignoreMissingBcls_(ignoreMissingBcls)
    {
    }

    virtual ~BclFileReader() { }
    virtual bool read(RawBclBufferGroup& outputBuffer) = 0;

protected:
    const boost::filesystem::path inputDir_;
    const layout::LaneInfo&       laneInfo_;
    common::CycleNumber           cycleNumber_;
    size_t                        cycleIndex_;
    bool                          ignoreMissingBcls_;
};

template<class BclFileType>
class BclFileReaderT : public BclFileReader
{
public:
    BclFileReaderT(const::boost::filesystem::path& inputDir,
                   const layout::LaneInfo&         laneInfo,
                   common::CycleNumber             cycleNumber,
                   size_t                          cycleIndex,
                   bool                            ignoreMissingBcls,
                   std::shared_ptr<BclFileType>    bclFile)
        : BclFileReader(inputDir,
                        laneInfo,
                        cycleNumber,
                        cycleIndex,
                        ignoreMissingBcls)
        , bclFile_(bclFile)
    {
    }

    virtual ~BclFileReaderT() { }

protected:
    std::shared_ptr<BclFileType> bclFile_;
};

} // namespace data
} // namespace bcl2fastq

#endif // BCL2FASTQ_DATA_BCL_FILE_READER_HH

