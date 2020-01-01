/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file FileExistenceVerifier.hh
 *
 * \brief Declaration of FileExistenceVerifier.
 *
 * \author Aaron Day
 */

#ifndef BCL2FASTQ_LAYOUT_FILE_EXISTENCE_VERIFIER_HH
#define BCL2FASTQ_LAYOUT_FILE_EXISTENCE_VERIFIER_HH

#include "common/Types.hh"

#include <boost/noncopyable.hpp>
#include <boost/filesystem/path.hpp>

namespace bcl2fastq
{
namespace layout
{

class LaneInfo;

class FileExistenceVerifier : private boost::noncopyable
{
public:

    static void verifyAllFilesExist(const boost::filesystem::path& inputDir,
                                    const boost::filesystem::path& intensitiesDir,
                                    const LaneInfo&                laneInfo,
                                    const common::TileFileMap&     tileFileNameMap,
                                    common::TileAggregationMode    aggregateTilesMode,
                                    bool                           isPatternedFlowcell,
                                    bool                           ignoreMissingBcls,
                                    bool                           ignoreMissingFilters,
                                    bool                           ignoreMissingPositions);

private:
    static void verifyFilesExist(const boost::filesystem::path& inputDir,
                                 const boost::filesystem::path& intensitiesDir,
                                 const LaneInfo&                laneInfo,
                                 const common::TileFileMap&     tileFileNameMap,
                                 common::TileAggregationMode    aggregateTilesMode,
                                 bool                           isPatternedFlowcell,
                                 bool                           ignoreMissingBcls,
                                 bool                           ignoreMissingFilters,
                                 bool                           ignoreMissingPositions,
                                 const common::TileNumber       tileNumber = 0);

    static void verifyFilterFileExists(const boost::filesystem::path& inputDir,
                                       common::TileAggregationMode    aggregateTilesMode,
                                       common::LaneNumber             laneNumber,
                                       common::TileNumber             tileNumber,
                                       common::ReadNumber             readNumber);

    static void verifyBclFileExists(const boost::filesystem::path& inputDir,
                                    common::TileAggregationMode    aggregateTilesMode,
                                    common::LaneNumber             laneNumber,
                                    common::CycleNumber            cycleNumber,
                                    common::TileNumber             tileNumber,
                                    const common::TileFileMap&     tileFileNameMap);

    static void verifyPositionsFileExists(const boost::filesystem::path& intensitiesDir,
                                          common::TileAggregationMode    aggregateTilesMode,
                                          bool                           isPatternedFlowcell,
                                          common::LaneNumber             laneNumber,
                                          common::TileNumber             tileNumber);

    static void throwException(const std::string&          fileType,
                               common::TileAggregationMode aggregateTilesMode,
                               common::LaneNumber          laneNumber,
                               common::TileNumber          tileNumber);
};

} // namespace layout 
} // namespace bcl2fastq

#endif // BCL2FASTQ_LAYOUT_FILE_EXISTENCE_VERIFIER_HH
