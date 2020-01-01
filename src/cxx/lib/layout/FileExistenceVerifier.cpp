/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file FileExistenceVerifier.cpp
 *
 * \brief Implementation of FileExistenceVerifier.
 *
 * \author Aaron Day
 */

#include "layout/FileExistenceVerifier.hh"

#include "layout/LaneInfo.hh"
#include "data/FilterFile.hh"
#include "data/PositionsFile.hh"
#include "data/BclFile.hh"
#include "data/CbclFile.hh"

namespace bcl2fastq
{
namespace layout
{

void FileExistenceVerifier::verifyAllFilesExist(const boost::filesystem::path& inputDir,
                                                const boost::filesystem::path& intensitiesDir,
                                                const LaneInfo&                laneInfo,
                                                const common::TileFileMap&     tileFileNameMap,
                                                common::TileAggregationMode    aggregateTilesMode,
                                                bool                           isPatternedFlowcell,
                                                bool                           ignoreMissingBcls,
                                                bool                           ignoreMissingFilters,
                                                bool                           ignoreMissingPositions)
{
    if (aggregateTilesMode == common::TileAggregationMode::AGGREGATED)
    {
        verifyFilesExist(inputDir,
                         intensitiesDir,
                         laneInfo,
                         tileFileNameMap,
                         aggregateTilesMode,
                         isPatternedFlowcell,
                         ignoreMissingBcls,
                         ignoreMissingFilters,
                         ignoreMissingPositions);
    }
    else
    {
        for (const auto& tileInfo : laneInfo.getTileInfos())
        {
            verifyFilesExist(inputDir,
                             intensitiesDir,
                             laneInfo,
                             tileFileNameMap,
                             aggregateTilesMode,
                             isPatternedFlowcell,
                             ignoreMissingBcls,
                             ignoreMissingFilters,
                             ignoreMissingPositions,
                             tileInfo.getNumber());
        }
    }
}

void FileExistenceVerifier::verifyFilesExist(const boost::filesystem::path& inputDir, 
                                             const boost::filesystem::path& intensitiesDir,
                                             const LaneInfo&                laneInfo,
                                             const common::TileFileMap&     tileFileNameMap,
                                             common::TileAggregationMode    aggregateTilesMode,
                                             bool                           isPatternedFlowcell,
                                             bool                           ignoreMissingBcls,
                                             bool                           ignoreMissingFilters,
                                             bool                           ignoreMissingPositions,
                                             const common::TileNumber       tileNumber /*= 0*/)
{
    for (const auto& readInfo : laneInfo.readInfos())
    {
        if (!ignoreMissingBcls)
        {
            for (const layout::CycleInfo &cycleInfo : readInfo.cycleInfos())
            {
                verifyBclFileExists(inputDir,
                                    aggregateTilesMode,
                                    laneInfo.getNumber(),
                                    cycleInfo.getNumber(),
                                    tileNumber,
                                    tileFileNameMap);
            }
        }

        if (!ignoreMissingFilters)
        {
            verifyFilterFileExists(inputDir,
                                   aggregateTilesMode,
                                   laneInfo.getNumber(),
                                   tileNumber,
                                   readInfo.getNumber());
        }
    }

    if (!ignoreMissingPositions)
    {
        verifyPositionsFileExists(intensitiesDir,
                                  aggregateTilesMode,
                                  isPatternedFlowcell,
                                  laneInfo.getNumber(),
                                  tileNumber);
    }
}

void FileExistenceVerifier::verifyBclFileExists(const boost::filesystem::path& inputDir,
                                                common::TileAggregationMode    aggregateTilesMode,
                                                common::LaneNumber             laneNumber,
                                                common::CycleNumber            cycleNumber,
                                                common::TileNumber             tileNumber,
                                                const common::TileFileMap&     tileFileNameMap)
{
    boost::filesystem::path fileName;
    switch (aggregateTilesMode)
    {
    case common::TileAggregationMode::AGGREGATED:
        data::BclFile::getAndVerifyFileName(inputDir,
                                            laneNumber,
                                            cycleNumber,
                                            false,
                                            fileName);
        break;
    case common::TileAggregationMode::NON_AGGREGATED:
        data::BclFile::getAndVerifyFileName(inputDir,
                                            laneNumber,
                                            tileNumber,
                                            cycleNumber,
                                            false,
                                            fileName);
        break;
    case common::TileAggregationMode::CBCL:
        fileName = data::CbclFileReader::getFileName(inputDir,
                                                     laneNumber,
                                                     cycleNumber,
                                                     tileNumber,
                                                     tileFileNameMap,
                                                     false);
        break;
    default:
        BCL2FASTQ_ASSERT_MSG(false, "Unexpected aggregation mode type");
    }
}

void FileExistenceVerifier::verifyFilterFileExists(const boost::filesystem::path& inputDir,
                                                   common::TileAggregationMode    aggregateTilesMode,
                                                   common::LaneNumber             laneNumber,
                                                   common::TileNumber             tileNumber,
                                                   common::ReadNumber             readNumber)
{
    size_t headerSize = 0;
    boost::filesystem::path filePath;
    if (!(data::FilterFile::doesFileExist(inputDir,
                                          aggregateTilesMode == common::TileAggregationMode::AGGREGATED,
                                          laneNumber,
                                          tileNumber,
                                          headerSize,
                                          filePath)))
    {
        throwException("filter",
                       aggregateTilesMode,
                       laneNumber,
                       tileNumber);
    }
}

void FileExistenceVerifier::verifyPositionsFileExists(const boost::filesystem::path& intensitiesDir,
                                                      common::TileAggregationMode    aggregateTilesMode,
                                                      bool                           isPatternedFlowcell,
                                                      common::LaneNumber             laneNumber,
                                                      common::TileNumber             tileNumber)
{
    size_t headerSize = 0;
    boost::filesystem::path posFilePath;
    if (!data::PositionsFileFactory::doesFileExist(intensitiesDir,
                                                   aggregateTilesMode == common::TileAggregationMode::AGGREGATED ||
                                                        aggregateTilesMode == common::TileAggregationMode::CBCL,
                                                   isPatternedFlowcell,
                                                   laneNumber,
                                                   tileNumber,
                                                   headerSize,
                                                   posFilePath))
    {
        throwException("positions",
                       aggregateTilesMode,
                       laneNumber,
                       tileNumber);
    }
}

void FileExistenceVerifier::throwException(const std::string&          fileType,
                                           common::TileAggregationMode aggregateTilesMode,
                                           common::LaneNumber          laneNumber,
                                           common::TileNumber          tileNumber)
{
    std::stringstream str;
    str << "Unable to find " << fileType << " file for lane: " << laneNumber;
    if (aggregateTilesMode == common::TileAggregationMode::NON_AGGREGATED)
    {
        str << " and tile: " << tileNumber;
    }
    str << std::endl;

    BOOST_THROW_EXCEPTION(common::IoError(ENOENT, str.str()));
}

} // namespace layout
} // namespace bcl2fastq

