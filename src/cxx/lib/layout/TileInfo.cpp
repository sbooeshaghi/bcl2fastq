/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file TileInfo.cpp
 *
 * \brief Implementation of tile metadata.
 *
 * \author Marek Balint
 */


#include "layout/TileInfo.hh"

#include "common/Logger.hh"

#include <boost/format.hpp>

namespace bcl2fastq
{
namespace layout
{


TileInfo::TileInfo()
: number_()
, index_(0)
, clustersCount_()
, haveClustersCount_(false)
, skippedTiles_(0)
, skippedClusters_(0)
{}

TileInfo::TileInfo(
    common::TileNumber number,
    size_t index,
    size_t skippedTiles
)
: number_(number)
, index_(index)
, clustersCount_()
, haveClustersCount_(false)
, skippedTiles_(skippedTiles)
, skippedClusters_(0)
{
    BCL2FASTQ_LOG(common::LogLevel::INFO) << "  Tile: " << this->getNumber()
                                          << " (index: " << this->getIndex()
                                          << (this->haveClustersCount() ? ( boost::format(", clusters count: %d")
                                                                                 % this->getClustersCount() ).str()
                                                                        : std::string(""))
                                          << ", skipped tiles: " << this->getSkippedTilesCount() << ")" << std::endl;
}

TileInfo::TileInfo(
    common::TileNumber number,
    size_t index,
    common::ClustersCount clustersCount,
    size_t skippedTiles,
    size_t skippedClusters
)
: number_(number)
, index_(index)
, clustersCount_(clustersCount)
, haveClustersCount_(true)
, skippedTiles_(skippedTiles)
, skippedClusters_(skippedClusters)
{
    BCL2FASTQ_LOG(common::LogLevel::INFO) << "  Tile: " << this->getNumber()
                                          << " (index: " << this->getIndex()
                                          << (this->haveClustersCount() ? ( boost::format(", clusters count: %d")
                                                                                 % this->getClustersCount() ).str()
                                                                        : std::string(""))
                                          << ", skipped tiles: " << this->getSkippedTilesCount()
                                          << ", skipped clusters: " << this->getSkippedClustersCount() << ")" << std::endl;
}


} // namespace layout
} // namespace bcl2fastq


