/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file TileInfo.hh
 *
 * \brief Declaration of tile metadata.
 *
 * \author Marek Balint
 */


#ifndef BCL2FASTQ_LAYOUT_TILEINFO_HH
#define BCL2FASTQ_LAYOUT_TILEINFO_HH


#include "common/Types.hh"


namespace bcl2fastq {


namespace layout {


/// \brief Tile metadata.
class TileInfo
{
public:

    /// \brief Default Constructor.
    TileInfo();

    /// \brief Constructor.
    /// \param number Tile number.
    /// \param index Tile index.
    /// \param skippedTiles Number of skipped tiles before this tile.
    explicit TileInfo(
        common::TileNumber number,
        size_t index,
        size_t skippedTiles
    );

    /// \brief Constructor.
    /// \param number Tile number.
    /// \param index Tile index.
    /// \param clustersCount Number of clusters.
    /// \param skippedTiles Number of skipped tiles before this tile.
    /// \param skippedClusters Number of skipped clusters before this tile.
    TileInfo(
        common::TileNumber number,
        size_t index,
        common::ClustersCount clustersCount,
        size_t skippedTiles,
        size_t skippedClusters
    );

public:

    /// \brief Get tile number.
    /// \return Tile number.
    common::TileNumber getNumber() const { return number_; }

    /// \brief Get tile index.
    /// \return Tile index.
    common::TileNumber getIndex() const { return index_; }

    /// \brief Get clusters count.
    /// \return Number of clusters on this tile.
    /// \pre <tt>this->haveClustersCount() == true</tt>
    common::TotalClustersCount getClustersCount() const { return static_cast<common::TotalClustersCount>(clustersCount_); }

    /// \brief Check whether clusters count is known.
    /// \retval true Clusters count is known.
    /// \retval false Clusters count is not known.
    bool haveClustersCount() const { return haveClustersCount_; }

    /// \brief Get number of skipped tiles.
    /// \return Number of tiles skipped before this tile.
    size_t getSkippedTilesCount() const { return skippedTiles_; }

    /// \brief Get number of skipped clusters.
    /// \return Number of clusters skipped before this tile.
    size_t getSkippedClustersCount() const { return skippedClusters_; }

    /// \brief Overloaded comparison operator.
    bool operator<( const TileInfo &tile ) const  { return index_ < tile.index_; }
private:

    /// \brief Tile number.
    ///
    /// Source: BCI file (Nova), config.xml (HiSeq/MiSeq)
    common::TileNumber number_;

    /// \brief Tile index.
    ///
    /// Source: calculated (index in container of tiles)
    size_t index_;

    /// \brief Number of clusters on the tile.
    ///
    /// Source: BCI file (Nova), config.xml (HiSeq/MiSeq)
    common::ClustersCount clustersCount_;

    /// \brief Clusters count known flag.
    ///
    /// Source: true only if BCI file is present
    bool haveClustersCount_;

    /// \brief Number of tiles skipped before this tile.
    ///
    /// Source: command line argument --tiles
    size_t skippedTiles_;

    /// \brief Number of clusters skipped before this tile.
    ///
    /// Source: command line argument --tiles
    size_t skippedClusters_;
};


} // namespace layout


} // namespace bcl2fastq


#endif // BCL2FASTQ_LAYOUT_TILEINFO_HH


