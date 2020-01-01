/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file BCIndex.hh
 *
 * \brief Declaration of BCI helper.
 *
 * \author Marek Balint
 */


#ifndef BCL2FASTQ_LAYOUT_BCINDEX_HH
#define BCL2FASTQ_LAYOUT_BCINDEX_HH


#include <boost/filesystem/path.hpp>

#include "common/Types.hh"


namespace bcl2fastq {


namespace layout {


/// \brief BCI file helper.
class BCIndex
{
public:

    /// \brief Raw data type definition.
    typedef std::vector<char> RawData;

    /// \brief Tile meta data.
    struct TileMetadata
    {
    public:

        /// \brief Tile number.
        common::TileNumber tileNumber_;

        /// \brief Clusters count.
        common::ClustersCount clustersCount_;
    };

    /// \brief Tile meta data container type definition.
    typedef std::vector<TileMetadata> TileMetadataContainer;

public:

    /// \brief Constructor.
    /// \param bciData BCI data.
    explicit BCIndex(const RawData &bciData);

public:

    /// \brief Get beginning of tile meta data.
    /// \return Iterator to beginning of tile meta data.
    BCIndex::TileMetadataContainer::const_iterator tileMetadataBegin() const;

    /// \brief Get end of tile meta data.
    /// \return Iterator to end of tile meta data.
    BCIndex::TileMetadataContainer::const_iterator tileMetadataEnd() const;

private:

    /// \brief Initialize tile meta data container.
    /// \param bciData Source BCI data.
    void initTileMetadata(const RawData &bciData);

private:

    /// \brief Tile meta data container.
    TileMetadataContainer tiles_;
};


/// \brief Parse BCI file.
/// \param bciFile Path to BCI file to be parsed.
/// \return Data strucutre containing data parsed out of the BCI file.
BCIndex::RawData parseBciFile(const boost::filesystem::path& bciFile);

/// \brief Check presence of BCI file.
/// \param inputDir Path to input directory.
/// \param laneNumber Lane number.
/// \retval true The BCI file is present.
/// \retval false There is no BCI file.
bool checkBciFile(boost::filesystem::path inputDir, common::LaneNumber laneNumber);

/// \brief Create BCIndex object.
/// \param inputDir Path to input directory.
/// \param laneNumber Lane number.
/// \return BCIndex object constructed with data from given input directory.
BCIndex createBCIndex(boost::filesystem::path inputDir, common::LaneNumber laneNumber);


} // namespace layout


} // namespace bcl2fastq


#endif // BCL2FASTQ_LAYOUT_BCINDEX_HH


