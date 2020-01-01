/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file ConfigXml.hh
 *
 * \brief Declaration of config.xml helper.
 *
 * \author Marek Balint
 */


#ifndef BCL2FASTQ_LAYOUT_CONFIG_HH
#define BCL2FASTQ_LAYOUT_CONFIG_HH


#include <boost/filesystem/path.hpp>
#include <boost/property_tree/ptree.hpp>

#include "common/Types.hh"
#include "layout/LaneInfo.hh"
#include "config/SampleSheetCsv.hh"


namespace bcl2fastq {


namespace layout {


/// \brief Confing.xml file helper.
class ConfigXml
{
public:

    /// \brief Tile meta data.
    struct TileMetadata
    {
    public:

        /// \brief Tile number.
        common::TileNumber tileNumber_;
    };

    /// \brief Tile meta data container type definition.
    typedef std::vector<TileMetadata> TileMetadataContainer;

public:

    /// \brief Constructor.
    /// \param configData Config XML data.
    explicit ConfigXml(boost::property_tree::ptree configData);

public:

    /// \brief Check if config.xml contains the lane number
    /// \param laneNumber Lane number
    /// \return True if config.xml contains the lane number.
    bool containsLaneNumber(common::LaneNumber laneNumber) const { return (laneNumber > 0 && laneNumber <= tiles_.size()); }

    /// \brief Get beginning of tile meta data.
    /// \param laneNumber Lane number.
    /// \return Iterator to beginning of tile meta data.
    ConfigXml::TileMetadataContainer::const_iterator tileMetadataBegin(common::LaneNumber laneNumber) const;

    /// \brief Get end of tile meta data.
    /// \param laneNumber Lane number.
    /// \return Iterator to end of tile meta data.
    ConfigXml::TileMetadataContainer::const_iterator tileMetadataEnd(common::LaneNumber laneNumber) const;

private:

    /// \brief Initialize tile meta data container.
    void initTileMetadata();

private:

    /// \brief Property tree holding the actual data.
    boost::property_tree::ptree ptree_;

    /// \brief Tile meta data container.
    std::vector<TileMetadataContainer> tiles_;
};


/// \brief Check presence of config.xml.
/// \param inputDir Path to input directory.
/// \retval true The config.xml file is present.
/// \retval false There is no config.xml.
bool checkConfigXml(const boost::filesystem::path& inputDir);


/// \brief Create ConfigXml object.
/// \param inputDir Path to input directory.
/// \return ConfigXml object constructed with data from given input directory.
ConfigXml createConfigXml(const boost::filesystem::path& inputDir);


} // namespace layout


} // namespace bcl2fastq


#endif // BCL2FASTQ_LAYOUT_CONFIG_HH


