/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file RunInfoXml.hh
 *
 * \brief Declaration of RunInfo.xml helper.
 *
 * \author Marek Balint
 * \author Mauricio Varea
 */


#ifndef BCL2FASTQ_LAYOUT_RUNINFOXML_HH
#define BCL2FASTQ_LAYOUT_RUNINFOXML_HH


#include <boost/filesystem/path.hpp>
#include <boost/property_tree/ptree.hpp>

#include "common/Types.hh"
#include "layout/LaneInfo.hh"
#include "layout/ReadMetadata.hh"

namespace bcl2fastq {


namespace layout {

/// \brief RunInfo.xml file helper.
class RunInfoXml
{
public:

    /// \brief Tile number data container type definition.
    typedef std::vector<common::TileNumber> TileNumbersContainer;

public:

    /// \brief Constructor.
    /// \param runInfoData RunInfo XML data.
    explicit RunInfoXml(boost::property_tree::ptree runInfoData);

public:

    /// \brief Get instrument.
    /// \return Instrument.
    std::string getInstrument() const;

    /// \brief Get run id.
    /// \return Run Id.
    std::string getRunId() const;

    /// \brief Get run number.
    /// \return Run number.
    std::string getRunNumber() const;

    /// \brief Get flowcell ID.
    /// \return Flowcell ID.
    std::string getFlowcellId() const;

    /// \brief Get lanes count.
    /// \return Lanes count
    common::LaneNumber getLanesCount() const;

    /// \brief Get beginning of read meta data.
    /// \return Iterator to beginning of read meta data.
    ReadMetadataContainer::const_iterator readMetadataBegin() const;

    /// \brief Get end of read meta data.
    /// \return Iterator to end of read meta data.
    ReadMetadataContainer::const_iterator readMetadataEnd() const;

    /// \brief Get beginning of tile numbers data.
    /// \return Iterator to beginning of tile numbers data.
    RunInfoXml::TileNumbersContainer::const_iterator tileNumbersBegin( common::LaneNumber lane ) const;

    /// \brief Get end of tile numbers data.
    /// \return Iterator to end of tile numbers data.
    RunInfoXml::TileNumbersContainer::const_iterator tileNumbersEnd( common::LaneNumber lane ) const;
private:

    /// \brief Initialize read meta data container.
    void initReadMetadata();

    /// \brief Initialize tile numbers data container.
    void initTileNumbers();

private:

    /// \brief Property tree holding the actual data.
    boost::property_tree::ptree ptree_;

    /// \brief Read meta data container.
    ReadMetadataContainer reads_;

    /// \brief Tile numbers, grouped by lanes
    std::vector< TileNumbersContainer > tiles_;
};


/// \brief Create RunInfoXml object.
/// \param runfolderDir Path to runfolder.
/// \return RunInfoXml object constructed with data from given runfolder.
RunInfoXml createRunInfoXml(boost::filesystem::path runfolderDir);


} // namespace layout


} // namespace bcl2fastq


#endif // BCL2FASTQ_LAYOUT_RUNINFOXML_HH


