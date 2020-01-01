/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file ConfigXml.cpp
 *
 * \brief Implementation of config.xml helper.
 *
 * \author Marek Balint
 */


#include <boost/optional.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

#include "common/Logger.hh"
#include "io/Xml.hh"
#include "layout/ConfigXml.hh"


namespace bcl2fastq {


namespace layout {


namespace detail {


/// \brief Generate path to config.xml file.
/// \param inputDir Path to input directory.
/// \return Path to config.xml file.
boost::filesystem::path generateConfigXmlFilePath(const boost::filesystem::path& inputDir)
{
    return boost::filesystem::path(inputDir / boost::filesystem::path("config.xml"));
}


} // namespace detail


ConfigXml::ConfigXml(boost::property_tree::ptree configData)
: ptree_(configData)
, tiles_()
{
    this->initTileMetadata();
}

ConfigXml::TileMetadataContainer::const_iterator ConfigXml::tileMetadataBegin(common::LaneNumber laneNumber) const
{
    BCL2FASTQ_ASSERT_MSG(laneNumber > 0, "Invalid lane number: " << laneNumber);
    BCL2FASTQ_ASSERT_MSG(laneNumber <= tiles_.size(), "Invalid lane number: " << laneNumber);
    return tiles_.at(laneNumber - 1).begin();
}

ConfigXml::TileMetadataContainer::const_iterator ConfigXml::tileMetadataEnd(common::LaneNumber laneNumber) const
{
    BCL2FASTQ_ASSERT_MSG(laneNumber > 0, "Invalid lane number: " << laneNumber);
    BCL2FASTQ_ASSERT_MSG(laneNumber <= tiles_.size(), "Invalid lane number: " << laneNumber);
    return tiles_.at(laneNumber - 1).end();
}

void ConfigXml::initTileMetadata()
{
    boost::optional<boost::property_tree::ptree &> tileSelectionElement = ptree_.get_child_optional("BaseCallAnalysis.Run.TileSelection");
    if (!tileSelectionElement) {
        BOOST_THROW_EXCEPTION(common::InputDataError("Required element 'BaseCallAnalysis.Run.TileSelection' missing in config.xml file"));
    }

    BOOST_FOREACH (boost::property_tree::ptree::value_type &laneElement, tileSelectionElement.get())
    {
        static const std::string laneElementName("Lane");
        if (laneElement.first != laneElementName)
        {
            continue;
        }

        boost::optional<common::LaneNumber> laneNumber = laneElement.second.get_optional<common::LaneNumber>("<xmlattr>.Index");
        if (!laneNumber)
        {
            BOOST_THROW_EXCEPTION(common::InputDataError("Required attrubute 'Index' missing in element element 'BaseCallAnalysis.Run.BaseCallParameters.TileSelection.Lane' in config.xml file"));
        }
        if (tiles_.size() <= laneNumber.get())
        {
            tiles_.resize(laneNumber.get());
        }
        BOOST_FOREACH (const boost::property_tree::ptree::value_type &tileElement, laneElement.second)
        {
            static const std::string tileElementName("Tile");
            if (tileElement.first != tileElementName)
            {
                continue;
            }

            TileMetadata tileMetadata;
            tileMetadata.tileNumber_ = tileElement.second.get_value<common::TileNumber>();
            tiles_.at(laneNumber.get() - 1).push_back(tileMetadata);
        }
    }
}


bool checkConfigXml(const boost::filesystem::path& inputDir)
{
    boost::filesystem::path configXmlFileName(detail::generateConfigXmlFilePath(inputDir));
    return boost::filesystem::ifstream(configXmlFileName).good() && !boost::filesystem::is_directory(configXmlFileName);
}

ConfigXml createConfigXml(const boost::filesystem::path& inputDir)
{
    boost::filesystem::path file = detail::generateConfigXmlFilePath(inputDir);
    BCL2FASTQ_LOG(common::LogLevel::INFO) << "config.xml: " << file << std::endl;
    return ConfigXml(io::parseXmlFile(file));
}


} // namespace layout


} // namespace bcl2fastq


