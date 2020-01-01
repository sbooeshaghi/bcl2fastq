/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file ConversionStatsXml.cpp
 *
 * \brief Xml Serialization of Conversion statistics.
 *
 * \author Mauricio Varea
 */

#include <boost/lexical_cast.hpp>

#include "stats/ConversionStatsXml.hh"

namespace bcl2fastq
{
namespace stats
{

ConversionStatsXml::ConversionStatsXml()
{
}

boost::property_tree::path ConversionStatsXml::getPrefix(
    const std::string &flowcellId,
    const std::string &projectName,
    const std::string &sampleName,
    const std::string &barcodeName,
    common::TileNumber tileNumber,
    const unsigned lane )
{
    return boost::property_tree::path( "Stats"
                                       "/<indexed>Flowcell/<flowcell-id>" + flowcellId
                                      +"/<indexed>Project/<name>" + projectName
                                      +"/<indexed>Sample/<name>" + sampleName
                                      +"/<indexed>Barcode/<name>" + barcodeName
                                      +"/<indexed>Lane/<number>" + boost::lexical_cast<std::string>(lane)
                                      +"/<indexed>Tile/<number>" + boost::lexical_cast<std::string>(tileNumber)
                                      , '/');
}

void ConversionStatsXml::addLaneTile(
    const std::string &flowcellId,
    const std::string &projectName,
    const std::string &sampleName,
    const std::string &barcodeName,
    common::TileNumber tileNumber,
    const unsigned lane,
    const TileBarcodeStats& tileStats)
{
    for(TileBarcodeStats::size_type i=0; i < tileStats.size(); i++)
    {
        add( getPrefix( flowcellId, projectName, sampleName, barcodeName, tileNumber, lane )
           / TileBarcodeStats::str_type[i] / "ClusterCount", tileStats[i].clustersCount );
    }
}

void ConversionStatsXml::addLaneTileRead(
    const std::string &flowcellId,
    const std::string &projectName,
    const std::string &sampleName,
    const std::string &barcodeName,
    common::TileNumber tileNumber,
    common::ReadNumber readNumber,
    const unsigned lane,
    const ReadBarcodeStats& tileStats)
{
    const boost::property_tree::path readNode( "<indexed>Read/<number>"
                                             + boost::lexical_cast<std::string>(readNumber), '/');

    for(TileBarcodeStats::size_type i=0; i < tileStats.size(); i++)
    {
        add( getPrefix( flowcellId, projectName, sampleName, barcodeName, tileNumber, lane )
           / TileBarcodeStats::str_type[i] / readNode / "Yield", tileStats[i].yield );
        add( getPrefix( flowcellId, projectName, sampleName, barcodeName, tileNumber, lane )
           / TileBarcodeStats::str_type[i] / readNode / "YieldQ30", tileStats[i].yieldQ30 );
        add( getPrefix( flowcellId, projectName, sampleName, barcodeName, tileNumber, lane )
           / TileBarcodeStats::str_type[i] / readNode / "QualityScoreSum", tileStats[i].qualityScoreSum );
    }
}


void ConversionStatsXml::addFlowcellLane(
    const std::string &flowcellId,
    const unsigned lane,
    const stats::BarcodeHits::Popular &topUnknownBarcodes)
{
    const boost::property_tree::path laneValuePrefix("Stats"
                                      "/<indexed>Flowcell/<flowcell-id>" + flowcellId
                                      +"/<indexed>Lane/<number>" + boost::lexical_cast<std::string>(lane)
                                      +"/TopUnknownBarcodes", '/');

    BOOST_FOREACH(const stats::BarcodeHits::Popular::value_type &unknownBarcode, topUnknownBarcodes)
    {
        add( laneValuePrefix / ("<indexed>Barcode/<sequence>" + unknownBarcode.first.toString()).c_str()
                             / "<xmlattr>"
                             / "count",
             unknownBarcode.second );
    }
}


} //namespace stats
} //namespace bcl2fastq

