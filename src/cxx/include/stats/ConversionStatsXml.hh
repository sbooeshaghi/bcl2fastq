/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file ConversionStatsXml.hh
 *
 * \brief Xml Serialization of Conversion statistics.
 *
 * \author Mauricio Varea
 */

#ifndef BCL2FASTQ_STATS_CONVERSION_STATS_XML_HH
#define BCL2FASTQ_STATS_CONVERSION_STATS_XML_HH

#include "common/Types.hh"
#include "io/Xml.hh"
#include "layout/Layout.hh"

#include "stats/TileStats.hpp"
#include "stats/BarcodeHits.hh"

namespace bcl2fastq
{
namespace stats
{

class ConversionStatsXml : public boost::property_tree::ptree
{
public:
    ConversionStatsXml();

    void addLaneTile(
        const std::string &flowcellId,
        const std::string &projectName,
        const std::string &sampleName,
        const std::string &barcodeName,
        common::TileNumber tileNumber,
        const unsigned lane,
        const stats::TileBarcodeStats& tileStats);

    void addLaneTileRead(
        const std::string &flowcellId,
        const std::string &projectName,
        const std::string &sampleName,
        const std::string &barcodeName,
        common::TileNumber tileNumber,
        common::ReadNumber readNumber,
        const unsigned lane,
        const stats::ReadBarcodeStats& tileStats);

    void addFlowcellLane(
        const std::string &flowcellId,
        const unsigned lane,
        const stats::BarcodeHits::Popular &topUnknownBarcodes);

private:
    static boost::property_tree::path getPrefix(
        const std::string &flowcellId,
        const std::string &projectName,
        const std::string &sampleName,
        const std::string &barcodeName,
        common::TileNumber tileNumber,
        const unsigned lane );
};

inline std::ostream &operator<< (std::ostream &os, const ConversionStatsXml &tree)
{
    return bcl2fastq::io::serializeAsXml(os, tree);
}

inline std::istream &operator>> (std::istream &is, ConversionStatsXml &tree)
{
    return bcl2fastq::io::parseAsXml(is, dynamic_cast< boost::property_tree::ptree& >(tree) );
}


} //namespace stats
} //namespace bcl2fastq

#endif //BCL2FASTQ_STATS_CONVERSION_STATS_XML_HH
