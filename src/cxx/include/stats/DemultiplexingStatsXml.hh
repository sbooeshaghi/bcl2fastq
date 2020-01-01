/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file DemultiplexingStatsXml.hh
 *
 * \brief Xml Serialization of Demultiplexing statistics.
 *
 * \author Mauricio Varea
 */

#ifndef BCL2FASTQ_STATS_DEMULTIPLEXING_STATS_XML_HH
#define BCL2FASTQ_STATS_DEMULTIPLEXING_STATS_XML_HH

#include "io/Xml.hh"
#include "layout/Layout.hh"

#include "stats/BarcodeStats.hh"


namespace bcl2fastq
{
namespace stats
{

class DemultiplexingStatsXml : public boost::property_tree::ptree
{
public:
    DemultiplexingStatsXml();

    void addLaneBarcode(
        const std::string &flowcellId,
        const std::string &projectName,
        const std::string &sampleName,
        const std::string &barcodeName,
        const unsigned lane,
        const stats::BarcodeStats& barcodeStats);
};

inline std::ostream &operator<< (std::ostream &os, const DemultiplexingStatsXml &tree)
{
    return bcl2fastq::io::serializeAsXml(os, tree);
}

inline std::istream &operator>> (std::istream &is, DemultiplexingStatsXml &tree)
{
    return bcl2fastq::io::parseAsXml(is, dynamic_cast< boost::property_tree::ptree& >(tree) );
}


} //namespace stats
} //namespace bcl2fastq

#endif //BCL2FASTQ_STATS_DEMULTIPLEXING_STATS_XML_HH
