/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file \file DemultiplexingStatsXml.cpp
 *
 * \brief Xml Serialization of Demultiplexing statistics.
 *
 * \author Mauricio Varea
 */

#include <boost/lexical_cast.hpp>

#include "stats/DemultiplexingStatsXml.hh"

namespace bcl2fastq
{
namespace stats
{

DemultiplexingStatsXml::DemultiplexingStatsXml()
{
}


void DemultiplexingStatsXml::addLaneBarcode(
    const std::string &flowcellId,
    const std::string &projectName,
    const std::string &sampleName,
    const std::string &barcodeName,
    const unsigned lane,
    const stats::BarcodeStats& barcodeStats)
{
    std::string matchCountStrs[stats::BarcodeStats::MAX_MISMATCHES+1] =
        {"PerfectBarcodeCount",
         "OneMismatchBarcodeCount",
         "TwoMismatchesBarcodeCount",
         "ThreeMismatchesBarcodeCount",
         "FourMismatchesBarcodeCount",
         "FiveMismatchesBarcodeCount"};

    const boost::property_tree::path tileValuePrefix("Stats"
                                      "/<indexed>Flowcell/<flowcell-id>" + flowcellId
                                      +"/<indexed>Project/<name>" + projectName
                                      +"/<indexed>Sample/<name>" + sampleName
                                      +"/<indexed>Barcode/<name>" + barcodeName
                                      +"/<indexed>Lane/<number>" + boost::lexical_cast<std::string>(lane)
                                      , '/');

    add(tileValuePrefix / "BarcodeCount", barcodeStats.getBarcodeCount());

    auto pos = std::find_if(barcodeStats.getBarcodeMismatchCounts().rbegin(),
                            barcodeStats.getBarcodeMismatchCounts().rend(),
                            [](common::ClustersCount count) {return count != 0; });

    size_t numMismatches = 0;
    for (auto iter = barcodeStats.getBarcodeMismatchCounts().begin(); iter != pos.base(); ++iter)
    {
        add(tileValuePrefix / matchCountStrs[numMismatches++].c_str(), *iter);
    }
}


} //namespace stats
} //namespace bcl2fastq

