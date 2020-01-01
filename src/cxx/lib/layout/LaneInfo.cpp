/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file LaneInfo.cpp
 *
 * \brief Implementation of lane metadata.
 *
 * \author Marek Balint
 */


#include <sstream>
#include <iterator>
#include <numeric>
#include <utility>

#include <boost/assign.hpp>
#include <boost/format.hpp>
#include <boost/foreach.hpp>

#include "common/Debug.hh"
#include "layout/LaneInfo.hh"


namespace bcl2fastq {


namespace layout {


LaneInfo::LaneInfo(common::LaneNumber number)
: number_(number)
, sampleInfos_()
, tileInfos_()
, readInfos_()
, minimumTrimmedReadLength_(0)
{
}

void LaneInfo::maskBarcode(common::ReadNumber indexNumber)
{
    for (auto& sampleInfo : sampleInfos_)
    {
        sampleInfo.maskBarcode(indexNumber);
    }
}

void LaneInfo::removeDefaultSample()
{
    sampleInfos_.erase(
        std::remove_if(sampleInfos_.begin(),
                       sampleInfos_.end(),
                       [](const SampleInfo& sample) { return sample.getNumber() == 0; }));
}

void LaneInfo::sortSampleInfos()
{
    std::sort(sampleInfos_.begin(),
              sampleInfos_.end(),
              [] (const SampleInfo& sample1,
                  const SampleInfo& sample2)
                 { return sample1.getNumber() < sample2.getNumber(); });
}

common::TotalClustersCount LaneInfo::getClustersCount() const
{
    return std::accumulate(
        tileInfos_.begin(),
        tileInfos_.end(),
        common::TotalClustersCount(0),
        boost::bind(
            std::plus<common::TotalClustersCount>(),
            _1,
            boost::bind(&TileInfo::getClustersCount, _2)
        )
    );
}

bool LaneInfo::haveClustersCount() const
{
    for (const auto &tileInfo : tileInfos_)
    {
        if (!tileInfo.haveClustersCount())
        {
            return false;
        }
    }
    return true;
}

common::CycleNumber LaneInfo::getNumCyclesToLoad() const
{
    return std::accumulate(
        readInfos_.begin(),
        readInfos_.end(),
        common::CycleNumber(0),
        [](common::CycleNumber sum, const layout::ReadInfo& elem) { return sum + elem.getNumCyclesToLoad(); });
}



} // namespace layout


} // namespace bcl2fastq


