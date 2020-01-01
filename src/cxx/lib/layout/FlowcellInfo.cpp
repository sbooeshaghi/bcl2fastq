/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file FlowcellInfo.cpp
 *
 * \brief Implementation of flowcell metadata.
 *
 * \author Marek Balint
 * \author Mauricio Varea
 */


#include <sstream>

#include <boost/format.hpp>

#include "layout/FlowcellInfo.hh"


namespace bcl2fastq {


namespace layout {


void FlowcellInfo::setInstrument(const std::string &instrument)
{
    instrument_ = instrument;
}

void FlowcellInfo::setRunNumber(const std::string &runNumber)
{
    runNumber_ = runNumber;
}

void FlowcellInfo::setFlowcellId(const std::string &flowcellId)
{
    flowcellId_ = flowcellId;
}

common::TileAggregationMode FlowcellInfo::getAggregateTilesMode() const
{
    return aggregateTilesFlag_;
}

void FlowcellInfo::setAggregateTilesFlag(common::TileAggregationMode aggregateTilesFlag)
{
    aggregateTilesFlag_ = aggregateTilesFlag;
}

bool FlowcellInfo::isPatternedFlowcell() const
{
    return patternedFlowcell_;
}

void FlowcellInfo::setPatternedFlowcell()
{
	patternedFlowcell_ = true;
}

void FlowcellInfo::resetPatternedFlowcell()
{
	patternedFlowcell_ = false;
}


} // namespace layout


} // namespace bcl2fastq


