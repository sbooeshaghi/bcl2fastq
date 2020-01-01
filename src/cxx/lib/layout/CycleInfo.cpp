/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file CycleInfo.cpp
 *
 * \brief Implementation of cycle metadata.
 *
 * \author Marek Balint
 */


#include <sstream>

#include <boost/format.hpp>

#include "layout/CycleInfo.hh"


namespace bcl2fastq {


namespace layout {


CycleInfo::CycleInfo(common::CycleNumber number)
: number_(number)
{
}

common::CycleNumber CycleInfo::getNumber() const
{
    return number_;
}

} // namespace layout


} // namespace bcl2fastq


