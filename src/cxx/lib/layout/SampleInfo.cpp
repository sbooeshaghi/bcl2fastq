/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file SampleInfo.cpp
 *
 * \brief Implementation of sample metadata.
 *
 * \author Marek Balint
 */


#include "layout/SampleInfo.hh"


namespace bcl2fastq
{
namespace layout
{


const std::string SampleInfo::defaultProject = "default";
const std::string SampleInfo::defaultId      = "unknown";
const std::string SampleInfo::defaultName    = "Undetermined";

SampleInfo::SampleInfo(
    common::SampleNumber number,
    SampleId sampleId,
    SampleName sampleName,
    Project project
)
: number_(number)
, sampleId_(sampleId)
, sampleName_(sampleName)
, project_(project)
, barcodes_()
{
}

common::SampleNumber SampleInfo::getNumber() const
{
    return number_;
}

const Barcode& SampleInfo::getBarcode(BarcodesContainer::size_type idx) const
{
    return barcodes_[idx];
}

bool SampleInfo::hasId() const
{
    return SampleInfo::defaultId != sampleId_;
}

bool SampleInfo::hasName() const
{
    return SampleInfo::defaultName != sampleName_;
}

bool SampleInfo::hasProject() const
{
    return SampleInfo::defaultProject != project_;
}

bool SampleInfo::hasBarcodes() const
{
    return !barcodes_.empty();
}

void SampleInfo::addBarcode(const Barcode &barcode)
{
    barcodes_.push_back(barcode);
}

void SampleInfo::maskBarcode(common::ReadNumber indexNumber)
{
    bool noComponents = false;
    for (auto& barcode : barcodes_)
    {
        barcode.mask(indexNumber);
        if (barcode.getComponents().empty())
        {
            noComponents = true;
            break;
        }
    }

    if (noComponents)
    {
        barcodes_.clear();
    }
}

std::ostream& operator<<(std::ostream& os, const SampleInfo &sampleInfo)
{
    os << "#" << sampleInfo.number_
       << " '" << std::string(sampleInfo.sampleId_) << "' '" << std::string(sampleInfo.sampleName_)
       << "' [" << std::string(sampleInfo.project_) << "]";
    return os;
}


} // namespace layout
} // namespace bcl2fastq


