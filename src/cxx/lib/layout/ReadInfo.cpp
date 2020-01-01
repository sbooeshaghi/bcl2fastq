/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file ReadInfo.cpp
 *
 * \brief Implementation of read metadata.
 *
 * \author Marek Balint
 */


#include "layout/ReadInfo.hh"
#include "common/Exceptions.hh"

namespace bcl2fastq {


namespace layout {


ReadInfo::ReadInfo(common::ReadNumber        number,
                   common::ReadType          readType,
                   const common::CycleRange& cyclesToUse,
                   const common::CycleRange& unmaskedCycles,
                   const common::CycleRange& umiCycles)
: number_(number)
, readType_(readType)
, cycleInfos_()
, unmaskedCycleInfos_()
, cyclesToLoad_()
, maskAdapters_()
, trimAdapters_()
, umiCycleRange_(0, 0)
, bclBufferOffset_(0)
{
    BCL2FASTQ_LOG(common::LogLevel::INFO) << "Read: " << number_ << " " << getReadType() << std::endl;

    std::ostringstream cycles;
    for (common::CycleNumber cycleNumber = cyclesToUse.first; cycleNumber <= cyclesToUse.second; ++cycleNumber)
    {
        CycleInfo cycleInfo(cycleNumber);
        addCycle(cycleInfo);

        if (!cycles.str().empty())
        {
            cycles << ",";
        }
        cycles << cycleNumber;
    }
    for (common::CycleNumber cycleNumber = unmaskedCycles.first; cycleNumber <= unmaskedCycles.second; ++cycleNumber)
    {
        addUnmaskedCycle(CycleInfo(cycleNumber));
    }
    BCL2FASTQ_LOG(common::LogLevel::INFO) << "  Cycles: " << cycles.str() << std::endl;

    if (readType == common::ReadType::DATA)
    {
        addUmiBases(umiCycles);
    }

    if (!cycleInfos_.empty() && ! cyclesToLoad_.empty())
    {
        common::CycleNumber firstCycleToUse = cycleInfos_.front().getNumber();

        for (const auto& cycleToLoad : cyclesToLoad_)
        {
            if (cycleToLoad.getNumber() == firstCycleToUse)
            {
                break;
            }
            ++bclBufferOffset_;
        }
    }
}

void ReadInfo::addUmiBases(common::CycleRange umiCycleRange)
{
    if (umiCycleRange.second == 0)
    {
        return;
    }

    // Now holds relative start position and length
    umiCycleRange_ = umiCycleRange;

    if (!unmaskedCycleInfos_.empty())
    {
        // The CycleRange is the first and last cycle. Both read 1 and read 2
        // use 1 as the start number. Move to the actual cycle number for read 2.
        umiCycleRange.first += unmaskedCycleInfos_.front().getNumber() - 1;
        umiCycleRange.second += umiCycleRange.first - 1;
    }

    if (unmaskedCycleInfos_.empty() ||
        umiCycleRange.first < unmaskedCycleInfos_.front().getNumber() ||
        umiCycleRange.second > unmaskedCycleInfos_.back().getNumber())
    {
        BOOST_THROW_EXCEPTION(common::InputDataError("UMI bases are outside of the range of the bases for read " + boost::lexical_cast<std::string>(number_)));
    }

    if (umiCycleRange.second != 0)
    {
        BCL2FASTQ_LOG(common::LogLevel::INFO) << "Read: " << number_ << " " << "UMI" << std::endl;

        std::ostringstream cycles;
        for (common::CycleNumber cycle = umiCycleRange.first; cycle <= umiCycleRange.second; ++cycle)
        {
            if (!cycles.str().empty())
            {
                cycles << ",";
            }
            cycles << cycle;
            cyclesToLoad_.insert(layout::CycleInfo(cycle));
        }

        BCL2FASTQ_LOG(common::LogLevel::INFO) << "  Cycles: " << cycles.str() << std::endl;
    }

    // Everything up to here was validation and display.
    // Now we store the index into the bcl buffer
    common::CycleNumber firstUmiCycle = umiCycleRange_.first + unmaskedCycleInfos_.front().getNumber() - 1; 
    umiCycleRange_.first = 0;
    for (const auto& cycleToLoad : cyclesToLoad_)
    {
        if (cycleToLoad.getNumber() == firstUmiCycle)
        {
            break;
        }
        ++umiCycleRange_.first;
    }

    umiCycleRange_.second += umiCycleRange_.first;

    // umiCycleRange now holds the absolute index into the bcl buffer
}

void ReadInfo::addUnmaskedCycle(const CycleInfo &cycleInfo)
{
    unmaskedCycleInfos_.push_back(cycleInfo);
    if (isIndexRead())
    {
        // We need to load all the cycles from BCL file for index reads
        cyclesToLoad_.insert(cycleInfo);
    }
}

} // namespace layout


} // namespace bcl2fastq


