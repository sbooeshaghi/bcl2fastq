/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file UseBasesMask.cpp
 *
 * \brief Implementation of UseBasesMask.
 *
 * \author Aaron Day
 */

#include "layout/UseBasesMask.hh"

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>
#include <vector>
#include <functional>

namespace bcl2fastq {

namespace layout {

UseBasesMaskFormatError::UseBasesMaskFormatError(const std::string& message)
: common::RuntimeError(0, message)
{
}

UseBasesMaskFormatError::UseBasesMaskFormatError(const UseBasesMaskFormatError& other)
: common::RuntimeError(other)
{
}

UseBasesMask::UseBasesMask(std::string                           baseMasks,
                           ReadMetadataContainer::const_iterator readsBegin,
                           ReadMetadataContainer::const_iterator readsEnd) :
    reads_()
{
    if (baseMasks.empty()) return;

    try
    {
        boost::algorithm::to_lower(baseMasks);

        std::vector<std::string> masks;
        boost::split(masks, baseMasks, boost::is_any_of(","));

        common::ReadNumber indexReadsCounter = 1;
        common::ReadNumber dataReadsCounter = 1;
        common::CycleNumber cyclesCounter = 1;

        if (masks.size() != static_cast<size_t>(std::distance(readsBegin, readsEnd)))
        {
            BOOST_THROW_EXCEPTION(UseBasesMaskFormatError(
                "UseBasesMask formatting error. A mask must be specified for each read."
                " Number of reads: " + boost::lexical_cast<std::string>(std::distance(readsBegin, readsEnd))
                + " Base masks: '" + baseMasks + "'"));
        }

        size_t numCyclesRunInfoAllReads = 0;
        bool doesMatchRunInfo = true;
        ReadMetadataContainer::const_iterator readIter = readsBegin;
        BOOST_FOREACH(std::string& mask, std::make_pair(masks.begin(), masks.end()))
        {
            unsigned int totalCyclesRunInfoXml = readIter->lastCycle_ - readIter->firstCycle_ + 1;
            numCyclesRunInfoAllReads += totalCyclesRunInfoXml;

            // Get rid of white space
            boost::trim(mask);

            basicValidate(mask);

            expandNumbers(mask);

            parseAsteriskAndValidateSize(mask,
                                         doesMatchRunInfo,
                                         totalCyclesRunInfoXml);

            std::pair<int, int> numIgnoreCycles = parseIgnoreCycles(mask);
            cyclesCounter += numIgnoreCycles.first;
            addCycles(mask,
                      readIter->readType_,
                      dataReadsCounter,
                      indexReadsCounter,
                      cyclesCounter);
            cyclesCounter += numIgnoreCycles.second;

            // Go to the next read in RunInfo.xml
            ++readIter;
        }

        if (cyclesCounter-1 != numCyclesRunInfoAllReads)
        {
            BOOST_THROW_EXCEPTION(UseBasesMaskFormatError("The total number of cycles specified by"
                " the use-bases-mask does not match the number found in RunInfo.xml. Masks: " + baseMasks + "'"));
        }
    }
    catch(boost::bad_lexical_cast&)
    {
        BOOST_THROW_EXCEPTION(UseBasesMaskFormatError("UseBasesMask formatting error. Failed attempt at reading number of cycles"));
    }
}

ReadMetadataContainer::const_iterator UseBasesMask::readMetadataBegin() const
{
    return reads_.begin();
}

ReadMetadataContainer::const_iterator UseBasesMask::readMetadataEnd() const
{
    return reads_.end();
}

bool UseBasesMask::isMaskCharInvalid(char maskChar) const
{
    return maskChar != 'i' &&
           maskChar != 'y' &&
           maskChar != 'n';
}

// Throw an exception if the mask is not in a valid format
void UseBasesMask::basicValidate(const std::string& mask) const
{
    if (mask.size() < 2)
    {
        BOOST_THROW_EXCEPTION(UseBasesMaskFormatError(
            "UseBasesMask formatting error. Length too small."
            " Use base mask: '" + mask + "'"));
    }

    size_t invalidCharPos = mask.find_first_not_of("*yni0123456789");
    if (invalidCharPos != std::string::npos)
    {
        BOOST_THROW_EXCEPTION(UseBasesMaskFormatError(
               std::string("UseBasesMask formatting error. Invalid character: '")
               + mask[invalidCharPos] + std::string("' in mask: '") + mask + "'"));
    }

    if (isMaskCharInvalid(mask[0]))
    {
        BOOST_THROW_EXCEPTION(UseBasesMaskFormatError(
            "UseBasesMask formatting error. Use base mask must start with either a 'y', 'i', or 'n'."
            " Use base mask: '" + mask + "'"));
    }

    size_t asteriskPos = mask.find_first_of('*');
    if (asteriskPos != std::string::npos)
    {
        if (mask.find('*', asteriskPos+1) != std::string::npos)
        {
            BOOST_THROW_EXCEPTION(UseBasesMaskFormatError(
                "UseBasesMask formatting error. Use base mask must contain at most 1 '*'."
                " Use base mask: '" + mask + "'"));
        }

        if (isMaskCharInvalid(mask[asteriskPos-1]))
        {
            BOOST_THROW_EXCEPTION(UseBasesMaskFormatError(
                "UseBasesMask formatting error. An asterisk must be preceded by either a 'y', 'i', or 'n'."
                " Use base mask: '" + mask + "'"));
        }
            
    }
}

void UseBasesMask::expandNumbers(std::string& mask) const
{
    while (true)
    {
        size_t beginDigitPos = mask.find_first_of("0123456789");
        if (beginDigitPos == std::string::npos) return;

        if (beginDigitPos == 0 || mask[beginDigitPos] == '0')
        {
            BOOST_THROW_EXCEPTION(UseBasesMaskFormatError(
                "UseBasesMask formatting error. Error parsing number from mask."
                " Base mask: '" + mask + "'"));
        }

        size_t endDigitPos = mask.find_first_not_of("0123456789", beginDigitPos);
        size_t digitsLength = endDigitPos==std::string::npos ?
                                  mask.size()-beginDigitPos :
                                  endDigitPos-beginDigitPos;

        unsigned int numExpand = boost::lexical_cast<unsigned int>(mask.substr(beginDigitPos, digitsLength));

        mask.erase(beginDigitPos, digitsLength);

        // numExpand-1 because y2 should be replaced by yy (just 1 y added)
        mask.insert(beginDigitPos, numExpand-1, mask[beginDigitPos-1]);
    }
}

void UseBasesMask::parseAsteriskAndValidateSize(std::string& mask,
                                                bool&        doesMatchRunInfo,
                                                unsigned int totalCycles) const
{
    size_t asteriskPos = mask.find('*');

    if (asteriskPos == std::string::npos)
    {
        if (mask.size() != totalCycles)
        {
            doesMatchRunInfo = false;
        }

        return;
    }

    if (!doesMatchRunInfo)
    {
        BOOST_THROW_EXCEPTION(UseBasesMaskFormatError("One or more read masks do not match"
            " the number of cycles specified in RunInfo.xml, and a '*' was specified."
            " There is no way to determine the number of cycles implied by the '*'."
            " Base mask: " + mask + "'"));
    }

    if (mask.size() > totalCycles + 2)
    {
        BOOST_THROW_EXCEPTION(UseBasesMaskFormatError(
            "UseBasesMask formatting error. Mask size does not match number of cycles in RunInfo.xml."
            " RunInfo.xml cycles: " + boost::lexical_cast<std::string>(totalCycles) +
            " Base mask: '" + mask + "'"));
    }

    if (asteriskPos == 0 ||
        (mask[asteriskPos - 1] != 'i' &&
         mask[asteriskPos - 1] != 'n' &&
         mask[asteriskPos - 1] != 'y'))
    {
        BOOST_THROW_EXCEPTION(UseBasesMaskFormatError(
            "UseBasesMask formatting error.  Invalid use of '*'."
            " Base mask: '" + mask + "'"));
    }

    // If we got here, the format is valid.
    mask.erase(asteriskPos, 1);

    int numExpandCycles = totalCycles - mask.size();
    if (numExpandCycles >= 0)
    {
        mask.insert(asteriskPos, numExpandCycles, mask[asteriskPos-1]);
    }
    else
    {
        BCL2FASTQ_ASSERT_MSG(numExpandCycles == -1, "Use bases mask internal error");
        mask.resize(mask.size()-1);
    }
}

// Cut the n's off the mask. Return (# n's at beginning, #n's at end)
std::pair<int, int> UseBasesMask::parseIgnoreCycles(std::string& mask) const
{
    size_t firstNotN = mask.find_first_not_of('n');

    if (firstNotN == std::string::npos)
    {
        std::pair<int, int> ret(mask.size(), 0);
        mask = "";
        return ret;
    }

    size_t lastNotN = mask.find_last_not_of('n');

    std::string maskWithoutN = mask.substr(firstNotN, lastNotN-firstNotN+1);

    if (maskWithoutN.find('n') != std::string::npos)
    {
        BOOST_THROW_EXCEPTION(UseBasesMaskFormatError(
            "UseBasesMask formatting error.  Masking out cycles in the middle of the read is not allowed."
            " Base mask: '" + mask + "'"));
    }

    int numIgnoreAtEnd = mask.size() - lastNotN - 1;
    mask = maskWithoutN;

    return std::make_pair(firstNotN, numIgnoreAtEnd);
}


void UseBasesMask::addCycles(const std::string&   mask,
                             common::ReadType     readType,
                             common::ReadNumber&  dataReadsCounter,
                             common::ReadNumber&  indexReadsCounter,
                             common::CycleNumber& cyclesCounter)
{
    size_t nextReadBeginPos = 0;

    if (mask.find_first_not_of('n') == std::string::npos)
    {
        // The read is all 'n's, all bases are masked.
        switch (readType)
        {
            case common::ReadType::DATA:
                ++dataReadsCounter;
                break;
            case common::ReadType::INDEX:
                ++indexReadsCounter;
                break;
            default:
                BCL2FASTQ_ASSERT_MSG(false, "Unrecognized read type");
                break;
        }
        return;
    }

    while (nextReadBeginPos != std::string::npos)
    {
        size_t readBeginPos = nextReadBeginPos;

        reads_.push_back(ReadMetadata());
        ReadMetadata& readMetadata = reads_.back();

        size_t firstNotI = mask.find_first_not_of('i', readBeginPos);

        if (firstNotI == readBeginPos)
        {
            readMetadata.readNumber_ = dataReadsCounter++;
            readMetadata.readType_ = common::ReadType::DATA;
            nextReadBeginPos = mask.find_first_not_of('y', readBeginPos);
        }
        else
        {
            readMetadata.readNumber_ = indexReadsCounter++;
            readMetadata.readType_ = common::ReadType::INDEX;
            nextReadBeginPos = mask.find_first_not_of('i', readBeginPos);
        }

        readMetadata.firstCycle_ = cyclesCounter;
        readMetadata.lastCycle_ = readMetadata.firstCycle_ +
            (nextReadBeginPos == std::string::npos ?
             mask.size() - readBeginPos :
             nextReadBeginPos - readBeginPos) - 1;

        // I have no idea how to figure out which cycles are masked/unmasked
        // when the UseBasesMask could alter the read layout entirely.
        readMetadata.firstUnmaskedCycle_ = readMetadata.firstCycle_;
        readMetadata.lastUnmaskedCycle_ = readMetadata.lastCycle_;

        cyclesCounter += readMetadata.lastCycle_ - readMetadata.firstCycle_ + 1;
    }
}

}

}
