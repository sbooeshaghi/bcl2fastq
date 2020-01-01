/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 * 
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file UseBasesMask.hh
 *
 * \brief Declaration of UseBasesMask.
 *
 * \author Aaron Day
 */

#ifndef BCL2FASTQ_LAYOUT_USE_BASES_MASK_HH
#define BCL2FASTQ_LAYOUT_USE_BASES_MASK_HH

#include "common/Types.hh"

#include "layout/ReadMetadata.hh"
#include "common/Exceptions.hh"

#include <string>

namespace bcl2fastq {


namespace layout {

class UseBasesMaskFormatError : public common::RuntimeError
{
public:
    UseBasesMaskFormatError(const std::string& message);
    UseBasesMaskFormatError(const UseBasesMaskFormatError& other);
};

class UseBasesMask
{
public:
    UseBasesMask(std::string                           mask,
                 ReadMetadataContainer::const_iterator readsBegin,
                 ReadMetadataContainer::const_iterator readsEnd);

    bool isEmpty() const { return reads_.empty(); }

    /// \brief Get beginning of read meta data.
    /// \return Iterator to beginning of read meta data.
    ReadMetadataContainer::const_iterator readMetadataBegin() const;

    /// \brief Get end of read meta data.
    /// \return Iterator to end of read meta data.
    ReadMetadataContainer::const_iterator readMetadataEnd() const;
private:

    /// \brief Returns true if the character is invalid
    bool isMaskCharInvalid(char maskChar) const; 

    /// \brief Throw an exception if the mask does not meet basic validation criteria.
    void basicValidate(const std::string& mask) const;

    /// \brief Expand the numbers in the mask, i.e. y3 => yyy
    void expandNumbers(std::string& mask) const;

    /// \brief Parse the asterisk from the mask
    void parseAsteriskAndValidateSize(std::string& mask,
                                      bool&        doesMatchRunInfo,
                                      unsigned int totalCycles) const;

    /// \brief Parse the cycles to ignore.
    /// \return pair<# ignore beginning, # ignore end>
    std::pair<int, int> parseIgnoreCycles(std::string& mask) const;

    /// \brief Add the cycles to the reads_.
    void addCycles(const std::string&   mask,
                   common::ReadType     readType,
                   common::ReadNumber&  dataReadsCounter,
                   common::ReadNumber&  indexReadsCounter,
                   common::CycleNumber& cyclesCounter);

    /// \brief Read meta data container.
    ReadMetadataContainer reads_;
};

} // namespace layout

} // namespace bcl2fastq

#endif // BCL2FASTQ_LAYOUT_USE_BASES_MASK_HH

