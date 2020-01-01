/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file CycleInfo.hh
 *
 * \brief Declaration of cycle metadata.
 *
 * \author Marek Balint
 */


#ifndef BCL2FASTQ_LAYOUT_CYCLEINFO_HH
#define BCL2FASTQ_LAYOUT_CYCLEINFO_HH


#include <vector>

#include <boost/filesystem/path.hpp>

#include "common/Types.hh"


namespace bcl2fastq {


namespace layout {


/// \brief Cycle metadata.
class CycleInfo
{
public:

    /// \brief BCL file metadata container type definition.
    typedef std::vector<boost::filesystem::path> BclFilesContainer;

public:

    /// \brief Constructor.
    /// \param number Cycle number.
    explicit CycleInfo(common::CycleNumber number);

public:

    /// \brief Get cycle number.
    /// \return Cycle number.
    common::CycleNumber getNumber() const;

    /// \brief Equality operator.
    /// \return true if equal.
    bool operator==(const CycleInfo& other) const { return number_ == other.number_; }

    /// \brief Less than operator
    /// \return true if less.
    bool operator<(const CycleInfo& other) const { return number_ < other.number_; }

private:

    /// \brief Cycle number.
    /// Source: RunInfo.xml (FirstCycle..LastCycle or NumCycles)
    common::CycleNumber number_;
};


} // namespace layout


} // namespace bcl2fastq


#endif // BCL2FASTQ_LAYOUT_CYCLEINFO_HH


