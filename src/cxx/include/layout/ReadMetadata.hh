/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 * 
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file ReadMetadata.hh
 *
 * \brief Declaration of ReadMetadata.
 *
 * \author Marek Balint
 * \author Mauricio Varea
 * \author Aaron Day
 */

#ifndef BCL2FASTQ_LAYOUT_READ_METADATA_HH
#define BCL2FASTQ_LAYOUT_READ_METADATA_HH

#include "common/Types.hh"

#include <vector>

namespace bcl2fastq {


namespace layout {

    struct ReadMetadata
    {
        typedef std::pair<common::CycleNumber,common::CycleNumber> Range;

        ReadMetadata( Range r = Range(0, 0), common::ReadNumber readNumber = 0, common::ReadType readType = common::ReadType::DATA )
        : readNumber_(readNumber), readType_(readType), firstCycle_(r.first), lastCycle_(r.second), firstUnmaskedCycle_(r.first), lastUnmaskedCycle_(r.second) {}
    public:

        /// \brief Read number.
        common::ReadNumber readNumber_;

        /// \brief Read type (DATA, INDEX, or UMI).
        common::ReadType readType_;

        /// \brief Number of the first cycle.
        common::CycleNumber firstCycle_;

        /// \brief Number of the last cycle.
        common::CycleNumber lastCycle_;

        /// \brief Number of first unmasked cycle.
        common::CycleNumber firstUnmaskedCycle_;

        /// \brief Number of last unmasked cycle.
        common::CycleNumber lastUnmaskedCycle_;
    };

    /// \brief Read meta data container type definition.
    typedef std::vector<ReadMetadata> ReadMetadataContainer;

} // namespace layout

} // namespace bcl2fastq

#endif // BCL2FASTQ_LAYOUT_READ_METADATA_HH

