/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file FastqBuffer.hh
 *
 * \brief Declaration of FASTQ buffer.
 *
 * \author Marek Balint
 * \author Aaron Day
 */


#ifndef BCL2FASTQ_CONVERSION_FASTQBUFFER_HH
#define BCL2FASTQ_CONVERSION_FASTQBUFFER_HH

#include <stdint.h>
#include <vector>

namespace bcl2fastq {
namespace conversion {


/// \brief FASTQ buffer.
struct FastqBuffer
{
public:
    FastqBuffer();

    /// \brief FASTQ data container type definition.
    /// \note Container is indexed by sample, read and task
    /// (task is needed for parallelization).
    typedef std::vector<std::vector<std::vector<std::vector<char> > > > FastqsContainer;

public:

    /// \brief Buffers for FASTQ data.
    FastqsContainer fastqs_;
    uint32_t groupNumber_;

    uint32_t getGroupNumber() const { return groupNumber_; }
    void setGroupNumber(uint32_t groupNumber) { groupNumber_ = groupNumber; }
};

} // namespace conversion
} // namespace bcl2fastq


#endif // BCL2FASTQ_CONVERSION_FASTQBUFFER_HH

