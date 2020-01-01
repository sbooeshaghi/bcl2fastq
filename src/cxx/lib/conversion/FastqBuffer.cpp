/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file FastqBuffer.cpp
 *
 * \brief Implementation of FASTQ buffer.
 *
 * \author Marek Balint
 * \author Aaron Day
 */


#include "conversion/FastqBuffer.hh"


namespace bcl2fastq {
namespace conversion {


FastqBuffer::FastqBuffer()
: fastqs_()
, groupNumber_(0)
{
}


} // namespace conversion
} // namespace bcl2fastq

