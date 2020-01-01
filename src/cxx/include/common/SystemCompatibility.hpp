/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file SystemCompatibility.hpp
 *
 * \brief Interface layer for system-dependent functionalities.
 *
 * \author Marek Balint
 */


#ifndef BCL2FASTQ_COMMON_SYSTEMCOMPATIBILITY_HPP
#define BCL2FASTQ_COMMON_SYSTEMCOMPATIBILITY_HPP


#include <algorithm>
#include <iterator>

#include "common/SystemCompatibility.hh"


namespace bcl2fastq {


namespace common {


template<typename T>
T swapEndian(T value)
{
    union {
        T value;
        unsigned char bytes[sizeof(T)];
    } src, dst;

    src.value = value;
    std::copy(
        &src.bytes[0],
        &src.bytes[sizeof(T)],
        std::reverse_iterator<unsigned char *>(&dst.bytes[sizeof(T)-1])
    );

    return dst.value;
}


} // namespace common


} // namespace bcl2fastq


#endif // BCL2FASTQ_COMMON_SYSTEMCOMPATIBILITY_HPP


