/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file FastIo.hpp
 *
 * \brief Implementation of fast IO routines for integers and fixed width floating points.
 *
 * \author Come Raczy
 * \author Marek Balint
 */


#ifndef BCL2FASTQ_COMMON_FASTIO_HPP
#define BCL2FASTQ_COMMON_FASTIO_HPP


#include "common/FastIo.hh"


namespace bcl2fastq {


namespace common {


template<typename T>
void putUnsignedInteger(T value, std::vector<char>& copyBuffer)
{
    // generate a compilation error when instantiated on signed types
    BOOST_MPL_ASSERT_MSG(boost::is_unsigned<T>::value,
            SIGNED_TYPES_ARE_NOT_ALLOWED_FOR_putUnsignedInteger, (T));

    char begin[1 + std::numeric_limits<T>::digits10];
    char *buffer = begin;
    while (value >= 10)
    {
        *buffer++ = '0' + (value % 10);
        value /= 10;
    }
    *buffer++ = '0' + value;
    std::reverse(begin, buffer);

    size_t bufferSize = copyBuffer.size();
    copyBuffer.resize(bufferSize + buffer-begin);
    std::copy(begin, buffer, copyBuffer.begin() + bufferSize);
}


} // namespace common


} // namespace bcl2fastq


#endif // BCL2FASTQ_COMMON_FASTIO_HPP


