/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file Endianness.hh
 *
 * \brief Fast IO routines for integers and fixed width floating points.
 *
 * \author Come Raczy
 */

#ifndef BCL2FASTQ_COMMON_ENDIANNESS_HH
#define BCL2FASTQ_COMMON_ENDIANNESS_HH


namespace bcl2fastq
{
namespace common
{

template <typename T, typename IterT> IterT extractLittleEndian(const IterT p, T &result)
{
#ifdef BCL2FASTQ_IS_BIG_ENDIAN
    // TODO: untested
    BOOST_STATIC_ASSERT(sizeof(T) / 8 == 4);
    // this should work on any endian machine but if we're on little-endian already the alternative is probably faster
    result = T(*p) + T(*(p + 1)) * 256 + T(*(p + 2)) * 256 * 256 + T(*(p + 3)) * 256 * 256 * 256;
#else
    result = *reinterpret_cast<const T*>(&*p);
#endif
    return p + sizeof(T);
}

template <typename T, typename IterT> T extractLittleEndian(const IterT p)
{
#ifdef BCL2FASTQ_IS_BIG_ENDIAN
    // TODO: untested
    BOOST_STATIC_ASSERT(sizeof(T) / 8 == 4);
    // this should work on any endian machine but if we're on little-endian already the alternative is probably faster
    return T(*p) + T(*(p + 1)) * 256 + T(*(p + 2)) * 256 * 256 + T(*(p + 3)) * 256 * 256 * 256;
#endif
    return *reinterpret_cast<const T*>(&*p);
}

} // namespace common
} // namespace bcl2fastq

#endif // #ifndef BCL2FASTQ_COMMON_ENDIANNESS_HH
