/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file FastIo.hh
 *
 * \brief Definition fast IO routines for integers and fixed width floating points.
 *
 * \author Come Raczy
 * \author Marek Balint
 */


#ifndef BCL2FASTQ_COMMON_FASTIO_HH
#define BCL2FASTQ_COMMON_FASTIO_HH


#include <algorithm>
#include <limits>
#include <vector>

#include <boost/type_traits.hpp>
#include <boost/mpl/assert.hpp>


namespace bcl2fastq {


namespace common {


/// \brief Fast and portable output of an unsigned integer into a stream.
/// 
/// Drops all the usual formatting options to the benefit of speed.
/// 
/// Can be instantiated on unsigned versions of char, short, int and
/// long, and generally on any type supporting '<=', '%', '/=', '+'
/// and defining 'digits10 in std::numeric_limits.
template<typename T>
void putUnsignedInteger(T value, std::vector<char>& copyBuffer);


} // namespace common


} // namespace bcl2fastq


#include "common/FastIo.hpp"


#endif // BCL2FASTQ_COMMON_FASTIO_HH


