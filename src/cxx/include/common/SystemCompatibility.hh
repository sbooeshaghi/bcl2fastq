/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file SystemCompatibility.hh
 *
 * \brief Interface layer for system-dependent functionalities.
 *
 * \author Come Raczy
 * \author Marek Balint
 */


#ifndef BCL2FASTQ_COMMON_SYSTEMCOMPATIBILITY_HH
#define BCL2FASTQ_COMMON_SYSTEMCOMPATIBILITY_HH

#include <thread>

namespace bcl2fastq {
namespace common {


/// \brief Generate a core dump with a meaningful backtrace.
/// \note This function is guaranteed to never return.
void terminateWithCoreDump();

/// \brief Determine architecture's endianness.
/// \retval true Little endian.
/// \retval false Big endian.
bool isLittleEndian();

/// \brief Swap endianness.
/// \param value Value which endianness is to be swapped.
/// \return Value with swapped endianness.
template<typename T>
T swapEndian(T value);

/// \brief Sets the max number of file handles
void adjustMaxFileHandles();

} // namespace common
} // namespace bcl2fastq


#include "common/SystemCompatibility.hpp"


#endif // BCL2FASTQ_COMMON_SYSTEMCOMPATIBILITY_HH


