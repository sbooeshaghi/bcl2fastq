/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file Debug.hh
 *
 * \brief Various debugging-related helpers.
 *
 * \author Roman Petrovski
 * \author Marek Balint
 */


#ifndef BCL2FASTQ_COMMON_DEBUG_HH
#define BCL2FASTQ_COMMON_DEBUG_HH


#include "common/Logger.hh"
#include "common/SystemCompatibility.hh"


/// \brief Asserts condition (always, i.e. even if NDEBUG is set and so on).
/// \param condition Condition to be asserted.
/// \param message Error message.
/// \note Also uses ostream serialization which, unlike the standard assert,
/// has shown not to allocate the dynamic memory at the time when you least
/// expect this to happen.
#define BCL2FASTQ_ASSERT_MSG(condition, message)                \
    {                                                           \
        if (condition)                                          \
        {                                                       \
            /* nothing here */                                  \
        }                                                       \
        else                                                    \
        {                                                       \
            BCL2FASTQ_LOG(::bcl2fastq::common::LogLevel::FATAL) \
                << "***** Internal Program Error - assertion (" \
                << #condition                                   \
                << ") failed in "                               \
                << __FILE__ << ":"                              \
                << __LINE__ << ":"                              \
                << (BOOST_CURRENT_FUNCTION) << ": "             \
                << message                                      \
                << std::endl                                    \
            ;                                                   \
            ::bcl2fastq::common::terminateWithCoreDump();       \
        }                                                       \
    }


namespace bcl2fastq {


namespace common {


} // namespace common


} // namespace bcl2fastq


#endif // BCL2FASTQ_COMMON_DEBUG_HH


