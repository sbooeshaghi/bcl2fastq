/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file Timer.hh
 *
 * \brief Declaration of common data types
 *
 * \author Aaron Day
 */


#ifndef BCL2FASTQ_COMMON_TIMER_HH
#define BCL2FASTQ_COMMON_TIMER_HH


#include <chrono>

#include <boost/noncopyable.hpp>

namespace bcl2fastq {
namespace common {

class Timer : private boost::noncopyable
{
public:
    Timer(size_t& timeInMicroseconds) :
        timeInMicroseconds_(timeInMicroseconds),
        startTime_(std::chrono::steady_clock::now())
    {
    }

    ~Timer()
    {
        timeInMicroseconds_ += std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - startTime_).count();
    }

private:
    size_t& timeInMicroseconds_;
    std::chrono::time_point<std::chrono::steady_clock> startTime_;
};

} // namespace common
} // namespace bcl2fastq


#endif // BCL2FASTQ_DATA_TIMER_HH


