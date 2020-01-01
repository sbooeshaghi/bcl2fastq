/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file SystemCompatibility.cpp
 *
 * \brief Interface layer for system-dependent functionalities.
 *
 * \author Come Raczy
 * \author Marek Balint
 */


#include <signal.h>

#include "common/SystemCompatibility.hh"
#include "common/Logger.hh"

#ifdef WIN32
#include <Windows.h>
#endif

namespace bcl2fastq {
namespace common {


void terminateWithCoreDump()
{
    ::raise(SIGSEGV);
}


bool isLittleEndian()
{
    const unsigned long v = 0x0706050403020100;
    const unsigned char * const p = reinterpret_cast<const unsigned char *>(&v);
    for (unsigned i = 0; i < sizeof(v); ++i)
    {
        if (p[i] != i)
        {
            return false;
        }
    }
    return true;
}


void adjustMaxFileHandles()
{
#ifdef WIN32
    // increase the maximum number of file handles
    _setmaxstdio(2048);
#endif
}

} // namespace common
} // namespace bcl2fastq


