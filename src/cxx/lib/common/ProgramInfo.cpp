/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file ProgramInfo.cpp
 *
 * \brief Basic program information.
 *
 * \author Mauricio Varea
 * \author Marek Balint
 */


#include <iostream>
#include <iomanip>

#include "common/ProgramInfo.hh"
#include "config.h"


namespace bcl2fastq {


namespace common {


ProgramInfo::ProgramInfo(const std::string &cmd, unsigned int width)
: command(cmd)
, binaryName(this->command.filename().string())
, projectNameShort(BCL2FASTQ_NAME_SHORT)
, projectNameLong(BCL2FASTQ_NAME_LONG)
, copyright(BCL2FASTQ_COPYRIGHT)
, version(BCL2FASTQ_VERSION)
, width_(width)
{
}


namespace {


/// \brief Output helper.
/// \param os Output stream to print to.
/// \param str String to be printed out.
/// \param width Line width.
static void outputLine(std::ostream& os, const std::string &str, unsigned int width)
{
    os
        //<< std::setfill(' ') << std::setw((width - str.length()) / 2)
        << str << std::endl
    ;
}


} // anonymous namespace


std::ostream& operator<<(std::ostream& os, ProgramInfo pi)
{
    //os
        //<< std::setfill('=') << std::setw(pi.width_)
        //<< "=" << std::endl
    //;

    outputLine(os, pi.projectNameLong, pi.width_);
    outputLine(os, pi.binaryName + " v" + pi.version, pi.width_);
    outputLine(os, pi.copyright, pi.width_);

    //os
        //<< std::setfill('=') << std::setw(pi.width_)
        //<< "=" << std::endl
    //;

    return os;
}


} // namespace common


} // namespace bcl2fastq


