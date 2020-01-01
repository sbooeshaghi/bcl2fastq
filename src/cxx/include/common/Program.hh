/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file Program.hh
 *
 * \brief Declaration of the skeleton of all c++ programs.
 *
 * \author Come Raczy
 * \author Roman Petrovski
 * \author Marek Balint
 */


#ifndef BCL2FASTQ_COMMON_PROGRAM_HH
#define BCL2FASTQ_COMMON_PROGRAM_HH


#include <string>

#include "common/ProgramInfo.hh"


namespace bcl2fastq {


namespace common {


/// \brief Unified behavior of all programs.
template<typename O>
void run(int(*callback)(const ProgramInfo &, O &), int argc, const char *argv[]);


} // namespace common


} // namespace bcl2fastq


#include "common/Program.hpp"


#endif // BCL2FASTQ_COMMON_PROGRAM_HH


