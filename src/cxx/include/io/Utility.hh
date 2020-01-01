/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file Utility.hh
 *
 * \brief Declaration of various I/O related utilities.
 *
 * \author Roman Petorvski
 * \author Marek Balint
 */


#ifndef BCL2FASTQ_IO_UTILITY_HH
#define BCL2FASTQ_IO_UTILITY_HH


#include <ios>


namespace bcl2fastq {


namespace io {


/// \brief Translate ios access flags to stdio mode.
/// \param mode Access flags used by ios.
/// \return Mode string used by stdio functions.
const char * iosFlagsToStdioMode(std::ios_base::openmode mode);


} // namespace io


} // namespace bcl2fastq


#endif // BCL2FASTQ_IO_UTILITY_HH


