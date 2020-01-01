/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file ProgramInfo.hh
 *
 * \brief Basic program information.
 *
 * \author Mauricio Varea
 * \author Marek Balint
 */

#ifndef BCL2FASTQ_COMMON_PROGRAMINFO_HH
#define BCL2FASTQ_COMMON_PROGRAMINFO_HH


#include <string>
#include <iosfwd>

#include <boost/filesystem.hpp>


namespace bcl2fastq {


namespace common {


/// \brief Program information.
class ProgramInfo
{

public:

    /// \brief Constructor.
    /// \param cmd Command used to invoke the program (argv[0]).
    /// \param width Line width (used for output).
    ProgramInfo(const std::string &cmd, unsigned int width);

public:

    /// \brief Command used to invoke the program (argv[0]).
    const boost::filesystem::path command;

    /// \brief Name of binary executable.
    const std::string binaryName;

    /// \brief Short project name (one word).
    const std::string projectNameShort;

    /// \brief Long project name (descriptive sentence).
    const std::string projectNameLong;

    /// \brief Copyright string.
    const std::string copyright;

    /// \brief Version string (NN.YY.MM.DD)
    const std::string version;

private:

    /// \brief Line width (used for output).
    const unsigned int width_;

    friend std::ostream& operator<<(std::ostream& os, ProgramInfo pi);
};


/// \brief Output operator for ProgramInfo.
/// \param os Output stream.
/// \param pi Program info to be output.
/// \return Output stream.
std::ostream& operator<<(std::ostream& os, ProgramInfo pi);


} // namespace common


} // namespace bcl2fastq


#endif // BCL2FASTQ_COMMON_PROGRAMINFO_HH


