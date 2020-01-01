/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file FastqFile.hh
 *
 * \brief Declaration of FASTQ file.
 *
 * \author Marek Balint
 */


#ifndef BCL2FASTQ_DATA_FASTQFILE_HH
#define BCL2FASTQ_DATA_FASTQFILE_HH


#include <boost/shared_ptr.hpp>
#include <boost/filesystem/path.hpp>

#include "io/FileBufWithReopen.hh"


namespace bcl2fastq {


namespace data {


class FastqFileState;


/// \brief FASTQ file.
class FastqFile
{
public:

    /// \brief FASTQ file state smart pointer type definition.
    typedef boost::shared_ptr<FastqFileState> FastqFileStatePtr;

public:

    /// \brief Open FASTQ file.
    /// \param path Path to the FASTQ file to be opened.
    /// \param mode Access mode to be used.
    /// \return External state of opened FASTQ file.
    static FastqFile::FastqFileStatePtr openFile(const boost::filesystem::path &path,
                                                 std::ios_base::openmode        mode);

public:

    /// \brief Get file path.
    /// \param fileState FASTQ file external state.
    /// \return File path.
    static const boost::filesystem::path& getPath(const FastqFileState &fileState);

    /// \brief Write bytes from buffer to file.
    /// \param fileState FASTQ file external state.
    /// \param sourceBuffer Source buffer to read to.
    /// \param sourceSize Number of bytes to be written.
    static void write(FastqFileState &fileState, const char *sourceBuffer, std::streamsize sourceSize);
};


/// \brief FASTQ file internal state.
class FastqFileState
{
protected:

    /// \brief Default constructor.
    FastqFileState();

    /// \brief Virtual destructor.
    virtual ~FastqFileState() = 0;

protected:

    /// \brief Write bytes from buffer to file.
    /// \param os Output stream to write to.
    /// \param sourceBuffer Source buffer to read to.
    /// \param sourceSize Number of bytes to be written.
    virtual void write(std::ostream &os, const char *sourceBuffer, std::streamsize sourceSize) = 0;

private:

    /// \brief File path.
    boost::filesystem::path path_;

    /// \brief Filebuf.
    io::FileBufWithReopen fileBuf_;

private:

    friend class FastqFile;
};


} // namespace data


} // namespace bcl2fastq


#endif // BCL2FASTQ_DATA_FASTQFILE_HH


