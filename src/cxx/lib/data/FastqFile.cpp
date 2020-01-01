/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file FastqFile.cpp
 *
 * \brief Implementation of FASTQ file.
 *
 * \author Marek Balint
 */


#include <boost/format.hpp>

#include "common/Debug.hh"
#include "common/Logger.hh"
#include "common/SystemCompatibility.hh"
#include "data/FastqFile.hh"


namespace bcl2fastq {


namespace data {


FastqFileState::FastqFileState()
: path_()
, fileBuf_(std::ios_base::binary | std::ios_base::out | std::ios_base::app)
{
}

FastqFileState::~FastqFileState()
{
}

void FastqFileState::write(std::ostream &os, const char *sourceBuffer, std::streamsize sourceSize)
{
    os.write(sourceBuffer, sourceSize);
}


namespace detail {


/// \brief Default FASTQ file internal state.
class DefaultFastqFileState : public FastqFileState
{
private:

    virtual void write(std::ostream &os, const char *sourceBuffer, std::streamsize sourceSize);
};

void DefaultFastqFileState::write(std::ostream &os, const char *sourceBuffer, std::streamsize sourceSize)
{
    FastqFileState::write(os, sourceBuffer, sourceSize);
}


} // namespace detail


FastqFile::FastqFileStatePtr FastqFile::openFile(const boost::filesystem::path &path,
                                                 std::ios_base::openmode        mode)
{
    FastqFileStatePtr fileState(new detail::DefaultFastqFileState);
    fileState->path_ = path;

    io::FileBufWithReopen fileBuf(std::ios_base::out | mode);
    fileBuf.reopen(path, io::FileBufWithReopen::FadviseFlags::SEQUENTIAL_ONCE);
    int errnum = errno;
    if (!fileBuf.is_open())
    {
        BOOST_THROW_EXCEPTION(common::IoError(errnum, (boost::format("Unable to open FASTQ file '%s' for writing") % path.string()).str()));
    }
    fileState->fileBuf_.reopen(path, io::FileBufWithReopen::FadviseFlags::SEQUENTIAL_ONCE);
    errnum = errno;
    if (!fileState->fileBuf_.is_open())
    {
        BOOST_THROW_EXCEPTION(common::IoError(errnum, (boost::format("Unable to open FASTQ file '%s' for writing") % path.string()).str()));
    }
    BCL2FASTQ_LOG(common::LogLevel::INFO) << "Created FASTQ file '" << path << std::endl;

    return fileState;
}

const boost::filesystem::path& FastqFile::getPath(const FastqFileState &fileState)
{
    return fileState.path_;
}

void FastqFile::write(FastqFileState &fileState, const char *sourceBuffer, std::streamsize sourceSize)
{
    //fileBuf_.reopen(fileState.path_.native().c_str(), io::FileBufWithReopen::FadviseFlags::SEQUENTIAL_ONCE);
    //BCL2FASTQ_ASSERT_MSG(fileBuf_.is_open(), "Unable to reopen FASTQ file '" << fileState.path_.native() << "' for append");
    BCL2FASTQ_ASSERT_MSG(fileState.fileBuf_.is_open(), "FASTQ file '" << fileState.path_ << "' is not opened");
    std::ostream os(&fileState.fileBuf_);

    fileState.write(os, sourceBuffer, sourceSize);
}


} // namespace data


} // namespace bcl2fastq


