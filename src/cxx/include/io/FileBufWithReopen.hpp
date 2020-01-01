/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file FileBufWithReopen.hpp
 *
 * \brief Implementation of the file buffer with reopen.
 *
 * \author Roman Petorvski
 * \author Marek Balint
 */


#ifndef BCL2FASTQ_IO_FILEBUFWITHREOPEN_HPP
#define BCL2FASTQ_IO_FILEBUFWITHREOPEN_HPP


#include <iostream>
#include <fcntl.h>

#include "common/Exceptions.hh"
#include "common/Debug.hh"
#include "common/Logger.hh"
#include "io/Utility.hh"
#include "io/FileBufWithReopen.hh"


namespace bcl2fastq {


namespace io {


template<typename CharT, typename Traits>
BasicFileBufWithReopen<CharT, Traits>::BasicFileBufWithReopen(std::ios_base::openmode mode)
: std::basic_filebuf<CharT,Traits>()
, mode_(mode)
#ifndef WIN32
, file_(std::basic_filebuf<CharT, Traits>::_M_file) // implementation dependent
#endif
{
    if (!this->reserve())
    {
        BOOST_THROW_EXCEPTION(common::IoError(errno, "Failed to allocate a file handle"));
    }
}

template<typename CharT, typename Traits>
BasicFileBufWithReopen<CharT, Traits>* BasicFileBufWithReopen<CharT, Traits>::reopen(const boost::filesystem::path& path, typename FadviseFlags::value_type fadvise)
{
    const char * const openMode = iosFlagsToStdioMode(mode_);
    BCL2FASTQ_ASSERT_MSG(openMode, "Combination of open mode flags is invalid");
    //BCL2FASTQ_LOG(common::LogLevel::DEBUG) << "Reopening '" << path << "' with mode '" << openMode << "'" << std::endl;

#ifndef WIN32
    // I haven't figured out a way to reopen a file buffer on windows...
    BCL2FASTQ_ASSERT_MSG(this->is_open(), "The file must be already open for reopen to be possible");
#endif

    if ((mode_ & std::ios_base::out))
    {
        //flush any pending data before reopening a (usually) different file, reset eof flag
        int syncResult = this->pubsync();
        BCL2FASTQ_ASSERT_MSG(syncResult == 0, "pubsync failed: " << errno << ": " << strerror(errno));
    }
#ifndef WIN32
    /// \todo Refactoring: Spaghetti design: Move repeating code that calls posix_fadvise to separate function.

    if (
        (fadvise & FadviseFlags::NOREUSE)
        &&
        ::posix_fadvise(::fileno(file_.getCFile()), 0, 0, POSIX_FADV_NOREUSE)
        &&
        errno // reqired since POSIX_FADV_NOREUSE fails with errno 0 on /dev/null
    ) {
        BCL2FASTQ_LOG(common::LogLevel::WARNING)
            << "posix_fadvise failed for POSIX_FADV_NOREUSE: "
            << errno << ": " << strerror(errno) << ": "
            << path << std::endl
        ;
    }

    if (
        (fadvise & FadviseFlags::WILLNEED)
        &&
        posix_fadvise(::fileno(file_.getCFile()), 0, 0, POSIX_FADV_WILLNEED)
        &&
        errno // reqired since POSIX_FADV_WILLNEED fails with errno 0 on /dev/null 
    ) {
        BCL2FASTQ_LOG(common::LogLevel::WARNING)
            << "posix_fadvise failed for POSIX_FADV_WILLNEED: "
            << errno << ": " << strerror(errno) << ": "
            << path << std::endl
        ;
    }

    if (
        (fadvise & FadviseFlags::DONTNEED)
        &&
        posix_fadvise(::fileno(file_.getCFile()), 0, 0, POSIX_FADV_DONTNEED)
        &&
        errno // reqired since POSIX_FADV_DONTNEED fails with errno 0 on /dev/null 
    ) {
        BCL2FASTQ_LOG(common::LogLevel::WARNING)
            << "posix_fadvise failed for POSIX_FADV_DONTNEED: "
            << errno << ": " << strerror(errno) << ": "
            << path << std::endl
        ;
    }

    errno = 0;
    std::FILE *result = std::freopen(path.string().c_str(), openMode, file_.getCFile());
    int errnum = errno;
    if (result == NULL)
    {
        this->close();
        errno = errnum;
        return NULL;
    }
    BCL2FASTQ_ASSERT_MSG(file_.getCFile() == result, "According to specs, returned pointer must be the same as the one passed to freopen");
#else
    errno = 0;
    if (this->is_open())
    {
        close();
    }
    open(path.string().c_str(), mode_);
#endif

    if (!(mode_ & std::ios_base::app))
    {
        if (this->seekoff(0, std::ios_base::beg, mode_) != std::streampos(0))
        {
            BOOST_THROW_EXCEPTION(common::IoError(errno, path.string()));
        }
    }
#ifndef WIN32
    if ((fadvise & FadviseFlags::SEQUENTIAL) && posix_fadvise(::fileno(result), 0, 0, POSIX_FADV_SEQUENTIAL))
    {
        BOOST_THROW_EXCEPTION(common::IoError(errno, path.string()));
    }

    if ((fadvise & FadviseFlags::RANDOM) && posix_fadvise(::fileno(result), 0, 0, POSIX_FADV_RANDOM))
    {
        BOOST_THROW_EXCEPTION(common::IoError(errno, path.string()));
    }
#endif
    return this;
}

template<typename CharT, typename Traits>
std::ios_base::openmode BasicFileBufWithReopen<CharT, Traits>::getMode() const
{
    return mode_;
}

template<typename CharT, typename Traits>
void BasicFileBufWithReopen<CharT, Traits>::setMode(std::ios_base::openmode mode)
{
    this->reset();
    mode_ = mode;
}

template<typename CharT, typename Traits>
void BasicFileBufWithReopen<CharT, Traits>::flush()
{
    /// \todo Refactoring: Make some sensible error reporting.

    //flush any pending data before reopening a (usually) different file
    int result = this->pubsync();
    BCL2FASTQ_ASSERT_MSG(result == 0, "TODO: make some sensible error reporting...");
}


namespace detail {


/// \brief Get file that will exist on all systems.
/// \return String containing path to a file which always exists.
/// \todo Portability: Make this function Windows-compatible.
static const char * fileThatAlwaysExists()
{
#ifdef WIN32
    return "NUL";
#else
    return "/dev/null";
#endif
}


} // namespace detail


template<typename CharT, typename Traits>
void BasicFileBufWithReopen<CharT, Traits>::reset()
{
    if (this->is_open())
    {
        this->reopen(detail::fileThatAlwaysExists());
    }
    else
    {
        this->reserve();
    }
    BCL2FASTQ_ASSERT_MSG(this->is_open(), "Failed to reset file buffer");
}

template<typename CharT, typename Traits>
bool BasicFileBufWithReopen<CharT, Traits>::reserve()
{
    return !!this->open(detail::fileThatAlwaysExists(), mode_);
}


#ifndef WIN32
template<typename CharT, typename Traits>
BasicFileBufWithReopen<CharT, Traits>::File::File(file_type &file)
: file_(file)
{
}

template<typename CharT, typename Traits>
std::FILE * BasicFileBufWithReopen<CharT, Traits>::File::getCFile()
{
    return file_.file();
}
#endif

template<typename CharT>
std::streamsize read(BasicFileBufWithReopen<CharT> &fileBuf, char *targetBuffer, std::streamsize targetSize)
{
    std::basic_istream<CharT> is(&fileBuf);

    is.read(targetBuffer, targetSize);
    if (is)
    {
        return targetSize;
    }
    else
    {
        return is.gcount();
    }
}


} // namespace io


} // namespace bcl2fastq


#endif // BCL2FASTQ_IO_FILEBUFWITHREOPEN_HPP


