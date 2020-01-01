/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file SyncFile.cpp
 *
 * \brief Implementation of Synchronized file.
 *
 * \author Aaron Day
 */

#include "io/SyncFile.hh"

namespace bcl2fastq {
namespace io {

FileReaderBase::FileReaderBase()
: path_()
, ignoreErrors_(false)
{
}

FileReaderBase::FileReaderBase(const boost::filesystem::path& path,
                               bool                           ignoreErrors)
: path_(path)
, ignoreErrors_(ignoreErrors)
{
}

UnprocessedFile::UnprocessedFile(std::ios_base::openmode openFlags)
: FileReaderBase()
, fileBuf_(openFlags)
, fileSize_(0)
{
}

UnprocessedFile::UnprocessedFile(const boost::filesystem::path& path,
                                 bool                           ignoreErrors,
                                 std::ios_base::openmode        openFlags /*=std::ios_base::in | std::ios_base::binary*/)
: FileReaderBase(path, ignoreErrors)
, fileBuf_(openFlags)
, fileSize_(0)
{
    openFileBuf();
}

UnprocessedFile::~UnprocessedFile()
{
}

void UnprocessedFile::openFile(const boost::filesystem::path& path,
                               bool                           ignoreErrors,
                               size_t)
{
    openFileImpl(path, ignoreErrors);
}

void UnprocessedFile::openFileImpl(const boost::filesystem::path& path,
                                   bool                           ignoreErrors)
{
    path_ = path;
    ignoreErrors_ = ignoreErrors;

    openFileBuf();
}

bool UnprocessedFile::openFileBuf()
{
    try
    {
        errno = 0;

        fileBuf_.reopen(getPath(), io::FileBufWithReopen::FadviseFlags::SEQUENTIAL_ONCE);
        int errnum = errno;
        if (!fileBuf_.is_open())
        {
            fileBuf_.reset();
            boost::format errFormat(boost::format("Unable to open %s file '%s' for reading")
                                         % getFileTypeStr()
                                         % getPath().string());
            if (ignoreErrors_)
            {
                BCL2FASTQ_LOG(common::LogLevel::WARNING)
                    << errFormat
                    << ":" << std::strerror(errnum) << " (" << errnum << ")" << std::endl;

                return false;
            }
            else
            {
                BOOST_THROW_EXCEPTION(common::IoError(errnum, errFormat.str()));
            }
        }
        else
        {
            fileSize_ = fileBuf_.pubseekoff(0, std::ios_base::end);
            fileBuf_.pubseekoff(0, std::ios_base::beg);
        }
    }
    CATCH_AND_IGNORE(common::IoError, "Ignoring read failure for " + getFileTypeStr() + " file '")
    CATCH_AND_IGNORE(std::ios_base::failure, "Ignoring read failure for " + getFileTypeStr() + " file '")
    CATCH_AND_IGNORE_ALL("Ignoring read failure for " + getFileTypeStr() + " file '")

    return true;
}

bool UnprocessedFile::readBytes(std::vector<char>& buffer,
                                uint32_t bytes)
{
    BCL2FASTQ_ASSERT_MSG(buffer.size() == bytes, "Incorrect buffer size");

    auto bytesRead = io::read(fileBuf_, buffer.data(), bytes);
    if (bytesRead != bytes)
    {
        logError(bytesRead,
                 bytes);

        buffer.resize(bytesRead > 0 ? bytesRead : 0);

        return false;
    }

    return true;
}

bool UnprocessedFile::readEntireFile(common::RawDataBuffer& buffer)
{
    bool ret = false;

    buffer.path_ = path_;

    auto filePos = fileBuf_.pubseekoff(0, std::ios_base::cur);

    common::resizeAndRealloc(buffer, fileSize_ - filePos);
    ret = readBytes(buffer, buffer.size());
    if (!ret)
    {
        BCL2FASTQ_LOG(common::LogLevel::WARNING) << "Returning false!!!" << std::endl;
    }

    return ret;
}

void UnprocessedFile::logError(std::streamsize bytesRead,
                               std::streamsize bytesExpected) const
{
    int errnum = errno;

    boost::format errFormat(boost::format(
        getFileTypeStr() + " file '%s' truncated: bytes_read=%d bytes_expected=%d") %
        this->getPath().string() % bytesRead % bytesExpected);

    if (ignoreErrors_)
    {
        BCL2FASTQ_LOG(common::LogLevel::WARNING)
            << errFormat
            << ":" << std::strerror(errnum) << " (" << errnum << ")" << std::endl;
    }
    else
    {
        BOOST_THROW_EXCEPTION(common::InputDataError(errnum,
                                                     errFormat.str()));
    }
}

SyncFile::SyncFile(std::ios_base::openmode openFlags)
: UnprocessedFile(openFlags)
, mut_()
, cv_()
, currentTile_(0)
, buffer_()
{
}

SyncFile::SyncFile(const boost::filesystem::path& path,
                   bool                           ignoreErrors,
                   std::ios_base::openmode        openFlags)
: UnprocessedFile(path,
                  ignoreErrors,
                  openFlags)
, mut_()
, cv_()
, currentTile_(0)
, buffer_()
{
}

SyncFile::~SyncFile()
{
}

bool SyncFile::read(char* buffer,
                    std::streamsize bytes)
{
    auto bytesRead = io::read(fileBuf_, buffer, bytes);
    if (bytesRead != bytes)
    {
        logError(bytesRead,
                 bytes);

        return false;
    }

    return true;
}

bool SyncFile::read(std::vector<char>& buffer,
                    std::streamsize bytes)
{
    common::resizeAndRealloc(buffer, bytes);
    return readBytes(buffer, bytes);
}

bool SyncFile::seek(std::streamsize offset,
                    std::ios_base::seekdir way /*= std::ios_base::cur*/)
{
    return fileBuf_.pubseekoff(offset, way);
}

bool SyncFile::appendBgzfBlock(std::vector<char> &gzipData)
{
    static const unsigned headerSize = GzipCompressor::getHeaderSize();
    const size_t dataOffset = gzipData.size();

    std::streamsize currentPos = fileBuf_.pubseekoff(0, std::ios_base::cur);
    if (1 + currentPos + headerSize >= static_cast<int64_t>(fileSize_)) {
        return true;
    }

    unsigned bytesToRead = headerSize;
    common::resizeAndRealloc(gzipData, gzipData.size()+bytesToRead);

    auto bytesRead = io::read(fileBuf_, gzipData.data()+dataOffset, bytesToRead);
    if (bytesRead != bytesToRead) {
        common::resizeAndRealloc(gzipData, dataOffset);
        return false;
    }

    unsigned gzipBlockSize = 0;
    try {
        gzipBlockSize = io::GzipCompressor::getBlockSize(&gzipData[dataOffset], true);
    }
    catch (...) {
        common::resizeAndRealloc(gzipData, dataOffset);
        return false;
    }

    bytesToRead = gzipBlockSize - bytesToRead;
    common::resizeAndRealloc(gzipData, gzipData.size()+bytesToRead);

    bytesRead = io::read(fileBuf_, gzipData.data()+dataOffset+headerSize, bytesToRead);
    if (bytesRead != bytesToRead) {
        common::resizeAndRealloc(gzipData, dataOffset);
        return false;
    }

    return true;
}

bool SyncFile::readGzipBlocks(std::vector<char>& gzipData,
                              std::streamsize startPos,
                              std::streamsize nextStartPos)
{
    BCL2FASTQ_ASSERT_MSG(nextStartPos == 0 || startPos <= nextStartPos, "Invalid aggregated read block");

    fileBuf_.pubseekoff(startPos, std::ios_base::beg);

    gzipData.clear();
    const bool readToEOF = (nextStartPos == 0);
    nextStartPos = (nextStartPos <= 0) ? fileSize_ : nextStartPos;
    std::streamsize bytes = nextStartPos - startPos;

    try {
        common::resizeAndRealloc(gzipData, bytes);
    }
    catch (std::bad_alloc e) {
        BCL2FASTQ_LOG(common::LogLevel::FATAL) << "Current tile allocation exceed memory available. " << bytes << std::endl;
        return false;
    }

    if (bytes > 0) {

        auto bytesRead = io::read(fileBuf_, gzipData.data(), bytes);

        if (bytesRead != bytes)
        {
            logError(bytesRead, bytes);

            gzipData.resize(bytesRead > 0 ? bytesRead : 0);

            return false;
        }

        // Sanity check: gzipData must start with a gzip header.
        if (gzipData.size() > 0)
        {
            try
            {
                /*unsigned gzipBlockSize =*/ io::GzipCompressor::getBlockSize(&gzipData[0]);
            }
            catch (common::IoError e)
            {
                if (ignoreErrors_)
                {
                    BCL2FASTQ_LOG(common::LogLevel::WARNING) << "Identified corrupt GZIP data - " << getPath() << std::endl;
                    return false;
                }

                throw;
            }
        }
    }

    // Read one more gzip block
    if (!readToEOF) {
        appendBgzfBlock(gzipData);
    }

    return true;
}

} // namespace io
} // namespace bcl2fastq

