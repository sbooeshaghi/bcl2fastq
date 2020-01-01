/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file BclFile.cpp
 *
 * \brief Implementation of BCL file.
 *
 * \author Marek Balint
 * \author Mauricio Varea
 * \author Aaron Day
 */


#include "io/GzipDecompressor.hh"
#include "data/BclFile.hh"
#include "common/Debug.hh"
#include "common/Logger.hh"
#include "common/SystemCompatibility.hh"

#include <boost/filesystem/operations.hpp>
#include <boost/format.hpp>

namespace bcl2fastq {
namespace data {


BclFile::BclFile(const common::RawDataBuffer& data,
                 io::GzipDecompressor&        decompressor,
                 bool                         ignoreErrors,
                 bool                         parseHeader,
                 bool                         resetDecompressor, /* = true */
                 common::ClustersCount        defaultClustersCount /*= 0*/)
: FileReaderBase(data,
                 ignoreErrors)
, BinaryFileReader(data,
                   ignoreErrors,
                   defaultClustersCount)
, gzipDecompressor_(&decompressor)
, isCompressed_(data.path_.extension() == ".gz" || data.path_.extension() == ".bgzf")
{
    if (resetDecompressor)
    {
        gzipDecompressor_->flush();
        gzipDecompressor_->reset();
    }

    static const std::string errMsgBegin = "Ignoring file open failure for bcl file '";

    try
    {
        bool headerRead = false;
        if (parseHeader)
        {
            headerRead = true;
            if (!readHeader(clustersCount_))
            {
                return;
            }

            /// \todo Refactoring: Do not determine endianness in runtime, let CMake to do it.
            if (!common::isLittleEndian())
            {
                clustersCount_ = common::swapEndian(clustersCount_);
            }
        }

        if (headerRead) {
            BCL2FASTQ_LOG(common::LogLevel::INFO)
                << "Opened BCL file '" << this->getPath().string() << "' with "
                << clustersCount_ << " clusters" << std::endl;
        }
    }
    CATCH_AND_IGNORE(std::ios_base::failure, errMsgBegin)
    CATCH_AND_IGNORE(io::ZlibError, errMsgBegin)
    CATCH_AND_IGNORE_ALL(errMsgBegin)
}

bool BclFile::getAndVerifyFileName(const boost::filesystem::path& inputDir,
                                   common::LaneNumber             laneNumber,
                                   common::CycleNumber            cycleNumber,
                                   bool                           ignoreErrors,
                                   boost::filesystem::path&       fileName)
{
    fileName = boost::filesystem::path(
        inputDir
        /
        boost::filesystem::path((boost::format("L%03d") % laneNumber).str())
        /
        boost::filesystem::path((boost::format("%04d.bcl.bgzf") % cycleNumber).str())
    );

    if (!boost::filesystem::exists(fileName))
    {
        std::string errStr("Unable to find BCL file: '" + fileName.string() + "'");
        if (ignoreErrors)
        {
            BCL2FASTQ_LOG(common::LogLevel::WARNING) << errStr << std::endl;
            return false;
        }
        else
        {
            BOOST_THROW_EXCEPTION(common::IoError(ENOENT, errStr));
        }

        fileName = boost::filesystem::path();
    }

    return true;
}

bool BclFile::getAndVerifyFileName(const boost::filesystem::path& inputDir,
                                   common::LaneNumber             laneNumber,
                                   common::TileNumber             tileNumber,
                                   common::CycleNumber            cycleNumber,
                                   bool                           ignoreErrors,
                                   boost::filesystem::path&       fileName)
{
    std::string sample = (boost::format("s_%d_%d") % laneNumber % tileNumber).str();

    fileName = boost::filesystem::path(
        inputDir
        /
        boost::filesystem::path((boost::format("L%03d") % laneNumber).str())
        /
        boost::filesystem::path((boost::format("C%d.1") % cycleNumber).str())
        /
        boost::filesystem::path((boost::format("%s.bcl.gz") % sample).str())
    );

    if (!boost::filesystem::exists(fileName))
    {
        fileName.replace_extension();

        if (!boost::filesystem::exists(fileName))
        {
            std::string errStr("Unable to find BCL file for '" + sample + "' in: " + fileName.remove_filename().string());
            if (ignoreErrors)
            {
                BCL2FASTQ_LOG(common::LogLevel::WARNING) << errStr << std::endl;
                return false;
            }
            else
            {
                BOOST_THROW_EXCEPTION(common::IoError(ENOENT, errStr));
            }

            fileName = boost::filesystem::path();
        }
    }

    return true;
}

std::streamsize BclFile::read(std::istream &is, char *buffer, std::streamsize size)
{
    if( isCompressed_ )
    {
        auto ret = gzipDecompressor_->read(is, buffer, size);
        return ret;
    }
    else
    {
        if (!is.read(buffer, size)) return -1;

        return is.gcount();
    }

    return 0;
}

std::streamsize BclFile::read(char* targetBuffer, std::streamsize targetSize)
{
    static const std::string errMsgBegin = "Ignoring read failure for bcl file '";

    if (!validateCondition(this->isOpen(), "File is not open")) { return 0; }

    std::streamsize bytesRead = -1;
    try
    {
        bytesRead = read(istr_, targetBuffer, targetSize);
        while ((bytesRead < targetSize) && !istr_.eof() && bytesRead != -1)
        {
            std::streamsize bytesRead2 = read(istr_, targetBuffer+bytesRead, targetSize-bytesRead);
            if (bytesRead2 != -1)
            {
                bytesRead += bytesRead2;
            } else {
                break;
            }
        }
    }
    CATCH_AND_IGNORE(io::ZlibError, errMsgBegin)
    CATCH_AND_IGNORE(common::IoError, errMsgBegin)
    CATCH_AND_IGNORE(std::ios_base::failure, errMsgBegin)
    CATCH_AND_IGNORE_ALL(errMsgBegin)

    return bytesRead;
}

bool BclFile::seek(std::streamsize compressedOffset, std::streamsize uncompressedOffset)
{
    if (!this->isOpen() && ignoreErrors_) {
        return false;
    }

    static const std::string errMsgBegin = "Ignoring seek failure for: '";

    bool success = false;
    try
    {
        istr_.seekg(compressedOffset);
        istr_.sync();
        std::vector<char> targetBuffer(uncompressedOffset);

        gzipDecompressor_->reset();
        gzipDecompressor_->flush();
        std::streamsize bytesRead = gzipDecompressor_->read(istr_, &*(targetBuffer.begin()), uncompressedOffset);

        int errnum = errno;
        success = (bytesRead == uncompressedOffset);
        if (!success)
        {
            BOOST_THROW_EXCEPTION(common::InputDataError(errnum, (boost::format("Unable to seek to given position of BCL file '%s': compressed_offset=%d uncompressed_offset=%d bytes_read=%d") % this->getPath().string() % compressedOffset % uncompressedOffset % bytesRead).str()));
        }
    }
    CATCH_AND_IGNORE(common::IoError, errMsgBegin)
    CATCH_AND_IGNORE_ALL(errMsgBegin)

    return success;
}


} // namespace data
} // namespace bcl2fastq


