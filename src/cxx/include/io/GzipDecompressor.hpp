/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file GzipDecompressor.hpp
 *
 * \brief Implementation of GZip decompression filter.
 *
 * \author Marek Balint
 * \author Mauricio Varea
 */


#ifndef BCL2FASTQ_IO_GZIPDECOMPRESSOR_HPP
#define BCL2FASTQ_IO_GZIPDECOMPRESSOR_HPP


#include <boost/iostreams/read.hpp>
#include <boost/foreach.hpp>

#include "common/Debug.hh"
#include "common/Logger.hh"
#include "io/GzipDecompressor.hh"


namespace bcl2fastq {


namespace io {


template<typename Source>
void GzipDecompressor::close(Source &)
{
    this->reset();
}

template<typename Source>
std::streamsize GzipDecompressor::read(Source &compressedStream, char *targetBuffer, std::streamsize targetSize)
{
    std::streamsize bytesRead = readInternal(compressedStream, targetBuffer, targetSize);
    while ((bytesRead < targetSize) && !compressedStream.eof() && bytesRead != -1)
    {
        std::streamsize bytesRead2 = readInternal(compressedStream, targetBuffer+bytesRead, targetSize-bytesRead);
        if (bytesRead2 != -1)
        {
            bytesRead += bytesRead2;
        } else {
            break;
        }
    }

    return bytesRead;
}

template<typename Source>
std::streamsize GzipDecompressor::readInternal(Source &compressedStream, char *targetBuffer, std::streamsize targetSize)
{
    if (targetSize == 0)
    {
        return 0;
    }

    zstream_.next_out = reinterpret_cast<Bytef *>(targetBuffer);
    zstream_.avail_out = static_cast<uInt>(targetSize);
    this->processPendingBytes(!compressedStream.good());

    if (zstream_.avail_out)
    {
        if (compressedStream.good())
        {
            std::streamsize bytesToRead = (buffer_.end() - buffer_.begin()) - pendingBytes_;
            std::streamsize result = 0;
            if (bytesToRead > 0) 
            {
                result = boost::iostreams::read(compressedStream, &buffer_.begin()[pendingBytes_], bytesToRead);
            }

            if (result != -1 /* EOF */)
            {
                pendingBytes_ += result;
            }
        }
        if (!compressedStream.good() && !compressedStream.eof())
        {
            BOOST_THROW_EXCEPTION(common::IoError(errno, "Failed to read compressed data"));
        }

        this->processPendingBytes(!compressedStream.good());
    }

    std::streamsize ret = targetSize - zstream_.avail_out;

    if (ret == 0)
    {
        if (pendingBytes_ > 0)
        {
            this->reset();
            this->processPendingBytes(!compressedStream.good());
            ret = targetSize - zstream_.avail_out;
        } else if (!compressedStream.good()) {
            BOOST_THROW_EXCEPTION(common::IoError(errno, "Unexpected end of file"));
        }

        //BCL2FASTQ_ASSERT_MSG( pendingBytes_ == 0,
        //    "When no bytes come out of decompressor and the input stream is all finished,"
        //    " expecting the pendingBytes_ to be 0. Actual: " << pendingBytes_ );

        //BCL2FASTQ_LOG(common::LogLevel::TRACE) << "GzipDecompressor::read: finished" << std::endl;
        //return -1;
    }

    return ret ? ret : -1;
}


} // namespace io


} // namespace bcl2fastq


#endif // BCL2FASTQ_IO_GZIPDECOMPRESSOR_HPP


