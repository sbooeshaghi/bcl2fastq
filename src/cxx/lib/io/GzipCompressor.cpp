/**
* BCL to FASTQ file converter
* Copyright (c) 2007-2017 Illumina, Inc.
*
* This software is covered by the accompanying EULA
* and certain third party copyright/licenses, and any user of this
* source file is bound by the terms therein.
*
* \file GzipCompressor.cpp
*
* \brief Implementation of gzip compressor.
*
* \author Aaron Day
*/

#include "io/GzipCompressor.hh"

#include "common/Debug.hh"
#include "common/Exceptions.hh"

#include <algorithm>

namespace bcl2fastq
{
namespace io
{

GzipCompressor::GzipCompressor(std::vector<char>&       outputBuffer,
                               bool                     useBgzf,
                               const bios::gzip_params& gzipParams)
: outputBuffer_(outputBuffer)
, internalBuffer_()
, sink_(outputBuffer_)
, compressor_(gzipParams)
, uncompressed_in_(0)
, currentBlockBegin_(0)
, useBgzf_(useBgzf)
{
    internalBuffer_.reserve(bgzf_buffer_size_);
    outputBuffer_.clear();
}

GzipCompressor::~GzipCompressor()
{
    flush();
}

void GzipCompressor::prepareHeader()
{
    // save room for XFIELD
    outputBuffer_.insert(outputBuffer_.end(), sizeof(Header::XFIELD), 0);
}

std::streamsize GzipCompressor::write(const char* s, std::streamsize src_size)
{
    std::streamsize buffer_left(max_uncompressed_per_block_ - uncompressed_in_);

    // Jump through a hoop for the basespace uploader bug
    if (src_size > buffer_left)
    {
        buffer_left = 0;
    }

    const std::streamsize to_buffer(std::min<std::streamsize>(buffer_left, src_size));

    if (uncompressed_in_ == 0 && useBgzf_)
    {
        // Start of a new block, prepare the header.
        prepareHeader();
    }

    if (to_buffer)
    {
        internalBuffer_.insert(internalBuffer_.end(), s, s+to_buffer);
        uncompressed_in_ += to_buffer;
    }

    if (src_size != to_buffer)
    {
        flush();
        const std::streamsize more = src_size - to_buffer;
        const std::streamsize written = write(s + to_buffer, more);
        if (more != written)
        {
            return to_buffer + written;
        }
    }

    return src_size;
}

void GzipCompressor::flush()
{
    if (uncompressed_in_)
    {
        compressor_.write(sink_, internalBuffer_.data(), internalBuffer_.size());
        internalBuffer_.clear();
        compressor_.close(sink_, BOOST_IOS::out);
        if (useBgzf_)
        {
            rewriteHeader();
            currentBlockBegin_ = outputBuffer_.size();
        }
    }

    uncompressed_in_ = 0;
}

unsigned GzipCompressor::getBlockSize(char* buffer, bool failQuietly)
{
    Header *h = reinterpret_cast<Header*>(buffer);

    if (h->xfield.SI1 != static_cast<unsigned char>(66) || h->xfield.SI2 != static_cast<unsigned char>(67))
    {
        if (failQuietly)
            return 0;

        BCL2FASTQ_LOG(common::LogLevel::WARNING) << static_cast<unsigned short>(h->xfield.SI1) << ", " << static_cast<unsigned short>(h->xfield.SI2) << ", " << static_cast<unsigned short>(h->xfield.XLEN[0]) << ", " << static_cast<unsigned short>(h->xfield.XLEN[1]) << ", " << static_cast<unsigned short>(h->XFL) << ", " << static_cast<unsigned short>(h->OS) << std::endl;
        BOOST_THROW_EXCEPTION(common::IoError(ENOENT, "Corrupt GZIP block"));
    }

    return h->xfield.getBSIZE() + 1;
}

void GzipCompressor::rewriteHeader()
{
    memmove(&outputBuffer_[currentBlockBegin_], &outputBuffer_[currentBlockBegin_+sizeof(Header::XFIELD)], sizeof(Header) - sizeof(Header::XFIELD));
    Header *h(reinterpret_cast<Header*>(&outputBuffer_[currentBlockBegin_]));
    h->xfield = Header::XFIELD(outputBuffer_.size() - currentBlockBegin_ - 1); // why -1???
    h->FLG |= 0x04; // tell gzip that XLEN is in effect now.
}

GzipCompressor::Header::XFIELD::XFIELD(size_t compressedSize)
: SI1(66)
, SI2(67)
{
    XLEN[0] = sizeof(XFIELD) - sizeof(short);
    XLEN[1] = (sizeof(XFIELD) - sizeof(short)) >> 8;

    SLEN[0] = SUBFIELD_LENGTH;
    SLEN[1] = 0;

    static const size_t maxSize = ~((size_t)0x0) >> ((8-SUBFIELD_LENGTH)*8);
    BCL2FASTQ_ASSERT_MSG( compressedSize <= maxSize,
                          "Block size exceeds " << maxSize << " bytes when compressed: Got '" << compressedSize << "' bytes");

    static const size_t bitMask = 0x00000000000000FF;
    for (int i = 0; i < SUBFIELD_LENGTH; ++i)
    {
        BSIZE[i] = (compressedSize & (bitMask << (i*8))) >> (i*8);
    }
}


} // namespace io
} // namespace bcl2fastq

