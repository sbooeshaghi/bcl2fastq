/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file GzipDecompressor.cpp
 *
 * \brief Implementation of GZip decompression filter.
 *
 * \author Marek Balint
 * \author Mauricio Varea
 */


#include <boost/iostreams/read.hpp>
#include <boost/lexical_cast.hpp>

#include "common/Debug.hh"
#include "common/Logger.hh"
#include "io/GzipDecompressor.hh"


namespace bcl2fastq {


namespace io {


inline std::ostream &operator <<(std::ostream &os, const z_stream_serialization &zs)
{
    return os << "z_stream( " <<
        " next_in: " << (void*)zs.next_in <<
        ", avail_in: " << zs.avail_in <<
        ", total_in: " << zs.total_in <<
        ", next_out: " << (void*)zs.next_out <<
        ", avail_out: " << zs.avail_out <<
        ", total_out: " << zs.total_out <<
        ", msg: '" << (zs.msg ? (const char*)zs.msg : (const char*)("null"))<<
        "', state: " << (void*)zs.state <<
        ", zalloc: " << zs.zalloc <<
        ", zfree: " << zs.zfree <<
        ", opaque: " << zs.opaque <<
        ", data_type: " << zs.data_type <<
        ", adler: " << zs.adler <<
        ", reserved: " << zs.reserved <<
        " )";
}


ZlibError::ZlibError(int errorNumber, z_stream &zstream, const char *msg)
: Exception(errorNumber, (zstream.msg ? zstream.msg : msg ?
        std::string(msg) + " unknown error " + boost::lexical_cast<std::string>(errorNumber) :
        "unknown error " + boost::lexical_cast<std::string>(errorNumber)) + std::string("\n") +
        boost::lexical_cast<std::string>(z_stream_serialization(zstream)))
, std::runtime_error((zstream.msg ? zstream.msg : msg ?
        std::string(msg) + " unknown error " + boost::lexical_cast<std::string>(errorNumber) :
        "unknown error " + boost::lexical_cast<std::string>(errorNumber)) + std::string("\n") +
        boost::lexical_cast<std::string>(z_stream_serialization(zstream)))
{
}

ZlibError::ZlibError(const ZlibError &that)
: Exception(that)
, std::runtime_error(that)
{
}

const char * ZlibError::what() const throw()
{
    return this->getMessage().c_str();
}


GzipDecompressor::GzipDecompressor(buffer_type::size_type bufferSize)
: zstream_()
, memPool_()
, pendingBytes_(0)
, buffer_(bufferSize, '\0')
{
    this->clear();
}

GzipDecompressor::GzipDecompressor(const GzipDecompressor &that)
: zstream_()
, memPool_()
, pendingBytes_(that.pendingBytes_)
, buffer_(that.buffer_)
{
//    zstream_.state = nullptr;
    this->reset();
}

GzipDecompressor & GzipDecompressor::operator=(const GzipDecompressor &rhs)
{
    if (this != &rhs)
    {
        this->reset();
        pendingBytes_ = rhs.pendingBytes_;
        buffer_ = rhs.buffer_;
    }
    return *this;
}

GzipDecompressor::~GzipDecompressor()
{
    if (zstream_.state)
    {
        inflateEnd(&zstream_);
    }
}

void GzipDecompressor::init()
{
    zstream_.zalloc = zalloc;
    zstream_.zfree = zfree;
    zstream_.opaque = &memPool_;

    int ret = inflateInit2(&zstream_, 16+MAX_WBITS);
    if (ret != Z_OK)
    {
        BOOST_THROW_EXCEPTION(ZlibError(ret, zstream_));
    }
}


void GzipDecompressor::clear()
{
    if (zstream_.state)
    {
        inflateEnd(&zstream_);
    }
    memset(&zstream_, 0, sizeof(zstream_));
    this->init();
}

void GzipDecompressor::reset()
{
    if (zstream_.state)
    {
        inflateEnd(&zstream_);
    }
    this->init();
}

void GzipDecompressor::flush()
{
    zstream_.avail_in = 0;
    zstream_.next_in = Z_NULL;
    pendingBytes_ = 0;
}

void GzipDecompressor::processPendingBytes(bool endOfData)
{
    if (pendingBytes_ <= 0)
    {
        return;
    }
    zstream_.next_in = reinterpret_cast<Bytef *>(&buffer_.begin()[0]);
    zstream_.avail_in = pendingBytes_;

    BCL2FASTQ_ASSERT_MSG( PageSize >= zstream_.avail_in,
             "There are more bytes to process than a single page can take" );
    int err = inflate(&zstream_, Z_SYNC_FLUSH);
    if (err != Z_OK && err != Z_STREAM_END)
    {
        if (endOfData)
        {
            BOOST_THROW_EXCEPTION(ZlibError(err, zstream_, "Premature end of compressed stream reached. "));
        } else {
            BOOST_THROW_EXCEPTION(ZlibError(err, zstream_));
        }
    }
    const std::streamsize compressedWritten = pendingBytes_ - zstream_.avail_in;
    std::copy(
        buffer_.begin() + compressedWritten,
        buffer_.begin() + pendingBytes_,
        buffer_.begin()
    );
    pendingBytes_ = zstream_.avail_in;
}

voidpf GzipDecompressor::zalloc OF((voidpf opaque, uInt items, uInt size))
{
    voidpf ptr = static_cast<voidpf>(
        reinterpret_cast<mempool_type *>(opaque)->allocate(items*size)
    );

    if (ptr == NULL)
    {
        return Z_NULL;
    }
    return ptr;
}

void GzipDecompressor::zfree OF((voidpf opaque, voidpf address))
{
    reinterpret_cast<mempool_type *>(opaque)->deallocate(address);
}


} // namespace io


} // namespace bcl2fastq


