/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file GzipDecompressor.hh
 *
 * \brief Declaration of GZip decompression filter.
 *
 * \author Roman Petrovski
 * \author Marek Balint
 * \author Mauricio Varea
 */


#ifndef BCL2FASTQ_IO_GZIPDECOMPRESSOR_HH
#define BCL2FASTQ_IO_GZIPDECOMPRESSOR_HH


#include <iosfwd>
#include <vector>

#include <cstddef>

#include <boost/iostreams/concepts.hpp>

#include <zlib.h>

#include "common/Exceptions.hh"
#include "common/StaticMemPool.hh"


namespace bcl2fastq {


namespace io {


struct z_stream_serialization : public z_stream
{
    z_stream_serialization(const z_stream &that) : z_stream(that){}
};


/// \brief Exception for ZLib errors.
class ZlibError : public common::Exception, public std::runtime_error
{
public:

    /// \brief Constructor.
    /// \param errorNumber Error number.
    /// \param zstream ZLib stream.
    ZlibError(int errorNumber, z_stream &zstream, const char *msg = 0);

    /// \brief Copy constructor.
    /// \param that Other exception object to copy from.
    ZlibError(const ZlibError &that);

public:

    /// \brief Get error message.
    /// \par Exception guarantee:
    /// no-throw
    virtual const char * what() const throw();
};


/// \brief Same as boost::gzip_decompressor, but ensures zlib buffer is allocated during construction.
class GzipDecompressor : public boost::iostreams::multichar_input_filter
{
public:

    /// \brief Memory pool page count.
    static const std::size_t PageCount = 3;

    /// \brief Memory pool page size.
    static const std::size_t PageSize = 65536;

public:

    /// \brief Intermediate buffer type definition.
    typedef std::vector<char> buffer_type;

public:

    /// \brief Constructor.
    /// \param bufferSize Intermediate buffer size.
    explicit GzipDecompressor(buffer_type::size_type bufferSize = PageSize);

    /// \brief Copy constructor.
    /// \param that Other instance of GZip decompressor to copy from.
    GzipDecompressor(const GzipDecompressor &that);

    /// \brief Copy assignment operator.
    /// \param rhs Right-hand-side parameter.
    GzipDecompressor & operator=(const GzipDecompressor &rhs);

    /// \brief Destructor.
    ~GzipDecompressor();

public:

    /// \brief Process closure notification.
    template<typename Source>
    void close(Source &);

    /// \brief Completely clears stream state. Resets input and output buffer information.
    void clear();

    /// \brief Reset the stream state so that new compressed block can be processed. Does not touch input or output positions.
    void reset();

    /// \brief Clear internal buffer.
    void flush();

public:

    /// \brief Read filtered characters.
    /// \param compressedStream Source to read from.
    /// \param targetBuffer Target buffer to read to.
    /// \param targetSize Maximum number of characters to be read.
    /// \return Number of characters read.
    template<typename Source>
    std::streamsize read(Source &compressedStream, char *targetBuffer, std::streamsize targetSize);

private:

    /// \brief Read filtered characters.
    /// \param compressedStream Source to read from.
    /// \param targetBuffer Target buffer to read to.
    /// \param targetSize Maximum number of characters to be read.
    /// \return Number of characters read.
    template<typename Source>
    std::streamsize readInternal(Source &compressedStream, char *targetBuffer, std::streamsize targetSize);

    /// \brief Initializes decompression stream.
    void init();

    /// \brief Process pending bytes from intermediate buffer.
    void processPendingBytes(bool endOfData);

private:

    /// \brief Memory pool type definition.
    typedef common::StaticMemPool<PageSize, PageCount> mempool_type;

private:

    /// \brief Meomory allocation routine for ZLib.
    /// \param opaque Private data object.
    /// \param items Number of items.
    /// \param size Size of an item.
    /// \return Pointer to allocated memory, or @c Z_NULL in case of failure.
    static voidpf zalloc OF((voidpf opaque, uInt items, uInt size));

    /// \brief Memory deallocation routine for ZLib.
    /// \param opaque Private data object.
    /// \param address Pointer to memory previously allocated by @c zalloc.
    static void zfree OF((voidpf opaque, voidpf address));

private:

    /// \brief ZLib stream.
    z_stream zstream_;

    /// \brief Memory pool for ZLib.
    mempool_type memPool_;

    /// \brief Number of bytes in intermediate buffer.
    std::streamsize pendingBytes_;

    /// \brief Intermediate buffer.
    buffer_type buffer_;
};


} // namespace io


} // namespace bcl2fastq


#include "io/GzipDecompressor.hpp"


#endif // BCL2FASTQ_IO_GZIPDECOMPRESSOR_HH


