/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file GzipCompressor.hh
 *
 * \brief Declaration of Gzip Compressor.
 *
 * \author Roman Petrovski
 * \author Aaron Day
 */

#ifndef BCL2FASTQ_IO_GZIP_COMPRESSOR_HH
#define BCL2FASTQ_IO_GZIP_COMPRESSOR_HH

#include <vector>
#include <boost/noncopyable.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

namespace bcl2fastq
{
namespace io
{

namespace bios=boost::iostreams;

class GzipCompressor : private boost::noncopyable
{
public:
    GzipCompressor(std::vector<char>&       outputBuffer,
                   bool                     useBgzf,
                   const bios::gzip_params& gzipParams = bios::gzip::default_compression);

    ~GzipCompressor();

    std::streamsize write(const char* s, std::streamsize n);

    void flush();

    static const int SUBFIELD_LENGTH = 2;  // this needs to be 2 to comply with the BGZF specs

    static const unsigned MAX_BGZF_BLOCK_SIZE = 0xFFFF;

private:
#pragma pack(push, 1)
    struct Header
    {
        struct XFIELD
        {
            XFIELD(size_t compressedSize);

            unsigned char XLEN[2];

            unsigned char SI1;
            unsigned char SI2;
            unsigned char SLEN[2];
            unsigned char BSIZE[SUBFIELD_LENGTH];

            unsigned getXLEN() const {return unsigned(XLEN[0]) + XLEN[1] * 256;}
            unsigned getBSIZE() const {return unsigned(BSIZE[0]) + BSIZE[1] * 256;}

        };

        unsigned char ID1;
        unsigned char ID2;
        unsigned char CM;
        unsigned char FLG;
        unsigned char MTIME[4];
        unsigned char XFL;
        unsigned char OS;

        XFIELD xfield;
    };
#pragma pack(pop)

public:
    static size_t getHeaderSize() { return sizeof(Header); }

    static unsigned getBlockSize(char* buffer, bool failQuietly = false);

private:
    void prepareHeader();

    void rewriteHeader();

    // unverified number of bytes the compressor would add to data with compression level 0
    static const size_t gzip_junk = 41;
    // bgzf cannot handle blocks over 0xFFFF bytes long. assume 1/1 compression ratio for input data
    static const unsigned bgzf_buffer_size_ = MAX_BGZF_BLOCK_SIZE;
    static const unsigned short max_uncompressed_per_block_ = bgzf_buffer_size_ - gzip_junk;

    std::vector<char>&                            outputBuffer_;
    std::vector<char>                             internalBuffer_;
    bios::back_insert_device< std::vector<char> > sink_;
    boost::iostreams::gzip_compressor             compressor_;

    // as it is difficult to predict the size of compressed data, the buffering is based on
    // the amount of uncompressed data consumed. The assumption is that compressed data
    // will always require no more space than the uncompressed.
    size_t uncompressed_in_;

    size_t currentBlockBegin_;

    bool useBgzf_;
};

} // namespace io
} // namespace bcl2fastq


#endif // BCL2FASTQ_IO_GZIP_COMPRESSOR_HH
