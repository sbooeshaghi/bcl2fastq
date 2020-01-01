/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file BclFile.hh
 *
 * \brief Declaration of BCL file.
 *
 * \author Marek Balint
 * \author Mauricio Varea
 */


#ifndef BCL2FASTQ_DATA_BCLFILE_HH
#define BCL2FASTQ_DATA_BCLFILE_HH


#include <memory>

#include <boost/filesystem/path.hpp>

#include "common/Types.hh"
#include "data/FileReader.hh"


namespace bcl2fastq {

namespace io {
class GzipDecompressor;
}

namespace data {


/// \brief BCL file.
class BclFile : public BinaryFileReader
{
    typedef boost::error_info<struct tag_errmsg, std::string> errmsg_info;
public:

    /// \brief Default constructor.
    BclFile(const common::RawDataBuffer& data,
            io::GzipDecompressor&        decompressor,
            bool                         ignoreErrors,
            bool                         parseHeader,
            bool                         resetDecompressor = true,
            common::ClustersCount        defaultClustersCount = 0);

    virtual ~BclFile() { }

public:
    /// \brief Verify the bcl file exists and return the name. (single tile).
    /// \param inputDir Directory to look for the bcl file.
    /// \param laneNumber The lane number.
    /// \param tileNumber The tile number.
    /// \param cycleNumber The cycle number.
    /// \param ignoreErrors Supress errors opening file and/or reading its header.
    /// \param fileName Output parameter for the bcl file name.
    static bool getAndVerifyFileName(const boost::filesystem::path& inputDir,
                                     common::LaneNumber             laneNumber,
                                     common::TileNumber             tileNumber,
                                     common::CycleNumber            cycleNumber,
                                     bool                           ignoreErrors,
                                     boost::filesystem::path&       fileName);

    /// \brief Verify the bcl file exists and return the name. (aggregate tiles).
    /// \param inputDir Directory to look for the bcl file.
    /// \param laneNumber The lane number.
    /// \param cycleNumber The cycle number.
    /// \param ignoreErrors Supress errors opening file and/or reading its header.
    /// \param fileName Output parameter for the bcl file name.
    static bool getAndVerifyFileName(const boost::filesystem::path& inputDir,
                                     common::LaneNumber             laneNumber,
                                     common::CycleNumber            cycleNumber,
                                     bool                           ignoreErrors,
                                     boost::filesystem::path&       fileName);

    /// \brief Read bytes from file to buffer.
    /// \param targetBuffer Target buffer to read to.
    /// \param targetSize Maximum number of bytes to be read.
    /// \pre <tt>this->isOpen() == true</tt>
    /// \return Number of bytes read.
    std::streamsize read(char *targetBuffer, std::streamsize targetSize);

    /// \brief Seek to specified position.
    /// \param compressedOffset Offset in compressed data denoting start of BGZF block.
    /// \param uncompressedOffset Offset in uncompressed data of the BGZF block.
    /// \return true if seek was successful.
    bool seek(std::streamsize compressedOffset, std::streamsize uncompressedOffset);

private:
    std::streamsize read(std::istream &is, char *targetBuffer, std::streamsize targetSize);

    /// \brief Gets the file type string used for error messages.
    /// \return File type string
    virtual std::string getFileTypeStr() const { return "BCL"; }

    /// \brief GZip decompression filter.
    io::GzipDecompressor* gzipDecompressor_;

    /// \brief Whether we need to use the decompressor for reading or not.
    bool isCompressed_;
};

} // namespace data
} // namespace bcl2fastq

#endif // BCL2FASTQ_DATA_BCLFILE_HH

