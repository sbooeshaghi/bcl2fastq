/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file CycleBCIFile.hh
 *
 * \brief Declaration of Cycle BCI file.
 *
 * \author Marek Balint
 */


#ifndef BCL2FASTQ_DATA_CYCLEBCIFILE_HH
#define BCL2FASTQ_DATA_CYCLEBCIFILE_HH


#include <boost/filesystem/path.hpp>

#include "common/Types.hh"


namespace bcl2fastq {


namespace data {


/// \brief Cycle BCI file.
class CycleBCIFile
{
public:

    /// \brief Cycle BCI file record type definition.
#pragma pack(push, 1)
    struct Record
    {
    public:

        /// \brief Uncompressed offset.
        uint64_t uncompressedOffset : 16;

        /// \brief Compressed offset.
        uint64_t compressedOffset : 48;
    };
#pragma pack(pop)

public:

    /// \brief Construction
    CycleBCIFile(bool ignoreMissingBcls);

    /// \brief Open filter file.
    /// \param inputDir Path to input directory.
    /// \param laneNumber Lane number.
    /// \param cycleNumber Cycle number.
    bool openFile(
        const boost::filesystem::path &inputDir,
        common::LaneNumber laneNumber,
        common::CycleNumber cycleNumber
    );

    static bool getAndVerifyFileName(
        const boost::filesystem::path &inputDir,
        common::LaneNumber laneNumber,
        common::CycleNumber cycleNumber,
        bool ignoreErrors,
        boost::filesystem::path &fileName
    );

    /// \brief Check whether the file is open.
    /// \retval true File is open.
    /// \retval false File is not open.
    bool isOpen() const;

    /// \brief Get file path.
    /// \return File path.
    /// \pre <tt>this->isOpen() == true</tt>
    boost::filesystem::path getPath() const;

    /// \brief Get number of tiles.
    /// \return Number of tiles in the file.
    /// \pre <tt>this->isOpen() == true</tt>
    common::TileNumber getTilesCount() const;

    /// \brief Read records from the file to buffer.
    /// \param targetBuffer Target buffer to read to.
    /// \param targetSize Maximum number of records to be read.
    /// \return Number of records read.
    /// \pre <tt>this->isOpen() == true</tt>
    std::size_t read(
        CycleBCIFile::Record *targetBuffer,
        std::size_t targetSize
    );

    /// \brief Get record for given tile index.
    /// \param tileIndex Tile index.
    /// \return Record for given tile index.
    /// \pre <tt>this->isOpen() == true</tt>
    Record getRecord(std::size_t tileIndex) const;

private:

    /// \brief Cycle BCIndex V0 header type definition.
#pragma pack(push, 1)
    struct Header
    {
    public:

        /// \brief Value of the version field.
        static const uint32_t VERSION;

    public:

        /// \brief Version field.
        uint32_t version_;

        /// \brief Clusters count.
        uint32_t tilesCount_;

    };
#pragma pack(pop)

private:

    /// \brief Ignore missing bcl/bci data
    bool ignoreMissingBcls_ = false; 

    /// \brief File path.
    boost::filesystem::path path_;

    /// \brief Number of clusters.
    common::TileNumber tilesCount_;

    /// \brief Internal buffer.
    std::vector<char> buffer_;

    /// \brief Current position in internal buffer.
    std::vector<char>::const_iterator bufferPosition_;

};


} // namespace data


} // namespace bcl2fastq


#endif // BCL2FASTQ_LAYOUT_BCINDEX_HH


