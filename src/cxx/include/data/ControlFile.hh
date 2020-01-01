/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file ControlFile.hh
 *
 * \brief Declaration of control file.
 *
 * \author Aaron Day
 */


#ifndef BCL2FASTQ_DATA_CONTROLFILE_HH
#define BCL2FASTQ_DATA_CONTROLFILE_HH


#include <boost/filesystem/path.hpp>

#include "data/FileReader.hh"
#include "common/Types.hh"


namespace bcl2fastq {


namespace data {


/// \brief Control file.
class ControlFile : public BinaryAllClustersFileReader 
{
public:

    /// \brief Control file record type definition.
    struct Record
    {
    public:
        Record() : data_(0) { }

        Record(uint16_t data) : data_(data) { }

        /// \brief Control data.
        uint16_t data_;
    };

public:

    /// \brief Construct and open control file.
    /// \param filePath Path to control file.
    /// \param ignoreErrors Supress errors opening file and/or reading it.
    ControlFile(
        const common::RawDataBuffer& data,
        bool ignoreErrors = false
    );

    /// \brief Read records from the file to buffer.
    /// \param targetBuffer Target buffer to read to.
    /// \param targetSize Maximum number of records to be read.
    /// \return Number of records read.
    /// \pre <tt>this->isOpen() == true</tt>
    std::size_t read(
        ControlFile::Record *targetBuffer,
        std::size_t targetSize
    );

    static bool doesFileExist(const boost::filesystem::path& intensitiesDir,
                              common::LaneNumber             laneNumber,
                              common::TileNumber             tileNumber,
                              boost::filesystem::path&       fileName);

private:
     /// \brief Control V2 header type definition.
#pragma pack(push, 1)
    struct Header
    {
        /// \brief Control file clusters count is 32 bit
        typedef uint32_t ControlFileClustersCount;

        /// \brief Value of the version field.
        static const uint32_t VERSION;

        /// \brief Version field.
        uint32_t version_;

        /// \brief Clusters count.
        ControlFileClustersCount clustersCount_;
    };
#pragma pack(pop)

    /// \brief Read the header.
    virtual bool readHeader();

    /// \brief Validate the header. Throw an exception on failure.
    virtual bool validateHeader(const Header& header);

    /// \brief Return the number of bytes in a record.
    virtual std::size_t getRecordBytes() const { return sizeof(Record); }

    /// \brief Return the file type string. Used for logging message purposes.
    virtual std::string getFileTypeStr() const { return "Control"; }
};


} // namespace data
} // namespace bcl2fastq

#endif // BCL2FASTQ_DATA_CONTROLFILE_HH

