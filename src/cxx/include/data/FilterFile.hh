/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file FilterFile.hh
 *
 * \brief Declaration of filter file.
 *
 * \author Marek Balint
 * \author Aaron Day
 */


#ifndef BCL2FASTQ_DATA_FILTERFILE_HH
#define BCL2FASTQ_DATA_FILTERFILE_HH


#include "data/FileReader.hh"
#include "data/ControlFile.hh"
#include "common/Types.hh"


namespace bcl2fastq {


namespace data {


/// \brief Filter file.
class FilterFileBase : public BinaryAllClustersFileReader 
{
public:
    /// \brief Filter file record type definition.
    struct Record
    {
    public:
        /// \brief Default constructor
        Record() : data_(0) { }

        /// \brief Constructor
        /// \param data Data to store
        Record(char data) : data_(data) { }

        /// \brief Filter data.
        char data_;
    };

    /// \brief Construct and open filter file.
    /// \param inputDir RaW data from the filter file.
    /// \param ignoreErrors Supress errors opening file and/or reading it.
    FilterFileBase(
        const common::RawDataBuffer& data,
        bool ignoreErrors,
        bool skipHeader
    );

    /// \brief Return the file type string. Used for logging message purposes.
    virtual std::string getFileTypeStr() const { return "Filter"; }

protected:
    /// \brief Filter V3 header type definition.
#pragma pack(push, 1)
    struct Header
    {
    public:
        /// \brief Filter file cluster count is 32 bits.
        typedef uint32_t FilterFileClustersCount;

        /// \brief Version field.
        uint32_t version_;

        /// \brief Clusters count. 32 bit
        FilterFileClustersCount clustersCountIn_;
    };
#pragma pack(pop)

    /// \brief Return the expected header version.
    virtual uint32_t getHeaderVersion() const { return 3; }

    /// \brief Read the header.
    virtual bool readHeader();

    /// \brief Validate the header. Throw an exception on failure.
    virtual bool validateHeader(const Header& header);

    bool skipHeader_;
};

class FilterFile : public FilterFileBase
{
public:
    FilterFile(const common::RawDataBuffer& data,
               bool                         ignoreErrors,
               bool                         skipHeader);

    /// \brief Read records from the file to buffer.
    /// \param targetBuffer Target buffer to read to.
    /// \param targetSize Maximum number of records to be read.
    /// \return Number of records read.
    /// \pre <tt>this->isOpen() == true</tt>
    virtual std::size_t read(
        Record *targetBuffer,
        std::size_t targetSize
    );

    static bool doesFileExist(const boost::filesystem::path& inputDir,
                              bool                           aggregateTilesFlag,
                              common::LaneNumber             laneNumber,
                              common::TileNumber             tileNumber,
                              size_t&                        headerSize,
                              boost::filesystem::path&       filePath);

private:
    /// \brief Return the number of bytes in a record.
    virtual std::size_t getRecordBytes() const { return sizeof(Record); }
};


} // namespace data
} // namespace bcl2fastq


#endif // BCL2FASTQ_DATA_FILTERFILE_HH


