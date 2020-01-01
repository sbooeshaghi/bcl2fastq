/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file PositionsFile.cpp
 *
 * \brief Implementation of positions file.
 *
 * \author Marek Balint
 * \author Mauricio Varea
 * \author Aaron Day
 */


#include "data/PositionsFile.hh"

#include <vector>
#include <algorithm>
#include <numeric>
#include <utility>

#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/device/array.hpp>
#include <boost/algorithm/string.hpp>


namespace bcl2fastq {
namespace data {


PositionsFile::PositionsFile(const common::RawDataBuffer& data,
                             bool ignoreErrors /*= false*/
)
    : FileReaderBase(data,
                     ignoreErrors)
{
}

namespace detail {

///////////////////////////////////////////////////////////////////////////////
// positions file type: LOCS
///////////////////////////////////////////////////////////////////////////////

/// \brief LOCS file.
class LocsFile : public PositionsFile, public BinaryAllClustersFileReader
{

public:

    static const boost::filesystem::path EXTENSION;

    /// \brief Constructor.
    /// \param data Raw data from the LOCS file.
    /// \param ignoreErrors Suppress errors opening file and/or reading it.
    explicit LocsFile(
        const common::RawDataBuffer& data,
        bool ignoreErrors,
        bool skipHeader
    );

    /// \brief Read records from the file to buffer.
    /// \param targetBuffer Target buffer to read to.
    /// \param targetSize Maximum number of records to be read.
    /// \return Number of records read.
    virtual std::size_t read(
        std::vector<data::PositionsFile::Record>& targetBuffer,
        std::size_t targetSize
    );

private:

    /// \brief LOCS header type definition.
#pragma pack(push, 1)
    struct Header
    {
        /// \brief Positions file clusters count.
        typedef uint32_t PositionsFileClustersCount;

        /// \brief Value of the first constant field.
        static const uint32_t FIELD1;

        /// \brief Value of the second constant field.
        static const float FIELD2;

        /// \brief The first constant field.
        uint32_t field1_;

        /// \brief The second constant field.
        float field2_;

        /// \brief Clusters count.
        PositionsFileClustersCount clustersCount_;
    };

    /// \brief LOCS record type definition
    struct Record
    {
        /// \brief X-coordinate.
        float x_;

        /// \brief y-coordinate.
        float y_;

    };
#pragma pack(pop)

    /// \brief Get the number of bytes in a record
    /// \return Number of bytes in a record
    virtual std::size_t getRecordBytes() const { return sizeof(Record); }

    /// \brief Read the header
    virtual bool readHeader(bool ignoreErrors);

    /// \brief Validate the header. Throw an exception on failure.
    virtual bool validateHeader(Header& header, bool ignoreErrors);

    /// \brief Gets the file type string used for error messages.
    /// \return File type string
    virtual std::string getFileTypeStr() const { return "LOCS"; }

    /// \brief Transform the Record type into a PositionsFile::Record
    /// \param record Record to transform
    static PositionsFile::Record transformRecord(const Record &record);
};

const boost::filesystem::path LocsFile::EXTENSION(".locs");

LocsFile::LocsFile(
    const common::RawDataBuffer& data,
    bool ignoreErrors,
    bool skipHeader
)
: FileReaderBase(data,
                 ignoreErrors)
, PositionsFile(data,
                ignoreErrors)
, BinaryAllClustersFileReader(data,
                              ignoreErrors)
{
    if (!skipHeader)
    {
        readHeader(ignoreErrors);
    }
}

bool LocsFile::readHeader(bool ignoreErrors)
{
    Header header;
    if (!BinaryAllClustersFileReader::readHeader(header) ||
        !validateHeader(header, ignoreErrors))
    {
        return false;
    }

    clustersCount_ = static_cast<common::ClustersCount>(header.clustersCount_);

    return true;
}

bool LocsFile::validateHeader(LocsFile::Header& header, bool ignoreErrors)
{
    int errnum = errno;
    if (header.field1_ != Header::FIELD1)
    {
        const uint32_t field1 = header.field1_;
        const std::string message = 
            (boost::format("Corrupted header of LOCS file '%s': header_field_1=%d")
                        % this->getPath().string() % field1).str();
        if (ignoreErrors) {
            BCL2FASTQ_LOG(common::LogLevel::WARNING) << message << std::endl;
        }
        else {
            BOOST_THROW_EXCEPTION(common::InputDataError(errnum, message));
        }

        clustersCount_ = 0;
        return false;
    }
    if (header.field2_ != Header::FIELD2)
    {
        const float field2 = header.field2_;
        const std::string message = 
            (boost::format("Corrupted header of LOCS file '%s': header_field_2=%f")
                        % this->getPath().string() % field2).str();

        if (ignoreErrors) {
            BCL2FASTQ_LOG(common::LogLevel::WARNING) << message << std::endl;
        }
        else {
            BOOST_THROW_EXCEPTION(common::InputDataError(errnum, message));
        }

        clustersCount_ = 0;
        return false;
    }

    clustersCount_ = static_cast<common::ClustersCount>(header.clustersCount_);

    return true;
}

std::size_t LocsFile::read(
    std::vector<data::PositionsFile::Record>& targetBuffer,
    std::size_t targetSize
)
{
    targetBuffer.resize(targetSize);
    std::size_t ret = readBytes(reinterpret_cast<char*>(targetBuffer.data()), targetSize*sizeof(Record));

    const Record * const begin = reinterpret_cast<const Record * const>(targetBuffer.data());
    const Record * const end = begin + targetSize;

    std::transform(begin, end, targetBuffer.data(), &LocsFile::transformRecord);

    return ret;
}

const uint32_t LocsFile::Header::FIELD1 = 1;

const float LocsFile::Header::FIELD2 = 1.0;

PositionsFile::Record LocsFile::transformRecord(const Record &record)
{
    PositionsFile::Record ret;
    ret.x_ = static_cast<PositionsFile::Record::ClusterCoordinate>(record.x_ * 10.0 + 1000.5);
    ret.y_ = static_cast<PositionsFile::Record::ClusterCoordinate>(record.y_ * 10.0 + 1000.5);

    return ret;
}


///////////////////////////////////////////////////////////////////////////////
// positions file type: CLOCS
///////////////////////////////////////////////////////////////////////////////

/// \brief CLOCS file.
class ClocsFile : public PositionsFile, public BinaryAllClustersFileReader
{
public:

    static const boost::filesystem::path EXTENSION;

    /// \brief Constructor.
    /// \param data Raw data from the CLOCS file.
    /// \param ignoreErrors Suppress errors opening file and/or reading it.
    ClocsFile(
        const common::RawDataBuffer& data,
        bool ignoreErrors,
        bool skipHeader
    );

    /// \brief Read records from the file to buffer.
    /// \param targetBuffer Target buffer to read to.
    /// \param targetSize Maximum number of records to be read.
    /// \return Number of records read.
    virtual std::size_t read(
        std::vector<data::PositionsFile::Record>& targetBuffer,
        std::size_t targetSize
    );

private:

    /// \brief CLOCS header type definition.
#pragma pack(push, 1)
    struct Header
    {
    public:
        /// \brief supported version number
        static const uint8_t VERSION;

        /// \brief Version of the format.
        uint8_t version_;

        /// \brief Number of tiles
        uint32_t numTiles_;
    };

    struct Record
    {
    public:

        /// \brief x-offset
        uint8_t xOffset_;

        /// \brief y-offset
        uint8_t yOffset_;

    };
#pragma pack(pop)

private:

    /// \brief Get the number of bytes in a record
    /// \return Number of bytes in a record
    virtual std::size_t getRecordBytes() const { return sizeof(Record); }

    /// \brief Validate the header. Throw an exception on failure.
    virtual bool readHeader(bool ignoreErrors);

    /// \brief Validate the header. Throw an exception on failure.
    virtual bool validateHeader(Header& header, bool ignoreErrors);

//    /// \brief Read the number of clusters in the tile.
//    uint8_t readClusterCount();

    /// \brief Gets the file type string used for error messages.
    /// \return File type string
    virtual std::string getFileTypeStr() const { return "CLOCS"; }

    /// \brief Transform the Record type into a PositionsFile::Record
    /// \param record Record to transform
    PositionsFile::Record transformRecord(const Record &record);

private:

    /// \brief Number of tiles
    uint32_t numTiles_;

    uint32_t tileIndex_;
};

const boost::filesystem::path ClocsFile::EXTENSION(".clocs");

ClocsFile::ClocsFile(
    const common::RawDataBuffer& data,
    bool ignoreErrors,
    bool skipHeader
)
: FileReaderBase(data,
                 ignoreErrors)
, PositionsFile(data,
                ignoreErrors)
, BinaryAllClustersFileReader(data,
                              ignoreErrors)
, numTiles_(0)
, tileIndex_(0)
{
    if (!skipHeader)
    {
        readHeader(ignoreErrors);
    }
}

bool ClocsFile::readHeader(bool ignoreErrors)
{
    Header header;
    if (!BinaryAllClustersFileReader::readHeader(header) ||
        !validateHeader(header, ignoreErrors))
    {
        return false;
    }

    numTiles_ = header.numTiles_;

    return true;
}

bool ClocsFile::validateHeader(ClocsFile::Header& header,
                               bool ignoreErrors)
{
    int errnum = errno;
    if (header.version_ != Header::VERSION) {

        const std::string message = 
            (boost::format("Corrupted header of CLOCS file '%s': header_version=%d")
                    % this->getPath().string() % header.version_).str();

        if (ignoreErrors) {
           BCL2FASTQ_LOG(common::LogLevel::WARNING) << message << std::endl;
        }
        else {
            BOOST_THROW_EXCEPTION(common::InputDataError(errnum, message));
        }

        return false;
    }

    return true;
}

//uint8_t ClocsFile::readClusterCount()
//{
//    uint8_t clustersCount = 0;
//
//    std::streamsize clustersCountLength =
//        BinaryFileReader::read(reinterpret_cast<char *>(&clustersCount),
//                               sizeof(clustersCount));
//
//    if (clustersCountLength != sizeof(clustersCount))
//    {
//        logError(clustersCountLength,
//                 sizeof(clustersCount));
//    }
//
//    return clustersCount;
//}

std::size_t ClocsFile::read(
    std::vector<data::PositionsFile::Record>& targetBuffer,
    std::size_t targetSize
)
{
    targetBuffer.clear();

    size_t numClustersRead = 0;
    if (!validateCondition(numTiles_ > 0, "No tiles read from clocs file.")) { return 0; }

    // Skip first 5 bytes, which are the version number and number of bins in the file
    size_t byteIndex = 5;
    for (uint32_t binNumber = 0; binNumber < numTiles_; ++binNumber)
    {
        if (!validateCondition(data_->size() > byteIndex, "Malformed clocs file.")) { return 0; }
        unsigned char numClusters = static_cast<unsigned char>(data_->at(byteIndex++));

        if (numClusters == 0)
        {
            continue;
        }

        targetBuffer.resize(targetBuffer.size()+numClusters);
        if (!validateCondition(data_->size() >= byteIndex + numClusters*2, "Malformed clocs file.")) { return 0; }
        for (uint32_t i = 0; i < numClusters; ++i)
        {
            auto& record = targetBuffer[numClustersRead + i];
            record.x_ = static_cast<PositionsFile::Record::ClusterCoordinate>( static_cast<uint8_t>((*data_)[byteIndex++]) + (binNumber % 82)*250) + 1000.5;
            record.y_ = static_cast<PositionsFile::Record::ClusterCoordinate>( static_cast<uint8_t>((*data_)[byteIndex++]) + (binNumber / 82)*250) + 1000.5;
        }

        numClustersRead += numClusters;
    }

    return numClustersRead;
}


const uint8_t ClocsFile::Header::VERSION = 0x01;

PositionsFile::Record ClocsFile::transformRecord(const Record &record)
{
    PositionsFile::Record ret;
    ret.x_ = static_cast<PositionsFile::Record::ClusterCoordinate>(record.xOffset_ + (tileIndex_ % 82)*250) + 1000.5;
    ret.y_ = static_cast<PositionsFile::Record::ClusterCoordinate>(record.yOffset_ + (tileIndex_ / 82)*250) + 1000.5;

    return ret;
}


class PosFile : public PositionsFile, FileReader
{
public:
    static const boost::filesystem::path EXTENSION;

    PosFile(const common::RawDataBuffer& data,
            bool                         ignoreErrors);

    /// \brief Read records from file to buffer.
    /// \param targetBuffer Target buffer to read to.
    /// \param targetSize Maximum number of records to be read.
    /// \return Number of records read.
    virtual std::size_t read(
        std::vector<data::PositionsFile::Record>& targetBuffer,
        std::size_t targetSize
    );

private:
    /// \brief Gets the file type string used for error messages.
    /// \return File type string
    virtual std::string getFileTypeStr() const { return "POS"; }
};

const boost::filesystem::path PosFile::EXTENSION(".txt");

PosFile::PosFile(const common::RawDataBuffer& data,
                 bool                         ignoreErrors)
: FileReaderBase(data,
                 ignoreErrors)
, PositionsFile(data,
                ignoreErrors)
, FileReader(data,
             ignoreErrors)
{
}

std::size_t PosFile::read(
    std::vector<data::PositionsFile::Record>& targetBuffer,
    std::size_t targetSize
)
{
    targetBuffer.clear();
    targetBuffer.reserve(targetSize);

    std::size_t recordsRead = 0;
    static const std::string errMsgBegin = "Ignoring read failure for POS file '";
    try
    {
        std::string line;
        while (std::getline(istr_, line) && recordsRead < targetSize)
        {
            std::vector<std::string> strs;
            boost::split(strs, line, boost::is_any_of("\t "));

            if (strs.size() != 2)
            {
                std::string errMsg = "Formatting error for POS file: " + this->getPath().string();
                if (ignoreErrors_)
                {
                    BCL2FASTQ_LOG(common::LogLevel::WARNING) << errMsg << std::endl;
                    return recordsRead;
                }
                else
                {
                    int errnum = errno;
                    BOOST_THROW_EXCEPTION(common::InputDataError(errnum, errMsg));
                }
            }

            PositionsFile::Record record;
            record.x_ = boost::lexical_cast<PositionsFile::Record::ClusterCoordinate>(strs[0]);
            record.y_ = boost::lexical_cast<PositionsFile::Record::ClusterCoordinate>(strs[1]);
            targetBuffer.push_back(record);

            ++recordsRead;
        }
    }
    CATCH_AND_IGNORE(common::IoError, errMsgBegin)
    CATCH_AND_IGNORE(std::ios_base::failure, errMsgBegin)
    CATCH_AND_IGNORE_ALL(errMsgBegin)

    return recordsRead;
}

} // namespace detail


bool PositionsFileFactory::doesFileExist(
    const boost::filesystem::path& intensitiesDir,
    bool                           aggregateTilesFlag,
    bool                           isPatternedFlowcell,
    common::LaneNumber             laneNumber,
    common::TileNumber             tileNumber,
    size_t&                        headerSize,
    boost::filesystem::path&       positionsFilePath
)
{
    if (binaryFilePathExists(intensitiesDir,
                             aggregateTilesFlag,
                             isPatternedFlowcell,
                             laneNumber,
                             tileNumber,
                             detail::LocsFile::EXTENSION,
                             positionsFilePath))
    {
        headerSize = 12;
        return true;
    }

    if (binaryFilePathExists(intensitiesDir,
                             aggregateTilesFlag,
                             isPatternedFlowcell,
                             laneNumber,
                             tileNumber,
                             detail::ClocsFile::EXTENSION,
                             positionsFilePath))
    {
        headerSize = 5;
        return true;
    }

    if (posFilePathExists(intensitiesDir,
                          laneNumber,
                          tileNumber,
                          positionsFilePath))
    {
        headerSize = 0;
        return true;
    }

    return false;
}

bool PositionsFileFactory::binaryFilePathExists(
    const boost::filesystem::path& intensitiesDir,
    bool                           aggregateTilesFlag,
    bool                           isPatternedFlowcell,
    common::LaneNumber             laneNumber,
    common::TileNumber             tileNumber,
    const boost::filesystem::path& fileExtension,
    boost::filesystem::path&       positionsFilePath
)
{
    if( isPatternedFlowcell )
    {
        positionsFilePath = boost::filesystem::path(
            intensitiesDir
            /
            boost::filesystem::path("s" + fileExtension.string())
        );
    } else if( aggregateTilesFlag ) {
        positionsFilePath = boost::filesystem::path(
            intensitiesDir
            /
            boost::filesystem::path((boost::format("L%03d") % laneNumber).str())
            /
            boost::filesystem::path((boost::format("s_%d" + fileExtension.string()) % laneNumber).str())
        );
    } else {
        positionsFilePath = boost::filesystem::path(
            intensitiesDir
            /
            boost::filesystem::path((boost::format("L%03d") % laneNumber).str())
            /
            boost::filesystem::path((boost::format("s_%d_%d" + fileExtension.string()) % laneNumber % tileNumber).str())
        );
    }

    return boost::filesystem::exists(positionsFilePath);
}

bool PositionsFileFactory::posFilePathExists(
    const boost::filesystem::path& intensitiesDir,
    common::LaneNumber             laneNumber,
    common::TileNumber             tileNumber,
    boost::filesystem::path&       positionsFilePath
)
{
    positionsFilePath = boost::filesystem::path(
        intensitiesDir
        /
        boost::filesystem::path((boost::format("s_%d_%03d_pos.txt") % laneNumber % tileNumber).str())
    );

    return boost::filesystem::exists(positionsFilePath);
}

std::streamsize PositionsFileFactory::read(const common::RawDataBuffer& inputData,
                                           bool ignoreErrors,
                                           bool skipHeader,
                                           std::vector<data::PositionsFile::Record>& targetBuffer,
                                           std::size_t targetSize)
{
    auto extension = inputData.path_.extension();
    if (extension == detail::LocsFile::EXTENSION)
    {
        detail::LocsFile posFile(inputData,
                                 ignoreErrors,
                                 skipHeader);

        return posFile.read(targetBuffer, targetSize > 0 ? targetSize : posFile.getClustersCount());

    }
    else if (extension == detail::ClocsFile::EXTENSION)
    {
        detail::ClocsFile posFile(inputData,
                                  ignoreErrors,
                                  skipHeader);

        return posFile.read(targetBuffer, targetSize);
    }
    else if (extension == detail::PosFile::EXTENSION)
    {
        detail::PosFile posFile(inputData,
                                ignoreErrors);

        posFile.read(targetBuffer, targetSize);
    }

    return 0;
}

} // namespace data
} // namespace bcl2fastq

