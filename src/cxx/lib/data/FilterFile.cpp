/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file FilterFile.cpp
 *
 * \brief Implementation of filter file.
 *
 * \author Marek Balint
 * \author Mauricio Varea
 */


#include <boost/format.hpp>
#include <boost/filesystem/operations.hpp>

#include "io/FileBufWithReopen.hh"
#include "data/FilterFile.hh"


namespace bcl2fastq {


namespace data {


FilterFileBase::FilterFileBase(const common::RawDataBuffer& data,
                               bool                         ignoreErrors,
                               bool                         skipHeader)
    : FileReaderBase(data,
                     ignoreErrors)
    , BinaryAllClustersFileReader(data,
                                  ignoreErrors)
    , skipHeader_(skipHeader)
{
}

bool FilterFileBase::readHeader()
{
    if (skipHeader_)
    {
        return true;
    }

    if (!BinaryAllClustersFileReader::readHeader(clustersCount_))
    {
        return false;
    }

    if (clustersCount_ == 0)
    {
        Header header;
        if (!BinaryAllClustersFileReader::readHeader(header))
        {
            return false;
        }

        if (!validateHeader(header))
        {
            return false;
        }

        clustersCount_ = static_cast<common::ClustersCount>(header.clustersCountIn_);
    }

    return true;
}

bool FilterFileBase::validateHeader(const FilterFile::Header& header)
{
    int errnum = errno;
    if (header.version_ != getHeaderVersion())
    {
        const uint32_t version = header.version_;
        BOOST_THROW_EXCEPTION(common::InputDataError(errnum, (boost::format("Corrupted header of filter file '%s': header_version=%d") % getPath().string() % version).str()));

        return false;
    }

    return true;
}

FilterFile::FilterFile(const common::RawDataBuffer& data,
                       bool                         ignoreErrors,
                       bool                         skipHeader)
    : FileReaderBase(data,
                     ignoreErrors)
    , FilterFileBase(data,
                     ignoreErrors,
                     skipHeader)
{
    if (!skipHeader)
    {
        readHeader();
    }
}

std::size_t FilterFile::read(
    FilterFile::Record *targetBuffer,
    std::size_t targetSize
)
{
    if (!validateCondition(this->isOpen(), "File is not open")) { return 0; }

    return readBytes(reinterpret_cast<char*>(targetBuffer), targetSize);
}


bool FilterFile::doesFileExist(const boost::filesystem::path& inputDir,
                               bool                           aggregateTilesFlag,
                               common::LaneNumber             laneNumber,
                               common::TileNumber             tileNumber,
                               size_t&                        headerSize,
                               boost::filesystem::path&       filePath)
{
    headerSize = 12;
    if( aggregateTilesFlag ) {
        filePath =
            inputDir
            / boost::filesystem::path((boost::format("L%03d") % laneNumber).str())
            / boost::filesystem::path((boost::format("s_%d.filter") % laneNumber).str());
    } else {
        boost::filesystem::path fileName((boost::format("s_%d_%d.filter") % laneNumber % tileNumber).str());
        boost::filesystem::path fileNameZeroPadded((boost::format("s_%d_%04d.filter") % laneNumber % tileNumber).str());
        boost::filesystem::path laneDir((boost::format("L%03d") % laneNumber).str());

        filePath = inputDir / laneDir / fileName;

        if (!boost::filesystem::exists(filePath))
        {
            filePath = inputDir / fileName;
            headerSize = 4;
        }

        if (!boost::filesystem::exists(filePath))
        {
            filePath = inputDir / laneDir / fileNameZeroPadded;
        }

        if (!boost::filesystem::exists(filePath))
        {
            filePath = inputDir / fileNameZeroPadded;
            headerSize = 4;
        }
    }

    return boost::filesystem::exists(filePath);
}


} // namespace data
} // namespace bcl2fastq


