/**
* BCL to FASTQ file converter
* Copyright (c) 2007-2017 Illumina, Inc.
*
* This software is covered by the accompanying EULA
* and certain third party copyright/licenses, and any user of this
* source file is bound by the terms therein.
*
* \file ControlFile.cpp
*
* \brief Implementation of control file.
*
* \author Aaron Day
*/

#include "data/ControlFile.hh"

#include <boost/filesystem/operations.hpp>

namespace bcl2fastq
{
namespace data
{

ControlFile::ControlFile(const common::RawDataBuffer& data,
                         bool                         ignoreErrors /*= false*/)
    : FileReaderBase(data,
                     ignoreErrors)
    , BinaryAllClustersFileReader(data,
                                  ignoreErrors)
{
}

bool ControlFile::readHeader()
{
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

        clustersCount_ = static_cast<common::ClustersCount>(header.clustersCount_);
    }

    return true;
}

bool ControlFile::validateHeader(const ControlFile::Header& header)
{
    int errnum = errno;
    if (header.version_ != Header::VERSION)
    {
        const uint32_t version = header.version_;
        BOOST_THROW_EXCEPTION(common::InputDataError(errnum, (boost::format("Corrupted header of control file '%s': header_version=%d") % getPath().string() % version).str()));
    }

    return true;
}

std::size_t ControlFile::read(
    ControlFile::Record *targetBuffer,
    std::size_t targetSize
)
{
    if (!validateCondition(this->isOpen(), "File is not open")) { return 0; }

    return readBytes(reinterpret_cast<char*>(targetBuffer), targetSize);
}

const uint32_t ControlFile::Header::VERSION = 2;

bool ControlFile::doesFileExist(const boost::filesystem::path& intensitiesDir,
                                common::LaneNumber             laneNumber,
                                common::TileNumber             tileNumber,
                                boost::filesystem::path&       fileName)
{
    boost::filesystem::path laneDir((boost::format("L%03d") % laneNumber).str());
    boost::filesystem::path filePath = intensitiesDir / laneDir;

    fileName = filePath / (boost::format("s_%d_%d.control") % laneNumber % tileNumber).str();
    if (boost::filesystem::exists(fileName))
    {
        return true;
    }

    // zero padded
    fileName = filePath / (boost::format("s_%d_%04d.control") % laneNumber % tileNumber).str();
    if (boost::filesystem::exists(fileName))
    {
        return true;
    }

    return false;
}

} // namespace data
} // namespace bcl2fastq

