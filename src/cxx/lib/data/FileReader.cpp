/**
* BCL to FASTQ file converter
* Copyright (c) 2007-2017 Illumina, Inc.
*
* This software is covered by the accompanying EULA
* and certain third party copyright/licenses, and any user of this
* source file is bound by the terms therein.
*
* \file FileReader.cpp
*
* \brief Implementation of FileReader.
*
* \author Aaron Day
*/

#include "data/FileReader.hh"

#include "common/Debug.hh"
#include "common/Exceptions.hh"


namespace bcl2fastq
{
namespace data
{

FileReaderBase::FileReaderBase(const common::RawDataBuffer& data,
                               bool                         ignoreErrors,
                               common::ClustersCount        defaultClustersCount /*= 0*/)
    : data_(&data)
    , ignoreErrors_(ignoreErrors)
    , clustersCount_(defaultClustersCount)
    , inputSrc_(data_->data(), data_->size())
    , istr_(inputSrc_)
{
}

FileReaderBase::~FileReaderBase()
{
}

const boost::filesystem::path& FileReaderBase::getPath() const
{
    validateCondition(this->isOpen(), "File is not open");
    return data_->path_;
}

common::ClustersCount FileReaderBase::getClustersCount() const
{
    validateCondition(this->isOpen(), "File is not open");
    return clustersCount_;
}

bool FileReaderBase::validateCondition(bool condition, const std::string& warningMsg) const
{
    if (!condition)
    {
        if (ignoreErrors_)
        {
            BCL2FASTQ_LOG(common::LogLevel::WARNING) << warningMsg << std::endl;
        }
        else
        {
            BCL2FASTQ_ASSERT_MSG(condition, warningMsg);
        }
    }

    return condition;
}

FileReader::FileReader(const common::RawDataBuffer& data,
                       bool                         ignoreErrors,
                       common::ClustersCount        defaultClustersCount /*=0*/)
    : FileReaderBase(data,
                     ignoreErrors,
                     defaultClustersCount)
{
}

FileReader::~FileReader()
{
}

void FileReader::logError(std::streamsize bytesRead,
                          std::streamsize bytesExpected) const
{
    int errnum = errno;

    boost::format errFormat(boost::format(
        getFileTypeStr() + " file '%s' truncated: bytes_read=%d bytes_expected=%d") %
        this->getPath().string() % bytesRead % bytesExpected);

    if (ignoreErrors_)
    {
        BCL2FASTQ_LOG(common::LogLevel::WARNING)
            << errFormat
            << ":" << std::strerror(errnum) << " (" << errnum << ")" << std::endl;
    }
    else
    {
        BOOST_THROW_EXCEPTION(common::InputDataError(errnum,
                                                     errFormat.str()));
    }
}

std::streamsize FileReader::readBytes(char* buffer,
                                      uint32_t bytes)
{
    std::streamsize bytesRead = 0;
    istr_.read(buffer, bytes);
    bytesRead = istr_.gcount();
    if (!istr_ || bytesRead < bytes)
    {
        logError(bytesRead,
                 bytes);
    }

    return bytesRead;
}

BinaryFileReader::BinaryFileReader(const common::RawDataBuffer& data,
                                   bool                         ignoreErrors,
                                   common::ClustersCount        defaultClustersCount /*=0*/)
    : FileReaderBase(data,
                     defaultClustersCount)
    , FileReader(data,
                 ignoreErrors,
                 defaultClustersCount)
{
}

BinaryFileReader::~BinaryFileReader()
{
}

std::streamsize BinaryFileReader::read(char *targetBuffer, std::streamsize targetSize)
{
    if (istr_.read(targetBuffer, targetSize))
    {
        return targetSize;
    }

    return istr_.gcount();
}

BinaryAllClustersFileReader::BinaryAllClustersFileReader(const common::RawDataBuffer& data,
                                                         bool                         ignoreErrors)
    : FileReaderBase(data,
                     ignoreErrors)
    , BinaryFileReader(data,
                       ignoreErrors)
{
}

bool BinaryAllClustersFileReader::readClusters(std::vector<char>& buffer,
                                               common::ClustersCount clustersCount)
{
    std::streamsize bufferSize = getRecordBytes() * clustersCount;
    std::streamsize bytesRead = readBytes(buffer.data(), bufferSize);

    if (bufferSize != bytesRead)
    {
        buffer.resize(bytesRead > 0 ? bytesRead - bytesRead % getRecordBytes() : 0);
        return false;
    }

    return true;
}

void BinaryAllClustersFileReader::validateHeader()
{
    // Let derived classes decide if they want to implement this
}

} // namespace data
} // namespace bcl2fastq

