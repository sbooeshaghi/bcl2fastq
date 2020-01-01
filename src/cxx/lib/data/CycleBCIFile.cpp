/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file CycleBCIFile.cpp
 *
 * \brief Implementation of Cycle BCI file.
 *
 * \author Marek Balint
 */


#include <boost/format.hpp>

#include "io/FileBufWithReopen.hh"
#include "data/CycleBCIFile.hh"


namespace bcl2fastq {


namespace data {


CycleBCIFile::CycleBCIFile(bool ignoreMissingBcls) 
: ignoreMissingBcls_(ignoreMissingBcls)
, path_()
, tilesCount_()
, buffer_()
, bufferPosition_()
{
}

bool CycleBCIFile::getAndVerifyFileName(const boost::filesystem::path& inputDir,
                                        common::LaneNumber             laneNumber,
                                        common::CycleNumber            cycleNumber,
                                        bool                           ignoreErrors,
                                        boost::filesystem::path&       fileName)
{
    fileName = boost::filesystem::path(
        inputDir
        /
        boost::filesystem::path((boost::format("L%03d") % laneNumber).str())
        /
        boost::filesystem::path((boost::format("%04d.bcl.bgzf.bci") % cycleNumber).str())
    );

    if (!boost::filesystem::exists(fileName))
    {
        std::string errStr("Unable to find BCI file: '" + fileName.string() + "'");
        if (ignoreErrors)
        {
            BCL2FASTQ_LOG(common::LogLevel::WARNING) << errStr << std::endl;
            return false;
        }
        else
        {
            BOOST_THROW_EXCEPTION(common::IoError(ENOENT, errStr));
        }

        fileName = boost::filesystem::path();
    }

    return true;
}

bool CycleBCIFile::openFile(
    const boost::filesystem::path &inputDir,
    common::LaneNumber laneNumber,
    common::CycleNumber cycleNumber
)
{
    if (!getAndVerifyFileName(inputDir, laneNumber, cycleNumber, ignoreMissingBcls_, path_))
    {
        return false;
    }

    io::FileBufWithReopen fileBuf(std::ios_base::in | std::ios_base::binary);
    fileBuf.reopen(path_, io::FileBufWithReopen::FadviseFlags::SEQUENTIAL_ONCE);
    int errnum = errno;
    if (!fileBuf.is_open())
    {
        BOOST_THROW_EXCEPTION(common::IoError(errnum, (boost::format("Unable to open cycle BCI file '%s' for reading") % path_.string()).str()));
    }

    Header header;
    std::streamsize headerLength = io::read(fileBuf, reinterpret_cast<char *>(&header), sizeof(header));
    errnum = errno;
    if (headerLength != sizeof(header))
    {
        BOOST_THROW_EXCEPTION(common::InputDataError(errnum, (boost::format("Unable to read header of cycle BCI file '%s': bytes_read=%d bytes_expected=%d") % path_.string() % headerLength % sizeof(header)).str()));
    }
    if (header.version_ != Header::VERSION)
    {
        const uint32_t version = header.version_;
        BOOST_THROW_EXCEPTION(common::InputDataError(errnum, (boost::format("Corrupted header of cycle BCI file '%s': header_version=%d") % path_.string() % version).str()));
    }
    tilesCount_ = header.tilesCount_;

    std::size_t bufferSize = sizeof(Record) * tilesCount_;
    buffer_.resize(bufferSize);
    std::streamsize bytesRead = io::read(fileBuf, reinterpret_cast<char *>(&*buffer_.begin()), bufferSize);
    errnum = errno;
    if (bytesRead != static_cast<std::streamsize>(bufferSize))
    {
        BOOST_THROW_EXCEPTION(common::InputDataError(errnum, (boost::format("Cycle BCI file '%s' truncated: bytes_read=%d bytes_expected=%d") % path_.string() % bytesRead % bufferSize).str()));
    }
    bufferPosition_ = buffer_.begin();

    return true;
}

bool CycleBCIFile::isOpen() const
{
    return !path_.empty();
}

boost::filesystem::path CycleBCIFile::getPath() const
{
    BCL2FASTQ_ASSERT_MSG(this->isOpen(), "File is not open");
    return path_;
}

common::TileNumber CycleBCIFile::getTilesCount() const
{
    BCL2FASTQ_ASSERT_MSG(this->isOpen(), "File is not open");
    return tilesCount_;
}

std::size_t CycleBCIFile::read(
    CycleBCIFile::Record *targetBuffer,
    std::size_t targetSize
)
{
    BCL2FASTQ_ASSERT_MSG(this->isOpen(), "File is not open");

    BCL2FASTQ_ASSERT_MSG((buffer_.end() - bufferPosition_) >= 0, "Invalid buffer size.");
    BCL2FASTQ_ASSERT_MSG((buffer_.end() - bufferPosition_) % sizeof(Record) == 0, "Invalid buffer size.");
    const std::size_t recordsToRead = std::min(
        static_cast<std::size_t>(buffer_.end() - bufferPosition_) / sizeof(Record),
        targetSize
    );

    const Record * const begin = reinterpret_cast<const Record * const>(&*bufferPosition_);
    const Record * const end = begin + recordsToRead;

    std::copy(begin, end, targetBuffer);
    bufferPosition_ += (recordsToRead * sizeof(Record));

    return recordsToRead;
}

CycleBCIFile::Record CycleBCIFile::getRecord(std::size_t tileIndex) const
{
    BCL2FASTQ_ASSERT_MSG(this->isOpen(), "File is not open");
    const Record * const ret = reinterpret_cast<const Record * const>(&*(buffer_.begin() + tileIndex*sizeof(Record)));
    return *ret;
}


const uint32_t CycleBCIFile::Header::VERSION = 0;


} // namespace data


} // namespace bcl2fastq


