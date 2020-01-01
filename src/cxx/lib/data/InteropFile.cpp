/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file InteropFile.cpp
 *
 * \brief Implementation of InterOp file.
 *
 * \author Marek Balint
 */


#include <boost/format.hpp>

#include "common/Debug.hh"
#include "common/Logger.hh"
#include "common/SystemCompatibility.hh"
#include "data/InteropFile.hh"


namespace bcl2fastq {


namespace data {


InteropFileState::InteropFileState(const boost::filesystem::path& filePath)
: path_(filePath)
{
}

InteropFileState::~InteropFileState()
{
}

void InteropFileState::write(std::ostream &os, const char *sourceBuffer, std::streamsize sourceSize)
{
    os.write(sourceBuffer, sourceSize);
}


namespace detail {


/// \brief Default InterOp file internal state.
class DefaultInteropFileState : public InteropFileState
{
public:
    DefaultInteropFileState(const boost::filesystem::path& filePath) : InteropFileState(filePath) { }

private:
    virtual void write(std::ostream &os, const char *sourceBuffer, std::streamsize sourceSize);
};

void DefaultInteropFileState::write(std::ostream &os, const char *sourceBuffer, std::streamsize sourceSize)
{
    InteropFileState::write(os, sourceBuffer, sourceSize);
}


} // namespace detail


InteropFile::InteropFile(InteropFactory::Version version)
: fileBuf_(std::ios_base::out | std::ios_base::app | std::ios_base::binary)
, version_(version)
{
    interopRecordWriter_ = InteropFactory::createInteropRecordWriter(version);
}

bool InteropFile::open(const boost::filesystem::path &path)
{
    io::FileBufWithReopen fileBuf(std::ios_base::out | std::ios_base::trunc);
    // I don't like the addition of this try block. The entire interop file writing implementaion is
    // poorly designed. All we need to do is write a number to a single file. This is over kill.
    try
    {
        fileBuf.reopen(path.native().c_str(), io::FileBufWithReopen::FadviseFlags::SEQUENTIAL_ONCE);
    }
    catch (common::IoError)
    {
        //pass
    }

    if (!fileBuf.is_open())
    {
        BCL2FASTQ_LOG(common::LogLevel::WARNING) << "Unable to open Interop file '" << path.string() << "' for writing" << std::endl;
        return false;
    }

    interopFileState_ =  InteropFileStatePtr(new detail::DefaultInteropFileState(path));
    if (!interopFileState_) {

        BCL2FASTQ_LOG(common::LogLevel::WARNING) << "Failed to open Interop file '" << path.string() << "' for writing" << std::endl;
        return false;
    }

    BCL2FASTQ_LOG(common::LogLevel::INFO) << "Created InterOp file '" << path << std::endl;
    CharVector byteVector(1);
    byteVector[0] = (char)(version_);

    write(&*byteVector.begin(), 1);//sizeof(versionChar));

    return true;
}

boost::filesystem::path InteropFile::getPath() const
{
    if (!interopFileState_) {
        return boost::filesystem::path();
    }

    return interopFileState_->path_;
}

void InteropFile::writeInteropRecord(
    data::InteropRecord &interopData)
{
    CharVector byteVector;
    size_t byteCount = interopRecordWriter_->writeInteropRecord(byteVector, interopData);

    if (byteVector.size() != byteCount) {
         BCL2FASTQ_LOG(common::LogLevel::WARNING)
             << "Inconsistent byte count writting interop index metrics entry \""
             << interopData.laneNumber << ":" << interopData.tileNumber << ":" << interopData.readNumber << ":"
             << interopData.indexName << ":" << interopData.clustersCount << ":"
             << interopData.sampleId << ":" << interopData.projectName << "\"" << std::endl;
    }

    write(&*byteVector.begin(), byteVector.size());
}

void InteropFile::write(const char *sourceBuffer, std::streamsize sourceSize)
{
    fileBuf_.reopen(interopFileState_->path_.native().c_str(), io::FileBufWithReopen::FadviseFlags::SEQUENTIAL_ONCE);
    BCL2FASTQ_ASSERT_MSG(fileBuf_.is_open(), "Unable to reopen InterOp file '" << interopFileState_->path_ << "' for append");
    std::ostream os(&fileBuf_);

    interopFileState_->write(os, sourceBuffer, sourceSize);
}

InteropRecord::InteropRecord() {}
InteropRecord::InteropRecord(
    LaneNumber lane,
    TileNumber tile,
    ReadNumber read,
    ClustersCount cCount,
    const std::string &index,
    const std::string &id,
    const std::string &project)
: laneNumber(lane)
, tileNumber(tile)
, readNumber(read)
, clustersCount(cCount)
, indexName(index)
, sampleId(id)
, projectName(project)
{
}

data::IInteropRecordWriterPtr InteropFactory::createInteropRecordWriter(InteropFactory::Version version)
{
    if (version == V1) {
        return std::make_shared< InteropRecordWriter</*CLUSTERS_COUNT*/uint32_t, /*TILE*/uint16_t> >();
    }

    // V2
    return std::make_shared< InteropRecordWriter</*CLUSTERS_COUNT*/uint64_t, /*TILE*/uint32_t> >();
}

} // namespace data


} // namespace bcl2fastq


