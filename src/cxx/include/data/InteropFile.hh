/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file InteropFile.hh
 *
 * \brief Declaration of InterOp file.
 *
 * V2 (NovaSeq): http://illumina.github.io/interop/index_v2.html
 * Previous V1: http://illumina.github.io/interop/index_v1.html
 *
 * \author Marek Balint
 * \author Daniel Berard
 */


#ifndef BCL2FASTQ_DATA_INTEROPFILE_HH
#define BCL2FASTQ_DATA_INTEROPFILE_HH


#include <boost/shared_ptr.hpp>
#include <boost/filesystem/path.hpp>

#include "io/FileBufWithReopen.hh"


namespace bcl2fastq {


namespace data {


class InteropFileState;

typedef std::vector<char> CharVector;

class InteropRecord
{
public:

    typedef uint16_t LaneNumber;        // Lane number type definition.
    typedef uint32_t TileNumber;        // Tile number type definition.
    typedef uint16_t ReadNumber;        // Read number type definition.
    typedef uint64_t ClustersCount;     // Clusters count type definition.
    typedef uint16_t IndexNameLength;   // Index name length type definition.
    typedef uint16_t SampleIdLength;    // Sample id length type definition.
    typedef uint16_t ProjectNameLength; // Project name length type definition.

    InteropRecord();
    InteropRecord(
        LaneNumber lane,
        TileNumber tile,
        ReadNumber read,
        ClustersCount count,
        const std::string &index,
        const std::string &id,
        const std::string &project);

public:

    LaneNumber laneNumber;
    TileNumber tileNumber;
    ReadNumber readNumber;
    ClustersCount clustersCount;
    std::string indexName;
    std::string sampleId;
    std::string projectName;

}; // end class InteropRecord


/// \brief Base class for interop record writer.
class IInteropRecordWriter
{
public:
    IInteropRecordWriter() {}
    virtual ~IInteropRecordWriter() {}

    virtual size_t writeInteropRecord(CharVector &byteVector, const InteropRecord &interopRecord) = 0;

}; // end class IInteropRecordWriter


/// \brief Write binary interop IndexMetricsOut records to file.
/// \param CLUSTERS_COUNT Integer type used in IndexMetricsOut file format to store cluster clounts.
/// \param TILE Integer type used in IndexMetricsOut file format to store tile identifying number.
template<class CLUSTERS_COUNT, class TILE>
class InteropRecordWriter : public IInteropRecordWriter
{
public:

    InteropRecordWriter() : IInteropRecordWriter() {}
    virtual ~InteropRecordWriter() {}

    /// \brief Method to write binary interop IndexMetricsOut records to file.
    virtual size_t writeInteropRecord(CharVector &byteVector, const InteropRecord &interopRecord)
    {
        byteVector.clear();

        size_t byteCount = 0;
        byteCount += appendChars (byteVector, interopRecord.laneNumber);   // lane number
        byteCount += appendChars (byteVector, static_cast<TILE>(interopRecord.tileNumber)); // tile number
        byteCount += appendChars (byteVector, interopRecord.readNumber);   // read number
        byteCount += appendString(byteVector, interopRecord.indexName);    // barcode
        byteCount += appendChars (byteVector, static_cast<CLUSTERS_COUNT>(interopRecord.clustersCount)); // clusters count
        byteCount += appendString(byteVector, interopRecord.sampleId);     // sample name
        byteCount += appendString(byteVector, interopRecord.projectName);  // project name

        return byteCount;
    }

protected:

    static bool nonPrintable(const char c) {
        // True for special charactes before <space>, <delete> and after <delete>.
        return (c < 32 || c > 126);
    }

    template<class T>
    static size_t appendChars(CharVector& byteVector, T value) {
        byteVector.insert(byteVector.end(), reinterpret_cast<const char *>(&value), reinterpret_cast<const char *>(&value)+sizeof(value));
        return sizeof(value);
    }

    static size_t appendString(CharVector& byteVector, const std::string& value) {

        // sanitize string first
        //value.erase(std::remove_if(value.begin(), value.end(), nonPrintable), value.end());

        uint16_t length = static_cast<uint16_t>(value.size());
        size_t byteCount = appendChars(byteVector, length);
        byteVector.insert(byteVector.end(), value.begin(), value.end());

        byteCount += value.size();
        return byteCount;
    }

}; // end class InteropRecordWriter

typedef std::shared_ptr<IInteropRecordWriter> IInteropRecordWriterPtr;

class InteropFactory
{
public:

    enum Version { V1 = 1, V2 };
    static IInteropRecordWriterPtr createInteropRecordWriter(Version version);

}; // end class InteropFactory

/// \brief InterOp file.
class InteropFile
{
public:

    /// \brief InterOp file state smart pointer type definition.
    typedef boost::shared_ptr<InteropFileState> InteropFileStatePtr;

public:

    /// \brief constructor.
    InteropFile(InteropFactory::Version version = InteropFactory::V2);

public:

    /// \brief Open InterOp file.
    /// \param path Path to the InterOp file to be opened.
    /// \return true if successful
    bool open(const boost::filesystem::path &path);

public:

    /// \brief Get file path.
    /// \return File path.
    boost::filesystem::path getPath() const;

    /// \brief Write interop record.
    /// \param interop record data
    void writeInteropRecord(
        InteropRecord &interopData);

    InteropFactory::Version getVersion() const
    {
        return version_;
    }

private:

    /// \brief Write bytes from buffer to file.
    /// \param sourceBuffer Source buffer to read to.
    /// \param sourceSize Number of bytes to be written.
    void write(const char *sourceBuffer, std::streamsize sourceSize);

    /// \brief Filebuf.
    io::FileBufWithReopen fileBuf_;

    InteropFileStatePtr interopFileState_;
    InteropFactory::Version version_;
    IInteropRecordWriterPtr interopRecordWriter_;

}; // end class InteropFactory


/// \brief InterOp file internal state.
class InteropFileState
{
protected:

    /// \brief Constructor.
    InteropFileState(const boost::filesystem::path& filePath);

    /// \brief Virtual destructor.
    virtual ~InteropFileState() = 0;

protected:

    /// \brief Write bytes from buffer to file.
    /// \param os Output stream to write to.
    /// \param sourceBuffer Source buffer to read to.
    /// \param sourceSize Number of bytes to be written.
    virtual void write(std::ostream &os, const char *sourceBuffer, std::streamsize sourceSize) = 0;

private:

    /// \brief File path.
    boost::filesystem::path path_;

private:

    friend class InteropFile;
};

} // namespace data


} // namespace bcl2fastq


#endif // BCL2INTEROP_DATA_INTEROPFILE_HH


