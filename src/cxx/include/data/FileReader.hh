/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file FileReader.hh
 *
 * \brief Declaration of FileReader file.
 *
 * \author Aaron Day
 */

#ifndef BCL2FASTQ_DATA_FILEREADER_HH
#define BCL2FASTQ_DATA_FILEREADER_HH

#include "common/Types.hh"
#include "common/Exceptions.hh"
#include "io/SyncFile.hh"

#include <boost/noncopyable.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/format.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/device/array.hpp>

namespace bcl2fastq {

namespace data {


/// \brief Base class for all file readers
class FileReaderBase
{
public:
    /// \brief Constructor
    /// \param data Raw data from file.
    /// \param ignoreErrors
    /// \param defaultClustersCount Default number of clusters
    FileReaderBase(const common::RawDataBuffer&   data,
                   bool                           ignoreErrors,
                   common::ClustersCount          defaultClustersCount = 0);

    /// \brief Destructor
    virtual ~FileReaderBase() = 0;

    /// \brief Get file path.
    /// \return File path.
    /// \pre <tt>this->isOpen() == true</tt>
    virtual const boost::filesystem::path& getPath() const;

    /// \brief Get number of clusters.
    /// \return Number of clusters in the file.
    /// \pre <tt>this->isOpen() == true</tt>
    virtual common::ClustersCount getClustersCount() const;

    /// \brief Check whether the file is open.
    /// \retval true File is open.
    /// \retval false File is not open.
    virtual bool isOpen() const { return !data_->path_.empty(); }

protected:

    /// \brief Validate the condition is true.
    /// \param warningMsg The warning message to use if condition is false
    /// \retval condition
    bool validateCondition(bool condition, const std::string& warningMsg) const;

    /// \brief Raw data
    const common::RawDataBuffer* data_;

    /// \brief Ignore errors opening file and/or reading its header.
    bool ignoreErrors_;

    /// \brief Number of clusters.
    common::ClustersCount clustersCount_;

    boost::iostreams::basic_array_source<char> inputSrc_;
    boost::iostreams::stream<boost::iostreams::basic_array_source<char>> istr_;
};

/// \brief File reading class. Virtual inheritance because of diamond inheritance.
class FileReader : public virtual FileReaderBase, private boost::noncopyable
{
public:

    /// \brief Constructor. Resource acquisition is initialization.
    /// \param data Raw data from file.
    /// \param ignoreErrors Supress errors opening file and/or reading its header.
    /// \param defaultClustersCount Number of clusters to assume in case of error opening the file
    /// and/or reading its header.
    FileReader(const common::RawDataBuffer&   data,
               bool                           ignoreErrors,
               common::ClustersCount          defaultClustersCount = 0);

    /// \brief Destructor
    virtual ~FileReader() = 0;

protected:

    /// \brief Return the file type string. Used for logging message purposes.
    virtual std::string getFileTypeStr() const = 0;

    /// \brief Log an error message on file io failure.
    /// \param bytesRead Number of bytes read from file.
    /// \param bytesExpected Number of bytes expected to be read from file.
    virtual void logError(std::streamsize bytesRead,
                          std::streamsize bytesExpected) const;

    std::streamsize readBytes(char* buffer,
                              uint32_t bytes);
};

/// \brief Binary File Reader. Base class for reading binary files.
class BinaryFileReader : public FileReader
{
public:

    /// \brief Constructor. Resource acquisition is initialization.
    /// \param data Raw data from file.
    /// \param ignoreErrors Supress errors opening file and/or reading its header.
    /// \param defaultClustersCount Number of clusters to assume in case of error opening the file
    /// and/or reading its header.
    BinaryFileReader(const common::RawDataBuffer& data,
                     bool                         ignoreErrors,
                     common::ClustersCount        defaultClustersCount =0);

    /// \brief Destructor
    virtual ~BinaryFileReader() = 0;

protected:
    /// \brief Read bytes from file to buffer.
    /// \param targetBuffer Target buffer to read to.
    /// \param targetSize Maximum number of bytes to be read.
    /// \return Number of bytes read.
    virtual std::streamsize read(char*           targetBuffer,
                                 std::streamsize targetSize);

    /// \brief Read the header. Template allows derived classes to implement their own header
    /// class without inheritance, which would increase the size of an instance.
    /// \param header Header to read.
    template<typename HEADER>
    bool readHeader(HEADER& header);

    /// \brief Read the header. The template readHeader method only loads the data from file
    /// into the HEADER instance. This method is implemented by derived classes to validate
    /// and do something with the data.
    virtual bool readHeader() { return true; }
};

/// \brief Binary file class that reads all clusters from the file on initialization.
class BinaryAllClustersFileReader : public BinaryFileReader
{
public:

    /// \brief Constructor
    /// \param data Raw data from file
    /// \param ignoreErrors Supress errors opening file and/or reading its header.
    BinaryAllClustersFileReader(const common::RawDataBuffer& data,
                                bool                         ignoreErrors);

protected:

    /// \brief Read the clusters from file
    /// \param buffer Buffer to read data into
    /// \param clustersCount Number of clusters to read (size of buffer)
    /// \retval true on success, false on failure.
    virtual bool readClusters(std::vector<char>& buffer,
                              common::ClustersCount clustersCount);

    /// \brief Validate the header. Throw an exception on failure.
    virtual void validateHeader();

    /// \brief Return the number of bytes in a record.
    virtual std::size_t getRecordBytes() const = 0;

    /// \brief Read all records in the file.
    virtual void readRecords() { }
};

} // namespace data

} // namespace bcl2fastq


#include "data/FileReader.hpp"


#endif // BCL2FASTQ_DATA_FILEREADER_HH
