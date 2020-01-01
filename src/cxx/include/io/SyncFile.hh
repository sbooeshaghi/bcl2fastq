/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file SyncFile.hh
 *
 * \brief Declaration of Synchronized file.
 *
 * \author Aaron Day
 */

#ifndef BCL2FASTQ_IO_SYNC_FILE_HH
#define BCL2FASTQ_IO_SYNC_FILE_HH

#include "io/FileBufWithReopen.hh"
#include "io/GzipCompressor.hh"
#include "common/Types.hh"

#include <boost/format.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/noncopyable.hpp>
#include <condition_variable>
#include <atomic>
#include <mutex>
#include <string>
#include <vector>

namespace bcl2fastq {
namespace io {

// This interface class is here so that most code can use a pointer to this
// base class rather than use a class that touches the file system.
// This makes writing unit tests easier.
class FileReaderBase : private boost::noncopyable
{
public:
    FileReaderBase();

    FileReaderBase(const boost::filesystem::path& path,
                   bool                           ignoreErrors);

    virtual ~FileReaderBase() { }

    /// \brief Return the file type string. Used for logging message purposes.
    virtual std::string getFileTypeStr() const { return "Unprocessed"; }

    /// \brief Get file path.
    /// \return File path.
    /// \pre <tt>this->isOpen() == true</tt>
    virtual const boost::filesystem::path& getPath() const { return path_; }

    virtual void openFile(const boost::filesystem::path& path,
                          bool                           ignoreErrors,
                          size_t                         tileIdx) = 0;

    virtual bool isOpen() const = 0;

    virtual bool readEntireFile(common::RawDataBuffer& buffer) = 0;

protected:
   boost::filesystem::path path_;

   bool ignoreErrors_;
}; // end class FileReaderBase


class UnprocessedFile : public FileReaderBase
{
public:
    UnprocessedFile(std::ios_base::openmode openFlags);

    UnprocessedFile(const boost::filesystem::path& path,
                    bool                           ignoreErrors,
                    std::ios_base::openmode        openFlags = std::ios_base::in | std::ios_base::binary);

    virtual ~UnprocessedFile();

    virtual void openFile(const boost::filesystem::path& path,
                          bool                           ignoreErrors,
                          size_t                         tileIdx);

    virtual bool isOpen() const { return !path_.empty(); }

    /// \brief Return the file type string. Used for logging message purposes.
    virtual std::string getFileTypeStr() const { return "Unprocessed"; }

    virtual bool readEntireFile(common::RawDataBuffer& buffer);

protected:

    virtual bool openFileBuf();

    virtual void openFileImpl(const boost::filesystem::path& path,
                              bool                           ignoreErrors);

    virtual bool readBytes(std::vector<char>& buffer,
                           uint32_t bytes);

    virtual void logError(std::streamsize bytesRead,
                          std::streamsize bytesExpected) const;

    /// \brief Filebuf.
    io::FileBufWithReopen fileBuf_;

    size_t fileSize_;
};

// Note: It is unlikely that this synchronization is even necessary, but I'm
// doing it to ensure correctness. In general, 100s of files are read for each
// tile (1 BCL file for each cycle). Read tasks are popped off the queue in
// order. The only way we would need synchronization is if somehow 1 read
// task took long enough for all of the other files to be read and the next
// tile to begin processing.
class SyncFile : public UnprocessedFile
{
public:

    // SyncFileReader must be used when reading the file in order to ensure
    // that the tile number is incremented and that the lock is held.
    class SyncFileReader
    {
    public:
        SyncFileReader(SyncFile& syncFile,
                       size_t tileIdx)
            : syncFile_(syncFile),
              lock_(syncFile.mut_)
        {
            syncFile_.cv_.wait(lock_, [this, tileIdx] { return tileIdx == syncFile_.currentTile_; });
        }

        ~SyncFileReader()
        {
            ++syncFile_.currentTile_;
            lock_.unlock();
            syncFile_.cv_.notify_all();
        }

        const boost::filesystem::path& getPath() const
        {
            return syncFile_.getPath();
        }

        void openFile(const boost::filesystem::path& path,
                      bool                           ignoreErrors)
        {
            syncFile_.openFile(path, ignoreErrors, syncFile_.currentTile_);
        }

        bool readEntireFile(common::RawDataBuffer& buffer)
        {
            buffer.path_ = syncFile_.path_;
            return syncFile_.readEntireFile(buffer);
        }

        bool read(common::RawDataBuffer& buffer,
                  std::streamsize bytes)
        {
            buffer.path_ = syncFile_.path_;
            return syncFile_.read(buffer, bytes);
        }

        bool read(char* buffer,
                  std::streamsize bytes)
        {
            return syncFile_.read(buffer, bytes);
        }

        bool readGzipBlocks(common::RawDataBuffer& buffer,
                            std::streamsize startPos,
                            std::streamsize nextStartPos)
        {
            buffer.path_ = syncFile_.path_;
            return syncFile_.readGzipBlocks(buffer, startPos, nextStartPos);
        }

        bool seek(std::streamsize offset,
                  std::ios_base::seekdir way = std::ios_base::cur)
        {
            return syncFile_.seek(offset, way);
        }

    private:
        SyncFile& syncFile_;
        std::unique_lock<std::mutex> lock_;
    };

    SyncFile(std::ios_base::openmode openFlags);

    SyncFile(const boost::filesystem::path& path,
             bool                           ignoreErrors,
             std::ios_base::openmode        openFlags);

    virtual ~SyncFile();

    /// \brief Return the file type string. Used for logging message purposes.
    virtual std::string getFileTypeStr() const { return "Synchronized"; }

 private:

    bool appendBgzfBlock(std::vector<char> &gzipData);

private:

    virtual bool read(char* buffer,
                      std::streamsize bytes);

    virtual bool read(std::vector<char>& buffer,
                      std::streamsize bytes);

    virtual bool seek(std::streamsize offset,
                      std::ios_base::seekdir way = std::ios_base::cur);

    /// \brief Read enough bgzf gzip blocks to get the tile data.
    /// \param buffer Holds the gzip char data read from file.
    /// \param startPos The file position to the start of the gzip block containing the tile data.
    /// \param nextStartPos The file position to the start of the gzip block containing the next tile data. If this value is zero, read to the end of the file.
    virtual bool readGzipBlocks(std::vector<char>& buffer,
                                std::streamsize startPos,
                                std::streamsize nextStartPos);

    std::mutex mut_;
    std::condition_variable cv_;
    std::atomic<size_t> currentTile_;
    std::vector<char> buffer_;
};


#define CATCH_AND_IGNORE(EXCEPTION_TYPE, MSG_BEGIN)\
catch (EXCEPTION_TYPE& e)\
{\
    if (!ignoreErrors_)\
    {\
        BCL2FASTQ_LOG(common::LogLevel::ERROR_TYPE)\
            << " File: " << getPath().string() << "'." << std::endl;\
\
        BOOST_THROW_EXCEPTION(e);\
    }\
\
    BCL2FASTQ_LOG(common::LogLevel::WARNING)\
            << MSG_BEGIN << getPath().string() << "'." << std::endl\
            <<  e.what() << std::endl;\
}

#define CATCH_AND_IGNORE_ALL(MSG_BEGIN)\
catch (...)\
{\
    if (!ignoreErrors_)\
    {\
        throw;\
    }\
\
    BCL2FASTQ_LOG(common::LogLevel::WARNING)\
            << MSG_BEGIN << getPath().string() << ".'" << std::endl\
            <<  "Unexpected exception." << std::endl;\
}

} // namespace io
} // namespace bcl2fastq

#endif // BCL2FASTQ_IO_SYNC_FILE_HH

