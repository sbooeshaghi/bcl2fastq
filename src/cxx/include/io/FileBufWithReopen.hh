/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file FileBufWithReopen.hh
 *
 * \brief Declaration of the file buffer with reopen.
 *
 * \author Roman Petorvski
 * \author Marek Balint
 */


#ifndef BCL2FASTQ_IO_FILEBUFWITHREOPEN_HH
#define BCL2FASTQ_IO_FILEBUFWITHREOPEN_HH


#include <ios>
#include <fstream>

#include <cstdio>

#include <boost/filesystem/path.hpp>
#include <boost/noncopyable.hpp>


namespace bcl2fastq {


namespace io {


/// \brief Same as std::filebuf but has reopen which is useful in parts of the code that are not allowed to do resource allocations.
template<typename CharT, typename Traits = std::char_traits<CharT> >
class BasicFileBufWithReopen : public std::basic_filebuf<CharT, Traits>, private boost::noncopyable
{
public:

    /// \brief Fadvise flags.
    struct FadviseFlags
    {
        enum value_type
        {
            NORMAL = 0,                             ///< No advice (default). Sets readahead window to default size.
            SEQUENTIAL = 1,                         ///< Sequential access. Sets readahead widnow to double the default size.
            RANDOM = 2,                             ///< Random access. Disables readahead window.
            NOREUSE = 4,                            ///< Same as @c WILLNEED until kernel 2.6.18. No-op since then.
            WILLNEED = 8,                           ///< Data will be acessed in near future.
            DONTNEED = 16,                          ///< Data will not be accessed in near future.
            SEQUENTIAL_ONCE = SEQUENTIAL|DONTNEED,  ///< Sequential access with caching.
            SEQUENTIAL_OFTEN = SEQUENTIAL|WILLNEED  ///< Sequential access without caching.
        };
    };

public:

    /// \brief Constructor.
    /// \param mode Access mode to be used.
    explicit BasicFileBufWithReopen(std::ios_base::openmode mode);

private:

    /// \brief Copy assignment operator.
    /// \param rhs Right-hand-side parameter.
    BasicFileBufWithReopen & operator=(const BasicFileBufWithReopen &rhs);
    
public:

    /// \brief (Re)opens specified file.
    /// \param path File to be opened.
    /// \param fadvise Fadvise flags.
    /// \return Pointer to this.
    /// \note Access mode specified during construction of the object
    /// will be used.
    /// \note Call to this method supposedly does not fail due to memory
    /// allocation failure as it does not do any.
    BasicFileBufWithReopen* reopen(const boost::filesystem::path& path, typename FadviseFlags::value_type fadvise = FadviseFlags::NORMAL);

    /// \brief Set access mode.
    /// \param mode Access mode to be set.
    void setMode(std::ios_base::openmode mode);

    /// \brief Get access mode.
    /// \return Access mode.
    std::ios_base::openmode getMode() const;

    /// \brief Flush pending data.
    void flush();

    /// \brief Reset.
    void reset();

private:

    /// \brief Reserve a file handle.
    /// \retval true Success.
    /// \retval false Failure.
    bool reserve();

private:

#ifndef WIN32
    /// \brief Abstraction of implementation dependent external buffer of underlying file stream.
    class File
    {
    public:

        /// \brief File type definition.
        typedef typename BasicFileBufWithReopen<CharT, Traits>::__file_type file_type;

    public:

        /// \brief Constructor.
        /// \param file Implementation dependent buffer of underlying file stream.
        explicit File(file_type &file);

    public:

        /// \brief Get C-file.
        /// \return Extracted stdio's @c FILE structure.
        std::FILE * getCFile();

    private:

        /// \brief Reference to implementation dependent buffer of underlying file stream.
        file_type &file_;
    };
#endif

private:

    /// \brief File access mode.
    const std::ios_base::openmode mode_;

#ifndef WIN32
    /// \brief Implementation dependent buffer of underlying file stream.
    File file_;
#endif
};


/// \brief File buffer with reopen using 1-byte characters.
typedef BasicFileBufWithReopen<char> FileBufWithReopen;


/// \brief Read bytes from file buffer.
/// \param fileBuf File buffer to read from.
/// \param targetBuffer Target buffer to read to.
/// \param targetSize Maximum number of bytes to be read.
/// \return Number of bytes read.
template<typename CharT>
std::streamsize read(BasicFileBufWithReopen<CharT> &fileBuf, char *targetBuffer, std::streamsize targetSize);


} // namespace io


} // namespace bcl2fastq


#include "io/FileBufWithReopen.hpp"


#endif // BCL2FASTQ_IO_FILEBUFWITHREOPEN_HH


