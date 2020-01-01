/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file PositionsFile.hh
 *
 * \brief Declaration of positions file.
 *
 * \author Marek Balint
 * \author Mauricio Varea
 */


#ifndef BCL2FASTQ_DATA_POSITIONSFILE_HH
#define BCL2FASTQ_DATA_POSITIONSFILE_HH


#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/filesystem/path.hpp>

#include "data/FileReader.hh"
#include "common/Types.hh"


namespace bcl2fastq {
namespace data {


/// \brief Positions file.
class PositionsFile : public virtual FileReaderBase
{
public:

    /// \brief Positions file record type definition.
    struct Record
    {
    public:

        /// \brief Cluster coordinate type definition.
        typedef unsigned int ClusterCoordinate;

        /// \brief X-coordinate.
        ClusterCoordinate x_;

        /// \brief y-coordinate.
        ClusterCoordinate y_;
    };

    /// \brief Read records from the file to buffer.
    /// \param targetBuffer Target buffer to read to.
    /// \param targetSize Maximum number of records to be read.
    /// \return Number of records read.
    /// \pre <tt>this->isOpen() == true</tt>
    virtual std::size_t read(
        std::vector<data::PositionsFile::Record>& targetBuffer,
        std::size_t targetSize
    ) = 0;

protected:
    PositionsFile(const common::RawDataBuffer& data,
                  bool ignoreErrors = false
    );
};




/// \brief Factory class to instantiate the appropriate PositionsFile type
class PositionsFileFactory : private boost::noncopyable
{
public:
    /// \brief Return true if the positions file exists.
    /// \param intensitiesDir Path to intensities directory.
    /// \param aggregateTilesFlag All tiles aggregated into single file flag.
    /// \param isPatternedFlowcell True if a patterned flowcell was used.
    /// \param laneNumber Lane number.
    /// \param tileNumber Tile number.
    static bool doesFileExist(const boost::filesystem::path& intensitiesDir,
                              bool                           aggregateTilesFlag,
                              bool                           isPatternedFlowcell,
                              common::LaneNumber             laneNumber,
                              common::TileNumber             tileNumber,
                              size_t&                        headerSize,
                              boost::filesystem::path&       posFilePath);

    static std::streamsize read(const common::RawDataBuffer& inputData,
                                bool ignoreErrors,
                                bool skipHeader,
                                std::vector<data::PositionsFile::Record>& targetBuffer,
                                std::size_t targetSize);

private:
    /// \brief Check if the binary file path exists for the specified parameters.
    /// \param intensitiesDir Path to intensities directory.
    /// \param aggregateTilesFlag All tiles aggregated into single file flag.
    /// \param isPatternedFlowcell True if a patterned flowcell was used.
    /// \param laneNumber Lane number.
    /// \param tileNumber Tile number.
    /// \param fileExtension File extension.
    /// \param positionsFilePath Output param for the file path.
    /// \return true if the file exists.
    static bool binaryFilePathExists(
        const boost::filesystem::path& intensitiesDir,
        bool                           aggregateTilesFlag,
        bool                           isPatternedFlowcell,
        common::LaneNumber             laneNumber,
        common::TileNumber             tileNumber,
        const boost::filesystem::path& fileExtension,
        boost::filesystem::path&       positionsFilePath
    );

    /// \brief Check if the *_pos.txt file path exists for the specified parameters
    /// \param intensitiesDir Path to intensities directory.
    /// \param laneNumber Lane number.
    /// \param tileNumber Tile number.
    /// \param positionsFilePath Output param for the file path.
    /// \return true if the file exists.
    static bool posFilePathExists(
        const boost::filesystem::path& intensitiesDir,
        common::LaneNumber             laneNumber,
        common::TileNumber             tileNumber,
        boost::filesystem::path&       positionsFilePath
    );
};

} // namespace data
} // namespace bcl2fastq


#endif // BCL2FASTQ_DATA_POSITIONSFILE_HH

