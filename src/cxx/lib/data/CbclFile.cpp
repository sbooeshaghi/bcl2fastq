/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA^M
 * and certain third party copyright/licenses, and any user of this^M
 * source file is bound by the terms therein.^M
 *
 * \file CbclFile.cpp
 *
 * \brief Implementation of CBCL file.
 *
 * \author Eunho Noh
 * \author Aaron Day
 */


#include "data/CbclFile.hh"
#include "layout/TileInfo.hh"

#include <boost/filesystem/operations.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/device/array.hpp>
#include <boost/lexical_cast.hpp>




#include <boost/iostreams/filter/gzip.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/filtering_stream.hpp>

namespace bcl2fastq {
namespace data {



CbclFileReader::Header::Header()
: m_version(0)
, m_header_size(0)
, m_number_of_bits_per_basecall(0)
, m_number_of_bits_per_qscore(0)
, m_number_of_qscore_bins(0)
, m_qscore_bins()
, m_remapped_qscores()
, m_number_of_tiles(0)
, m_tile_ids()
, m_number_of_clusters()
, m_uncompressed_block_size()
, m_compressed_block_size()
, m_non_pf_filtered(0)
, m_start_positions()
{
}

bool CbclFileReader::Header::load(io::SyncFile::SyncFileReader &cbclFile, bool ignoreErrors) 
{
    try
    {
        // This can throw an exception if the header is corrupt.
        readHeader(cbclFile);
    }
    catch(...)
    {
        if (ignoreErrors)
        {
            static const std::string msg = "Error reading CBCL file header. All data for this cycle in unavailable and will be skipped. ";
            BCL2FASTQ_LOG(common::LogLevel::WARNING) << msg << std::endl;
            return false;
        }
        else {
            throw;
        }
    }

    const std::string errorMsg = validateHeaderInfo();
    if (!errorMsg.empty()) 
    {
        if (ignoreErrors)
        {
            BCL2FASTQ_LOG(common::LogLevel::WARNING) << errorMsg << "(" << cbclFile.getPath() << ") " << std::endl;
            return false;
        }
        else 
        {
            BOOST_THROW_EXCEPTION(std::runtime_error(errorMsg));
        }
    }

    return true;
}

void CbclFileReader::Header::readHeader(io::SyncFile::SyncFileReader& in)
{
    read_binary(in, m_version);
    read_binary(in, m_header_size);
    read_binary(in, m_number_of_bits_per_basecall);
    read_binary(in, m_number_of_bits_per_qscore);
    read_binary(in, m_number_of_qscore_bins);

    if (m_number_of_qscore_bins > 0)
    {
        m_qscore_bins.resize(m_number_of_qscore_bins);
        m_remapped_qscores.resize(m_number_of_qscore_bins);

        for (header_value_t i = 0; i < m_number_of_qscore_bins; ++i)
        {
            read_binary(in, m_qscore_bins[i]);
            read_binary(in, m_remapped_qscores[i]);
        }
    }

    read_binary(in, m_number_of_tiles);

    if (m_number_of_tiles > 0) 
    {
        m_tile_ids.resize(m_number_of_tiles);
        m_number_of_clusters.resize(m_number_of_tiles);
        m_uncompressed_block_size.resize(m_number_of_tiles);
        m_compressed_block_size.resize(m_number_of_tiles);
        for (header_value_t i = 0; i < m_number_of_tiles; ++i)
        {
            read_binary(in, m_tile_ids[i]);
            read_binary(in, m_number_of_clusters[i]);
            read_binary(in, m_uncompressed_block_size[i]);
            read_binary(in, m_compressed_block_size[i]);
        }
    }

    read_binary(in, m_non_pf_filtered);

    if (m_number_of_tiles > 0) 
    {
        m_start_positions.resize(m_number_of_tiles);
        m_start_positions[0] = m_header_size;
        for (size_t i = 1; i < m_number_of_tiles; ++i)
        {
            m_start_positions[i] = m_start_positions[i-1] + m_compressed_block_size[i-1];
        }
    }
}

void CbclFileReader::Header::reset()
{
    m_version = 0;
    m_header_size = 0;
    m_number_of_bits_per_basecall = 0;
    m_number_of_bits_per_qscore = 0;
    m_number_of_qscore_bins = 0;
    m_qscore_bins.resize(0);
    m_remapped_qscores.resize(0);
    m_number_of_tiles = 0;
    m_tile_ids.resize(0);
    m_number_of_clusters.resize(0);
    m_uncompressed_block_size.resize(0);
    m_compressed_block_size.resize(0);
    m_non_pf_filtered = 0;
    m_start_positions.resize(0);
}

CbclFileReader::Header::header_value_t CbclFileReader::Header::calculateHeaderSize(size_t number_of_qscore_bins, size_t number_of_tiles) const
{
    return sizeof(header_version_t) +                           //version number
           sizeof(header_value_t) +                             // header size
           sizeof(header_byte_t) +                              // number of bits per basecall
           sizeof(header_byte_t) +                              // number of bits per q-score
           sizeof(header_value_t) +                             // number of q-score bins
           number_of_qscore_bins * sizeof(header_value_t) +     // q-score bins
           number_of_qscore_bins * sizeof(header_value_t) +     // q-scores remapped
           sizeof(header_byte_t) +                              // non-pf clusters filtered/un-filtered
           sizeof(header_value_t) +                             // number of tiles
           number_of_tiles * sizeof(header_value_t) +           // tile ids
           number_of_tiles * sizeof(header_value_t) +           // number of clusters per tile
           number_of_tiles * sizeof(header_value_t) +           // uncompressed size
           number_of_tiles * sizeof(header_value_t);            // compressed size
}

std::string CbclFileReader::Header::validateHeaderInfo() const
{
    static std::string badFileMsg = "A corrupt CBCL file was identified and cannot be loaded. ";;
    if (m_version != 1)
    {
        return badFileMsg; // + "Unsupported CBCL file version. ";
    }
    else if(m_header_size != calculateHeaderSize(m_number_of_qscore_bins, m_number_of_tiles))
    {
        return badFileMsg; // + "The header size information does not match the expected size.";
    }
    else if(m_number_of_bits_per_basecall+m_number_of_bits_per_qscore != NUM_BITS_BCL_STORAGE)
    {
        return badFileMsg; // + "Unsupported basecall/qscore size";
    }
    else if(m_number_of_qscore_bins != m_qscore_bins.size() || m_number_of_qscore_bins != m_remapped_qscores.size() )
    {
        return badFileMsg; // + "The q-score information is invalid.";
    }
    else if(m_number_of_tiles != m_tile_ids.size() || m_number_of_tiles != m_number_of_clusters.size() ||
       m_number_of_tiles != m_uncompressed_block_size.size() || m_number_of_tiles != m_compressed_block_size.size())
    {
        return badFileMsg; // + "The tile information is invalid.";
    }
    for(size_t tile = 0; tile < m_number_of_tiles; ++tile)
    {
        if(m_number_of_clusters[tile] > m_uncompressed_block_size[tile] * BASECALLS_PER_BYTE)
        {
            return badFileMsg; // + "The tile information is invalid.";
        }
    }

    static const std::string noError;
    return noError;
}


CbclFileReader::CbclFileReader(const::boost::filesystem::path& inputDir,
                               const layout::LaneInfo&         laneInfo,
                               common::CycleNumber             cycleNumber,
                               size_t                          cycleIndex,
                               bool                            ignoreMissingBcls,
                               const common::TileFileMap&      tileFileMap,
                               std::shared_ptr<io::SyncFile> bclFile)
: BclFileReaderT<io::SyncFile>(inputDir,
                               laneInfo,
                               cycleNumber,
                               cycleIndex,
                               ignoreMissingBcls,
                               bclFile)
, tileFileMap_(tileFileMap)
, header_()
{
}

bool CbclFileReader::doesFileExist(const boost::filesystem::path& inputDir,
                                   common::LaneNumber             laneNumber,
                                   common::CycleNumber            cycleNumber,
                                   common::TileNumber             tileNumber,
                                   const common::TileFileMap&     tileFileMap,
                                   bool                           ignoreErrors)
{
    boost::filesystem::path fileName = getFileName(inputDir,
                                                   laneNumber,
                                                   cycleNumber,
                                                   tileNumber,
                                                   tileFileMap,
                                                   ignoreErrors);

    return boost::filesystem::exists(fileName); 
}

void CbclFileReader::markCycleDataInvalid(PerCycleData &cycleData) 
{
    // reset data
    cycleData.bcls_.resize(0);
    cycleData.uncompressedBlockSize_ = 0;

    if (header_) { header_->reset(); }
}

bool CbclFileReader::read(RawBclBufferGroup& outputBuffer)
{
    for (auto& outputBufferForTile : outputBuffer)
    {
        io::SyncFile::SyncFileReader fileReader(*bclFile_,
                                                std::distance(laneInfo_.getTileInfos().begin(),
                                                              outputBufferForTile.tileInfo_));

        auto& cycleData = outputBufferForTile.cycleData_[cycleIndex_];

        boost::filesystem::path fileName = getFileName(inputDir_,
                                                       laneInfo_.getNumber(),
                                                       cycleNumber_,
                                                       outputBufferForTile.tileInfo_->getNumber(),
                                                       tileFileMap_,
                                                       ignoreMissingBcls_);

        if (ignoreMissingBcls_ && !boost::filesystem::exists(fileName))
        {
            // Method getFileName issues a warning.
            cycleData.bcls_.path_ = fileName;
            markCycleDataInvalid(cycleData);
            continue;
        } 
        // else let it fail below

        if (header_ && header_->getHeaderSize() == 0) {
            continue;
        }

        if (!header_ || !header_->containsTile(outputBufferForTile.tileInfo_->getNumber()))
        {
            // header_ = new
            fileReader.openFile(fileName,
                                ignoreMissingBcls_);

            header_ = std::unique_ptr<Header>(new Header());
            if (!header_->load(fileReader, ignoreMissingBcls_)) 
            {
                cycleData.bcls_.path_ = fileName;
                markCycleDataInvalid(cycleData);
                continue;
            }
        }

        cycleData.includeNonPf_ = !header_->isNonPfFiltered();

        if (header_->isQscoreDeflated())
        {
            cycleData.numBitsPerQscore_ = header_->getNumberOfBitsPerQscore();
            cycleData.remappedQscores_ = header_->getRemappedQscores();
        }

        auto tileIndex = header_->getTileIndex(outputBufferForTile.tileInfo_->getNumber());
        const size_t startPosition = header_->getBlockStartPosition(tileIndex);
        const uint32_t compressedBlockSize = header_->getCompressedBlockSize(tileIndex);

        cycleData.uncompressedBlockSize_ = header_->getUncompressedBlockSize(tileIndex);

        fileReader.seek(startPosition,
                        std::ios_base::beg);

        cycleData.bcls_.resize(compressedBlockSize);
        fileReader.read(cycleData.bcls_,
                        compressedBlockSize);
    }

    return true;
}

boost::filesystem::path CbclFileReader::getFileName(const boost::filesystem::path& inputDir,
                                              common::LaneNumber             laneNumber,
                                              common::CycleNumber            cycleNumber,
                                              common::TileNumber             tileNumber,
                                              const common::TileFileMap&     tileFileMap,
                                              bool                           ignoreErrors)
{
    boost::filesystem::path cbclPath;
    auto pos = tileFileMap.at(laneNumber).find(tileNumber);
    if (pos != tileFileMap.at(laneNumber).end())
    {
        cbclPath =
            inputDir /
            boost::filesystem::path((boost::format("L%03d") % laneNumber).str()) /
            boost::filesystem::path((boost::format("C%d.1") % cycleNumber).str()) /
            pos->second;
    }

    if (pos == tileFileMap.at(laneNumber).end() || !boost::filesystem::exists(cbclPath))
    {
        std::string errStr("Unable to find CBCL file: '" + cbclPath.string() + "'");
        if (ignoreErrors)
        {
            BCL2FASTQ_LOG(common::LogLevel::WARNING) << errStr << std::endl;
        }
        else
        {
            BOOST_THROW_EXCEPTION(common::IoError(ENOENT, errStr));
        }
    }

    return cbclPath;
}

} // namespace data
} // namespace bcl2fastq

