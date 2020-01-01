/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file CbclFile.hh
 *
 * \brief Declaration of CBCL file.
 *
 * \author Aaron Day
 * \author Eunho Noh
 */


#ifndef BCL2FASTQ_DATA_CBCLFILE_HH
#define BCL2FASTQ_DATA_CBCLFILE_HH

#include "data/BclFileReader.hh"
#include "data/BclFile.hh"
#include "data/BclBuffer.hh"

#include <boost/noncopyable.hpp>
#include <istream>

namespace bcl2fastq {

namespace layout {
class TileInfo;
}

namespace data {

class CbclFileReader : public BclFileReaderT<io::SyncFile>
{
public:
    CbclFileReader(const::boost::filesystem::path& inputDir,
                   const layout::LaneInfo&         laneInfo,
                   common::CycleNumber             cycleNumber,
                   size_t                          cycleIndex,
                   bool                            ignoreMissingBcls,
                   const common::TileFileMap&      tileFileMap,
                   std::shared_ptr<io::SyncFile>   bclFile);

    virtual ~CbclFileReader() { }

    virtual bool read(RawBclBufferGroup& outputBuffer);

    static bool doesFileExist(const boost::filesystem::path& inputDir,
                              common::LaneNumber             laneNumber,
                              common::CycleNumber            cycleNumber,
                              common::TileNumber             tileNumber,
                              const common::TileFileMap&     tileFileMap,
                              bool                           ignoreErrors);

    static boost::filesystem::path getFileName(const boost::filesystem::path& inputDir,
                                               common::LaneNumber             laneNumber,
                                               common::CycleNumber            cycleNumber,
                                               common::TileNumber             tileNumber,
                                               const common::TileFileMap&     tileFileMap,
                                               bool                           ignoreErrors);

class Header : private boost::noncopyable
{
private:
    typedef uint16_t header_version_t;
    typedef uint8_t header_byte_t;
    typedef uint32_t header_value_t;
    typedef std::vector<header_value_t> header_value_vector_t;
    typedef char stream_byte_t;
    typedef std::vector<stream_byte_t> header_t;

public:
    Header();

    header_byte_t getNumberOfBitsPerBasecall() const { return m_number_of_bits_per_basecall; }
    header_byte_t getNumberOfBitsPerQscore() const { return m_number_of_bits_per_qscore; }
    header_value_t getNumberOfTiles() const { return m_number_of_tiles; }
    const header_value_vector_t& getTileIds() const { return m_tile_ids; }

    // Used for invalid cbcl file.
    void reset();
    bool load(io::SyncFile::SyncFileReader& cbclFile, bool ignoreErrors);

    /** Returns the index of the tile number of interest
        *
        * @param tile_id the tile to find
         */
    size_t getTileIndex(size_t tile_id) const
    {
        auto pos = std::find(m_tile_ids.begin(), m_tile_ids.end(), tile_id);
        if(pos == m_tile_ids.end())
        {
            throw std::out_of_range("Invalid tile id");
        }
        return pos - m_tile_ids.begin();
    }

    bool containsTile(size_t tile_id) const
    {
        return std::find(m_tile_ids.begin(), m_tile_ids.end(), tile_id) != m_tile_ids.end();
    }

    /** Returns the tile number
        *
        * @param tile_index the index to the m_tileNumbers vector
         */
    size_t getTileId(size_t tile_index) const
    {
        if(tile_index >= static_cast<size_t>(m_tile_ids.size()))
        {
            throw std::out_of_range("Invalid tile index");
        }
        return m_tile_ids[tile_index];
    }

    const header_value_vector_t& getRemappedQscores() const { return m_remapped_qscores; }
    // header_value_t getRemappedQscore(header_value_t q_score_bin) const { return m_remapped_qscores[q_score_bin]; }

    size_t getBlockStartPosition(size_t tile_index) const { return m_start_positions[tile_index]; }
    size_t getNumberOfClusters(size_t tile_index) const { return m_number_of_clusters[tile_index]; }
    uint32_t getUncompressedBlockSize(size_t tile_index) const { return m_uncompressed_block_size[tile_index]; }
    size_t getCompressedBlockSize(size_t tile_index) const { return m_compressed_block_size[tile_index]; }
    uint8_t isNonPfFiltered() const { return m_non_pf_filtered; }
    bool isQscoreDeflated() const { return m_number_of_qscore_bins > 0; }

    header_value_t getHeaderSize() { return m_header_size; }

private:
    void readHeader(io::SyncFile::SyncFileReader& in);
    header_value_t calculateHeaderSize(size_t number_of_qscore_bins, size_t number_of_tiles) const;
    std::string validateHeaderInfo() const;

        /** Read an array of data from an input stream
             *
             * @param in input stream
             * @param buffer array of data
             * @param n number of elements in the array
             */
    template<class T>
    void read_binary(io::SyncFile::SyncFileReader& in, T* buffer, const size_t n)
    {
        bool validRead = in.read(reinterpret_cast<char *>(buffer), n*sizeof(T));
        if (!validRead) {
            BOOST_THROW_EXCEPTION(std::ios_base::failure("Corrupt CBCL header cannot be loaded. "));
        }

    }
    /** Read an value of binary data from an input stream
     *
     * @param in input stream
     * @param buffer value
     */
    template<class T>
    void read_binary(io::SyncFile::SyncFileReader& in, T& buffer)
    {
        read_binary(in, &buffer, 1);
    }

    const size_t BITS_PER_BYTE = 8;
    const size_t NUM_BITS_BCL_STORAGE = 4;
    const size_t BASECALLS_PER_BYTE = BITS_PER_BYTE / NUM_BITS_BCL_STORAGE;

    header_version_t m_version;
    header_value_t m_header_size;
    header_byte_t m_number_of_bits_per_basecall;
    header_byte_t m_number_of_bits_per_qscore;
    header_value_t m_number_of_qscore_bins;
    header_value_vector_t m_qscore_bins;
    header_value_vector_t m_remapped_qscores;
    header_value_t m_number_of_tiles;
    header_value_vector_t m_tile_ids;
    header_value_vector_t m_number_of_clusters;
    header_value_vector_t m_uncompressed_block_size;
    header_value_vector_t m_compressed_block_size;
    header_byte_t m_non_pf_filtered;
    std::vector<size_t> m_start_positions;
}; // end class Header

    void markCycleDataInvalid(PerCycleData &cycleData);

    const common::TileFileMap& tileFileMap_;
    std::unique_ptr<Header> header_;
};

} // namespace data
} // namespace bcl2fastq

#endif // BCL2FASTQ_DATA_CBCLFILE_HH

