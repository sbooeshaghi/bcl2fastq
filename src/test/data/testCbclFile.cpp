// data/testCbclFile.cpp

#include <vector>
#include <string>
#include <boost/format.hpp>
#include <iostream>
#include <ios>

#include "common/Debug.hh"
#include "common/SystemCompatibility.hh"
#include "common/Options.hh"
#include "config/Bcl2FastqOptions.hh"
#include "io/SyncFile.hh"
#include "common/Logger.hh"
#include "data/CbclFile.hh"
#include "gtest/gtest.h"
#include "../bgzfData.hh"

using namespace bcl2fastq;

namespace {

    const bool verbose = true;

    typedef uint32_t header_value_t;
    typedef std::vector<header_value_t> header_value_vector_t;
    typedef uint16_t header_version_t;
    typedef uint8_t header_byte_t;
    typedef uint32_t header_value_t;
    typedef std::vector<header_value_t> header_value_vector_t;
    typedef char stream_byte_t;
    typedef std::vector<stream_byte_t> header_t;

    template<class T>
    bool readBinary(std::ifstream& in, T& buffer)
    {
        in.read(reinterpret_cast<char *>(&buffer), sizeof(T));
        return (in) ? true : false;
    }

    template<class T>
    void writeBinary(std::ofstream& out, T& buffer)
    {
       for (int i = 0; i < sizeof(T); ++i) {
          out << (unsigned char)((buffer >> (8*i)) & 0xff);
       }
    }

    class CbclFileTester {
    public:

        // This MUST match the header code in bcl2fastq...
        // This is copied from bcl2fastq CbclFile.hh/cpp.

        header_version_t m_version;
        header_value_t m_header_size;
        header_byte_t m_number_of_bits_per_basecall;
        header_byte_t m_number_of_bits_per_qscore;
        header_value_t m_number_of_qscore_bins = 0;
        header_value_vector_t m_qscore_bins;
        header_value_vector_t m_remapped_qscores;
        header_value_t m_number_of_tiles;
        header_value_vector_t m_tile_ids;
        header_value_vector_t m_number_of_clusters;
        header_value_vector_t m_uncompressed_block_size;
        header_value_vector_t m_compressed_block_size;
        header_byte_t m_non_pf_filtered;
        std::vector<size_t> m_start_positions;

        std::ifstream readCbclHeader(const std::string &filename, bool verbose) {

            std::ifstream in(filename.c_str(), std::ios_base::in | std::ios_base::binary);
            readHeader(in, verbose);
        }

        void readHeader(std::ifstream& in, bool verbose)
        {
            in.seekg(0,in.beg);
            readBinary(in, m_version);

            readBinary(in, m_header_size);
            readBinary(in, m_number_of_bits_per_basecall);
            readBinary(in, m_number_of_bits_per_qscore);
            readBinary(in, m_number_of_qscore_bins);

            if (verbose) {
                std::cout << "Version: " << m_version << std::endl;
                std::cout << "Header size: " << m_header_size << std::endl;
                std::cout << "bpb: " << (int)m_number_of_bits_per_basecall << std::endl;
                std::cout << "bpq: " << (int)m_number_of_bits_per_qscore << std::endl;
                std::cout << "qscore bins: " << m_number_of_qscore_bins << std::endl;
            }

            if (m_number_of_qscore_bins > 0)
            {
                m_qscore_bins.resize(m_number_of_qscore_bins);
                m_remapped_qscores.resize(m_number_of_qscore_bins);

                for (header_value_t i = 0; i < m_number_of_qscore_bins; ++i)
                {
                    readBinary(in, m_qscore_bins[i]);
                    readBinary(in, m_remapped_qscores[i]);
                    if (verbose) {
                        std::cout << i << ": "  << m_qscore_bins[i] << "; " << m_remapped_qscores[i] << std::endl;
                    }
                }
            }

            readBinary(in, m_number_of_tiles);
            if (verbose)
                std::cout << "Tile count: " << m_number_of_tiles << std::endl;

            m_tile_ids.resize(m_number_of_tiles);
            m_number_of_clusters.resize(m_number_of_tiles);
            m_uncompressed_block_size.resize(m_number_of_tiles);
            m_compressed_block_size.resize(m_number_of_tiles);
            for (header_value_t i = 0; i < m_number_of_tiles; ++i)
            {
                readBinary(in, m_tile_ids[i]);
                readBinary(in, m_number_of_clusters[i]);
                readBinary(in, m_uncompressed_block_size[i]);
                readBinary(in, m_compressed_block_size[i]);

                if (verbose)
                    std::cout << i << ": "  << m_tile_ids[i] << "; " << m_number_of_clusters[i] << "; " << m_uncompressed_block_size[i] << std::endl;
            }


            readBinary(in, m_non_pf_filtered);
            if (verbose)
                std::cout << "non-pf: " << m_non_pf_filtered << std::endl;

            //m_start_positions.resize(m_number_of_tiles);
            //m_start_positions[0] = m_header_size;
            //for (size_t i = 1; i < m_number_of_tiles; ++i)
            //{
            //    m_start_positions[i] = m_start_positions[i-1] + m_compressed_block_size[i-1];
            //}
        }

        void writeHeader(std::ofstream& out)
        {
           writeBinary(out, m_version);
           writeBinary(out, m_header_size);
           writeBinary(out, m_number_of_bits_per_basecall);
           writeBinary(out, m_number_of_bits_per_qscore);
           writeBinary(out, m_number_of_qscore_bins);

           if (m_number_of_qscore_bins > 0)
           {
               for (header_value_t i = 0; i < m_number_of_qscore_bins; ++i)
               {
                   writeBinary(out, m_qscore_bins[i]);
                   writeBinary(out, m_remapped_qscores[i]);
               }
           }

           writeBinary(out, m_number_of_tiles);

           for (header_value_t i = 0; i < m_number_of_tiles; ++i)
           {
               writeBinary(out, m_tile_ids[i]);
               writeBinary(out, m_number_of_clusters[i]);
               writeBinary(out, m_uncompressed_block_size[i]);
               writeBinary(out, m_compressed_block_size[i]);
           }

           writeBinary(out, m_non_pf_filtered);

           ////wrong
           //writeBinary(out, m_version);
        }

    }; // end class CbclFileTester

} // end unnamed namespace


//=========================================
//
// TESTING A
//
//=========================================

TEST(testBasic, ValidateTypes) {

    EXPECT_EQ(8, sizeof(uint64_t));
    EXPECT_EQ(4, sizeof(uint32_t));
    EXPECT_EQ(4, sizeof(unsigned));
    EXPECT_EQ(8, sizeof(size_t));
    EXPECT_EQ(1, sizeof(char));
}

//TEST(testCbclFile, CheckHeader) {
//
//    const std::string filename = "/illumina/scratch/OncoBU/users/lteng/NovaSeq_TST170/Data/Intensities/BaseCalls/L001/C1.1/L001_1.cbcl";
//    std::ifstream in(filename.c_str(), std::ios_base::in | std::ios_base::binary);
//
//    CbclFileTester cbcl;
//    cbcl.readHeader(in, verbose);
//
//    EXPECT_EQ(1, (int)cbcl.m_version);
//}

//TEST(CbclFileture/btof1110Util, Truncate) {
//
//    const std::string filename = "/illumina/scratch/Voyager/170517_AP2-12_0173_A075kdnhvy/Data/Intensities/BaseCalls/L001/C6.1/L001_1.cbcl";
//    std::ifstream in(filename.c_str(), std::ios_base::in | std::ios_base::binary);
//
//    CbclFileTester cbcl;
//    cbcl.readHeader(in, verbose);
//
//    const std::string outFilename = "trunc_L001_1.cbcl";
//    std::ofstream out(outFilename.c_str(), std::ios_base::out | std::ios_base::binary);
//    cbcl.writeHeader(out);
//
//}

