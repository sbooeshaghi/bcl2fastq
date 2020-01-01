/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file testGzipCompressor.cpp
 *
 * \brief GzipCompressor cppunit test declarations.
 *
 * \author Aaron Day
 */


#include "RegistryName.hh"
#include "testGzipCompressor.hh"

#include "io/GzipCompressor.hh"

#include <boost/foreach.hpp>

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(TestGzipCompressor, registryName("GzipCompressor"));

void TestGzipCompressor::setUp()
{
}

void TestGzipCompressor::tearDown()
{
}

void TestGzipCompressor::writeCompressedBlock(const std::string& uncompressedData,
                                              bool               useBgzf,
                                              std::vector<char>& compressedData)
{
    bcl2fastq::io::GzipCompressor(compressedData, useBgzf).write(uncompressedData.c_str(), uncompressedData.size());
}

void TestGzipCompressor::decompressData(const std::vector<char>::const_iterator& compressedData,
                                        size_t                                   compressedSize,
                                        std::vector<char>&                       decompressedData)
{
    boost::iostreams::filtering_ostream os;
    os.push(boost::iostreams::gzip_decompressor());
    os.push(boost::iostreams::back_inserter(decompressedData));
    boost::iostreams::write(os, &(*compressedData), compressedSize);
}

void TestGzipCompressor::writeAllCompressedBlocks(std::vector<std::string>&         uncompressedData,
                                                  bool                              useBgzf,
                                                  std::vector< std::vector<char> >& compressedData,
                                                  std::vector<char>&                allCompressedBlocks)
{
    std::string fullString;
    BOOST_FOREACH(const std::string& uncompressedString, std::make_pair(uncompressedData.begin(), uncompressedData.end()))
    {
        compressedData.push_back(std::vector<char>());
        writeCompressedBlock(uncompressedString,
                             useBgzf,
                             compressedData.back());

        allCompressedBlocks.insert(allCompressedBlocks.end(),
                                   compressedData.back().begin(),
                                   compressedData.back().end());
    }
}

void TestGzipCompressor::createUncompressedData(std::vector<std::string>& uncompressedData,
                                                std::string&              fullString)
{
    uncompressedData.push_back("Hello. This is a test. ");
    uncompressedData.push_back("Here is a second message to compress.                                                               " 
                               "                                                                                                    "
                               "                                                                       Long string                 x" );
    uncompressedData.push_back("And one more.");

    BOOST_FOREACH(const std::string& uncompressedString, std::make_pair(uncompressedData.begin(), uncompressedData.end()))
    {
        fullString += uncompressedString;
    }
}

void TestGzipCompressor::locateAndDecompressBlocks(const std::vector<std::string>&         uncompressedData,
                                                   const std::vector< std::vector<char> >& compressedData,
                                                   const std::vector<char>&                allCompressedBlocks)
{
    int offset = 0;
    for (size_t i = 0; i < uncompressedData.size(); ++i)
    {
        // 4th byte should be the flags field.
        CPPUNIT_ASSERT(allCompressedBlocks.at(offset + 3) & 0x04);

        // 11th and 12th bytes are the length of the extra field.
        CPPUNIT_ASSERT_EQUAL(static_cast<int>(allCompressedBlocks.at(offset + 10)), 4 + (int)bcl2fastq::io::GzipCompressor::SUBFIELD_LENGTH);
        CPPUNIT_ASSERT_EQUAL(static_cast<int>(allCompressedBlocks.at(offset + 11)), 0);

        // 13th and 14th bytes are the subfield identifier.
        CPPUNIT_ASSERT_EQUAL(static_cast<int>(allCompressedBlocks.at(offset + 12)), 66);
        CPPUNIT_ASSERT_EQUAL(static_cast<int>(allCompressedBlocks.at(offset + 13)), 67);

        // 15th and 16th bytes are the length of the subfield data.
        // This is completely redundant with 11th and 12th bytes, but that's the expected format.
        CPPUNIT_ASSERT_EQUAL(static_cast<int>(allCompressedBlocks.at(offset + 14)), (int)bcl2fastq::io::GzipCompressor::SUBFIELD_LENGTH);
        CPPUNIT_ASSERT_EQUAL(static_cast<int>(allCompressedBlocks.at(offset + 15)), 0);

        // The next 8 bytes give the size of the compressed block.
        size_t compressedSize = 0;
        for (size_t j = (16 + bcl2fastq::io::GzipCompressor::SUBFIELD_LENGTH - 1); j >= 16; --j)
        {
            compressedSize <<= 8;
            compressedSize |= allCompressedBlocks.at(offset + j);
        }

        // bgzf uses size-1, so we will too.
        compressedSize += 1;
        CPPUNIT_ASSERT_EQUAL(compressedSize, compressedData[i].size());

        // Decompress a single block.
        std::vector<char> decompressedData;
        decompressData(allCompressedBlocks.begin()+offset,
                       compressedSize,
                       decompressedData);
        CPPUNIT_ASSERT_EQUAL(std::string(decompressedData.begin(), decompressedData.end()), uncompressedData[i]);

        offset += compressedSize;
    }
}

void TestGzipCompressor::testBgzf()
{
    std::vector<std::string> uncompressedData;
    std::string fullString;
    createUncompressedData(uncompressedData,
                           fullString);

    std::vector< std::vector<char> > compressedData;
    std::vector<char> allCompressedBlocks;
    writeAllCompressedBlocks(uncompressedData,
                             true,
                             compressedData,
                             allCompressedBlocks);

    // Try decompressing all blocks at once
    std::vector<char> decompressedAll;
    decompressData(allCompressedBlocks.begin(),
                   allCompressedBlocks.size(),
                   decompressedAll);
    CPPUNIT_ASSERT_EQUAL(std::string(decompressedAll.begin(), decompressedAll.end()), fullString);

    // Test our ability to locate individual compressed blocks
    locateAndDecompressBlocks(uncompressedData,
                              compressedData,
                              allCompressedBlocks);
}

void TestGzipCompressor::testSingleBlock()
{
    std::vector<std::string> uncompressedData;
    std::string fullString;
    createUncompressedData(uncompressedData,
                           fullString);

    std::vector< std::vector<char> > compressedData;
    std::vector<char> allCompressedBlocks;
    writeAllCompressedBlocks(uncompressedData,
                             false,
                             compressedData,
                             allCompressedBlocks);

    // Try decompressing all blocks at once
    std::vector<char> decompressedAll;
    decompressData(allCompressedBlocks.begin(),
                   allCompressedBlocks.size(),
                   decompressedAll);
    CPPUNIT_ASSERT_EQUAL(std::string(decompressedAll.begin(), decompressedAll.end()), fullString);
}
