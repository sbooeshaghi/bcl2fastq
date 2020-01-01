/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file testGzipCompressor.hh
 *
 * \brief GzipCompressor cppunit test declarations.
 *
 * \author Aaron Day
 */


#ifndef BCL2FASTQ_IO_TEST_GZIP_COMPRESSOR_HH
#define BCL2FASTQ_IO_TEST_GZIP_COMPRESSOR_HH


#include <cppunit/extensions/HelperMacros.h>

#include <string>
#include <vector>


/// \brief Test suite for GzipCompressor.
class TestGzipCompressor : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE(TestGzipCompressor);
    CPPUNIT_TEST(testBgzf);
    CPPUNIT_TEST(testSingleBlock);
    CPPUNIT_TEST_SUITE_END();

public:
    TestGzipCompressor() {}

    void setUp();
    void tearDown();
    void testBgzf();
    void testSingleBlock();

private:
    void writeCompressedBlock(const std::string& uncompressedData,
                              bool               useBgzf,
                              std::vector<char>& compressedData);

    void decompressData(const std::vector<char>::const_iterator& compressedData,
                        size_t                                   compressedSize,
                        std::vector<char>&                       decompressedData);

    void writeAllCompressedBlocks(std::vector<std::string>&         uncompressedData,
                                  bool                              useBgzf,
                                  std::vector< std::vector<char> >& compressedData,
                                  std::vector<char>&                allCompressedBlocks);

    void createUncompressedData(std::vector<std::string>& uncompressedData,
                                std::string&              fullString);

    void locateAndDecompressBlocks(const std::vector<std::string>&         uncompressedData,
                                   const std::vector< std::vector<char> >& compressedData,
                                   const std::vector<char>&                allCompressedBlocks);
};


#endif // BCL2FASTQ_IO_TEST_GZIP_COMPRESSOR_HH

