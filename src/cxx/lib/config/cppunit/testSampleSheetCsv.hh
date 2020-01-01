/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file testSampleSheetCsv.hh
 *
 * \brief Sample sheet helper cppunit test declarations.
 *
 * \author Marek Balint
 */


#ifndef BCL2FASTQ_LAYOUT_TEST_SAMPLESHEETCSV_HH
#define BCL2FASTQ_LAYOUT_TEST_SAMPLESHEETCSV_HH


#include <string>

#include <cppunit/extensions/HelperMacros.h>

#include "config/SampleSheetCsv.hh"


/// \brief Test suite for SampleSheetCsv.
class TestSampleSheetCsv : public CppUnit::TestFixture
{
private:

    CPPUNIT_TEST_SUITE(TestSampleSheetCsv);
    CPPUNIT_TEST(testNova);
    CPPUNIT_TEST(testInvalidCharacters);
    CPPUNIT_TEST(testDuplicateBarcode);
    CPPUNIT_TEST_SUITE_END();

public:

    void setUp();
    void tearDown();
    void testNova();
    void testInvalidCharacters();
    void testDuplicateBarcode();
};


#endif // #ifndef BCL2FASTQ_LAYOUT_TEST_SAMPLESHEETCSV_HH


