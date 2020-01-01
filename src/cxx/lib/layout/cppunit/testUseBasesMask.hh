/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file testUseBasesMask.hh
 *
 * \brief UseBasesMask cppunit test declarations.
 *
 * \author Aaron Day
 */


#ifndef BCL2FASTQ_LAYOUT_TEST_USE_BASES_MASK_HH
#define BCL2FASTQ_LAYOUT_TEST_USE_BASES_MASK_HH


#include <string>

#include "layout/ReadMetadata.hh"

#include <cppunit/extensions/HelperMacros.h>


/// \brief Test suite for UseBasesMask.
class TestUseBasesMask : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE(TestUseBasesMask);
    CPPUNIT_TEST(testAllButOne);
    CPPUNIT_TEST(testMaskingR1);
    CPPUNIT_TEST(testMaskingR2);
    CPPUNIT_TEST(testMixedDataIndex);
    CPPUNIT_TEST(testExpansion);
    CPPUNIT_TEST(testMixedDataIndexSingleRead);
    CPPUNIT_TEST(testIsomorphic);
    CPPUNIT_TEST(testMismatchingCyclesInRead);
    CPPUNIT_TEST(testZero);
    CPPUNIT_TEST(testMissingRead);
    CPPUNIT_TEST(testMissingAsteriskNone);
    CPPUNIT_TEST_SUITE_END();

private:
    std::vector< bcl2fastq::layout::ReadMetadataContainer > readMetadata_;

public:
    TestUseBasesMask() : readMetadata_(0) {}
    void setUp();
    void tearDown();
    void testAllButOne();
    void testMaskingR1();
    void testMaskingR2();
    void testMixedDataIndex();
    void testExpansion();
    void testMixedDataIndexSingleRead();
    void testIsomorphic();
    void testMismatchingCyclesInRead();
    void testZero();
    void testMissingRead();
    void testMissingAsteriskNone();
};


#endif // BCL2FASTQ_LAYOUT_TEST_USE_BASES_MASK_HH


