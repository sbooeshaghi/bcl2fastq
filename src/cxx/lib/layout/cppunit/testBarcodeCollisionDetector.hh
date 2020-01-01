/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file barcodeCollisionDetector.hh
 *
 * \brief BarcodeCollisionDetector cppunit test declarations.
 *
 * \author Aaron Day
 */


#ifndef BCL2FASTQ_LAYOUT_TEST_BARCODE_COLLISION_DETECTOR_HH
#define BCL2FASTQ_LAYOUT_TEST_BARCODE_COLLISION_DETECTOR_HH


#include <cppunit/extensions/HelperMacros.h>

#include <vector>

/// \brief Test suite for BarcodeCollisionDetector.
class TestBarcodeCollisionDetector : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE(TestBarcodeCollisionDetector);
    CPPUNIT_TEST(testSuccess);
    CPPUNIT_TEST(testMismatchLength);
    CPPUNIT_TEST(testNumComponents);
    CPPUNIT_TEST(testCollisionSameSample);
    CPPUNIT_TEST(testCollision);
    CPPUNIT_TEST_SUITE_END();

private:
    std::vector<std::size_t> componentMaxMismatches_;

public:
    TestBarcodeCollisionDetector() : componentMaxMismatches_() {}
    void setUp();
    void tearDown();
    void testSuccess();
    void testMismatchLength();
    void testNumComponents();
    void testCollisionSameSample();
    void testCollision();
};


#endif // BCL2FASTQ_LAYOUT_TEST_BARCODE_COLLISION_DETECTOR_HH
