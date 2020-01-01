/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file testAdapterLocator.hh
 *
 * \brief AdapterLocator cppunit test declarations.
 *
 * \author Aaron Day
 */


#ifndef BCL2FASTQ_CONVERSION_TEST_ADAPTER_LOCATOR_HH
#define BCL2FASTQ_CONVERSION_TEST_ADAPTER_LOCATOR_HH

#include "data/BclBuffer.hh"
#include "conversion/AdapterLocator.hh"

#include <string>
#include <vector>

#include <cppunit/extensions/HelperMacros.h>


/// \brief Test suite for AdapterLocaator.
class TestAdapterLocator : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE(TestAdapterLocator);
    CPPUNIT_TEST(testAll);
    CPPUNIT_TEST_SUITE_END();

public:
    TestAdapterLocator() : bclContainer_(), expectedResults_() { bclContainer_.bcls_.resize(100); }

    void setUp();
    void tearDown();
    void testAll();

private:

    void reset();

    void addToBclContainer(const std::string& read,
                           size_t             expectedResult);

    void testLongAdapter();
    void testShortAdapter();
    void testBeginningAdapter();
    void test2Mismatch();

    void verifyAdapterLocator(bcl2fastq::conversion::AdapterLocator<bcl2fastq::common::NumBasesPerByte::ONE>& adapterLocator);

    bcl2fastq::data::BclBuffer bclContainer_;
    std::vector<long>          expectedResults_;
};


#endif // BCL2FASTQ_CONVERSION_TEST_ADAPTER_LOCATOR_HH


