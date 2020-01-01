/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file testRunInfoXml.hh
 *
 * \brief RunInfo.xml helper cppunit test declarations.
 *
 * \author Marek Balint
 */


#ifndef BCL2FASTQ_LAYOUT_TEST_RUNINFOXML_HH
#define BCL2FASTQ_LAYOUT_TEST_RUNINFOXML_HH


#include <string>

#include <cppunit/extensions/HelperMacros.h>

#include "layout/RunInfoXml.hh"


/// \brief Test suite for RunInfoXml.
class TestRunInfoXml : public CppUnit::TestFixture
{
private:

    CPPUNIT_TEST_SUITE(TestRunInfoXml);
    CPPUNIT_TEST(testNova);
    CPPUNIT_TEST_SUITE_END();

public:

    void setUp();
    void tearDown();
    void testNova();
};


#endif // #ifndef BCL2FASTQ_LAYOUT_TEST_RUNINFOXML_HH


