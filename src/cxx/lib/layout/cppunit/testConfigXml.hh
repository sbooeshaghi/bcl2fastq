/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file testConfigXml.hh
 *
 * \brief Config.xml helper cppunit test declarations.
 *
 * \author Marek Balint
 */


#ifndef BCL2FASTQ_LAYOUT_TEST_CONFIGXML_HH
#define BCL2FASTQ_LAYOUT_TEST_CONFIGXML_HH


#include <string>

#include <cppunit/extensions/HelperMacros.h>

#include "layout/ConfigXml.hh"


/// \brief Test suite for ConfigXml.
class TestConfigXml : public CppUnit::TestFixture
{
private:

    CPPUNIT_TEST_SUITE(TestConfigXml);
    CPPUNIT_TEST(testHcs);
    CPPUNIT_TEST_SUITE_END();

public:

    void setUp();
    void tearDown();
    void testHcs();
};


#endif // #ifndef BCL2FASTQ_LAYOUT_TEST_CONFIGXML_HH


