/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file testExceptions.hh
 *
 * \brief Commmon exception mechanism cppunit test declarations.
 *
 * \author Roman Petrovski
 * \author Marek Balint
 */

                                                                    
#ifndef BCL2FASTQ_COMMON_TEST_EXCEPTIONS_HH
#define BCL2FASTQ_COMMON_TEST_EXCEPTIONS_HH


#include <string>

#include <cppunit/extensions/HelperMacros.h>

#include "common/Exceptions.hh"


/// \brief Test suite for exceptions.
class TestExceptions : public CppUnit::TestFixture
{
private:

    CPPUNIT_TEST_SUITE(TestExceptions);
    CPPUNIT_TEST(testErrorNumber);
    CPPUNIT_TEST_SUITE_END();

public:

    void setUp();
    void tearDown();
    void testErrorNumber();
};


#endif // #ifndef BCL2FASTQ_COMMON_TEST_EXCEPTIONS_HH


