/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file testExceptions.cpp
 *
 * \brief Commmon exception mechanism cppunit test declarations.
 *
 * \author Roman Petrovski
 * \author Marek Balint
 */


#include <iostream>
#include <sstream>
#include <string>

#include "RegistryName.hh"
#include "testExceptions.hh"


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(TestExceptions, registryName("Exceptions"));


void TestExceptions::setUp()
{
}

void TestExceptions::tearDown()
{
}

void TestExceptions::testErrorNumber()
{
    CPPUNIT_ASSERT(true);
}
