/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file testAdapterLocator.cpp
 *
 * \brief AdapterLocator cppunit test declarations.
 *
 * \author Aaron Day
 */


#include <string>

#include "RegistryName.hh"
#include "testAdapterLocator.hh"


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(TestAdapterLocator, registryName("AdapterLocator"));

void TestAdapterLocator::reset()
{
    expectedResults_.clear();
    bclContainer_.bcls_.clear();
    bclContainer_.bcls_.resize(100);
}

void TestAdapterLocator::addToBclContainer(const std::string& read,
                                           size_t             expectedResult)
{
    BOOST_ASSERT(bclContainer_.bcls_.size() == read.size());

    expectedResults_.push_back(expectedResult);

    for (size_t i = 0; i < read.size(); ++i)
    {
        switch (read[i])
        {
            case 'A':
                bclContainer_.bcls_[i].bcls_.push_back('@');
                break;
            case 'C':
                bclContainer_.bcls_[i].bcls_.push_back('A');
                break;
            case 'G':
                bclContainer_.bcls_[i].bcls_.push_back('B');
                break;
            case 'T':
                bclContainer_.bcls_[i].bcls_.push_back('C');
                break;
            case 'N':
                bclContainer_.bcls_[i].bcls_.push_back('\0');
                break;
            default:
                BOOST_ASSERT(false);
                break;
        }
    }
}

void TestAdapterLocator::setUp()
{
}

void TestAdapterLocator::tearDown()
{
}

void TestAdapterLocator::testAll()
{
    testShortAdapter();

    reset();

    testLongAdapter();

    reset();

    testBeginningAdapter();

    reset();

    // This test should pass:
    // test2Mismatch();
}

void TestAdapterLocator::testShortAdapter()
{
    addToBclContainer(
        "NNNNNNNNNNNNNNNNNNNNNNNNNGGCCGGGAGTTGGGCGAGTACGGGCTGCAGGCATACACTGAAGTGAAAACTGCTGTCTCTTATTACACATCTCCG",
        77);

    addToBclContainer(
        "NNNNNNNNNNNNNNNNNNNNNNNNNTATCCTCTTTCTTCTCCTTTCCCTTTTCTTCCCCCTTNNNNNNNNCNNNNNNNNNNNNNNCTGTCTCTTATACAC",
        85);

    addToBclContainer(
        "NNNNNNNNNNNNNNNNNNNNNNNNNNACGTAGGTATTCTTTCTTTCTTCTTTTTCCTTTTCTTCCTTNCTCCTNNNNNNNNNNNNNNNNNNNNNNNNNNN",
        100);

    addToBclContainer(
        "NNNNNNNNNNNNNNNNNNNNNNNNNCGTCCACACCCTCCCCGGCCAGCCATGCACCGTTCTCTTATACACATCTCCGAGCCCACGAGACAGGGAGAAATT",
         54);

    addToBclContainer(
        "NNNNNNNNNNNNNNNNNNNNNNNNNAAGCACAGGGCAATCATATTGTGCTTTTCACTGTCTCCAGGACATGTCTCTTATACACATCTCCGAGCCCACGAG",
        67);

    addToBclContainer(
        "NNNNNNNNNNNNNNNNNNNNNNNGCCTCATTGCCACAGCTGTTCTCTTATACACATCTCCGAGCCCACGAGACTAAGGCGAATCTCGTATGCCGTCTTCT",
        38);

    addToBclContainer(
        "NNNNNNNNNNNNNNNNNNNNNNNAAAACCACTTTGTCCTCATCTCATTCCCCAAAACCGGCAGGCCCTCTCAGTTCTTTTTCCTGTCTTTACTTCTCTTA",
        91);

    addToBclContainer(
        "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNC",
        99);

    bcl2fastq::conversion::AdapterLocatorWithIndels<bcl2fastq::common::NumBasesPerByte::ONE> adapterLocator(
        "CTGTCTCTTATACACATCT",
        .8,
        30);

    verifyAdapterLocator(adapterLocator);
}

void TestAdapterLocator::testLongAdapter()
{
    addToBclContainer(
        "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNC",
        99);

    addToBclContainer(
        "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCTGTCTCTTATACACATCTACTGCTAGTCGCGTTACGGTTNNNNNNNNNNNNNNN",
        45);

    bcl2fastq::conversion::AdapterLocatorWithIndels<bcl2fastq::common::NumBasesPerByte::ONE> adapterLocator(
        "CTGTCTCTTATACACATCTACTGCTAGTCGCGTTACGGTT",
        .8,
        30);

    verifyAdapterLocator(adapterLocator);
}

void TestAdapterLocator::testBeginningAdapter()
{
    addToBclContainer(
        "ATCGGAAGACACACGTCTGAACTCCAGTCACATTCCTCGATCTCGTATGCCGTTTTCTGCTTGAAAAAACCAACAACACCTGGGAACAAGCGAATCACCN",
        0);

    bcl2fastq::conversion::AdapterLocatorWithIndels<bcl2fastq::common::NumBasesPerByte::ONE> adapterLocator(
        "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
        .9,
        35);

    verifyAdapterLocator(adapterLocator);
}

void TestAdapterLocator::test2Mismatch()
{
    addToBclContainer(
        "TGTAAACATTTCTTTAGTAGAATCTGCAAGTTGATATTTAGATAGCTAGGAAGATTTCCTTGGAAACGGGAATATCTTAATATAAACTCTAGACGGAAGC",
        100);

    bcl2fastq::conversion::AdapterLocatorWithIndels<bcl2fastq::common::NumBasesPerByte::ONE> adapterLocator(
        "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
        .9,
        35);

    verifyAdapterLocator(adapterLocator);
}

void TestAdapterLocator::verifyAdapterLocator(bcl2fastq::conversion::AdapterLocator<bcl2fastq::common::NumBasesPerByte::ONE>& adapterLocator)
{
    for (unsigned int i = 0; i < bclContainer_.bcls_.front().bcls_.size(); ++i)
    {
        bcl2fastq::conversion::FastqConstIterator<bcl2fastq::common::NumBasesPerByte::ONE> begin(bclContainer_, i);
        bcl2fastq::conversion::FastqConstIterator<bcl2fastq::common::NumBasesPerByte::ONE> end(bclContainer_, i);
        end += bclContainer_.bcls_.size();

        bcl2fastq::conversion::FastqConstIterator<bcl2fastq::common::NumBasesPerByte::ONE> pos = adapterLocator.identifyAdapter(begin, end);

        CPPUNIT_ASSERT_EQUAL(std::distance(begin, pos), expectedResults_[i]);
    }
}

