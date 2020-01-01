/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file testUseBasesMask.cpp
 *
 * \brief UseBasesMask cppunit test declarations.
 *
 * \author Aaron Day
 */


#include <string>

#include <boost/foreach.hpp>
#include <boost/assign/list_of.hpp>

#include "RegistryName.hh"
#include "testUseBasesMask.hh"

#include "layout/UseBasesMask.hh"

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(TestUseBasesMask, registryName("UseBasesMask"));


void TestUseBasesMask::setUp()
{
    using boost::assign::list_of;

    readMetadata_.push_back( list_of(std::make_pair(1,51))
                                    (std::make_pair(52,58))
                                    (std::make_pair(59,65))
                                    (std::make_pair(66,116)) );
    readMetadata_.back()[2].readType_ = bcl2fastq::common::ReadType::INDEX;
    readMetadata_.back()[3].readType_ = bcl2fastq::common::ReadType::INDEX;

    readMetadata_.push_back( list_of(std::make_pair(1,52))
                                    (std::make_pair(53,63)) );

    readMetadata_.push_back( list_of(std::make_pair(1,12))
                                    (std::make_pair(13,19)) );

    readMetadata_.push_back( list_of(std::make_pair(1,21))
                                    (std::make_pair(22,42)) );

    readMetadata_.push_back( list_of(std::make_pair(1,100)) );

    readMetadata_.push_back( list_of(std::make_pair(5,5)) );
}

void TestUseBasesMask::tearDown()
{
}

void TestUseBasesMask::testAllButOne()
{
    bcl2fastq::layout::UseBasesMask useBases("y50n,I6n,I6n,y50n",
                                             readMetadata_[0].begin(),
                                             readMetadata_[0].end());

    bcl2fastq::layout::ReadMetadataContainer::const_iterator iter = useBases.readMetadataBegin();

    CPPUNIT_ASSERT(iter != useBases.readMetadataEnd());
    CPPUNIT_ASSERT_EQUAL(iter->readNumber_, 1UL);
    CPPUNIT_ASSERT_EQUAL(iter->readType_, bcl2fastq::common::ReadType::DATA);
    CPPUNIT_ASSERT_EQUAL(iter->firstCycle_, 1UL);
    CPPUNIT_ASSERT_EQUAL(iter->lastCycle_, 50UL);

    ++iter;
    CPPUNIT_ASSERT(iter != useBases.readMetadataEnd());
    CPPUNIT_ASSERT_EQUAL(iter->readNumber_, 1UL);
    CPPUNIT_ASSERT_EQUAL(iter->readType_, bcl2fastq::common::ReadType::INDEX);
    CPPUNIT_ASSERT_EQUAL(iter->firstCycle_, 52UL);
    CPPUNIT_ASSERT_EQUAL(iter->lastCycle_, 57UL);

    ++iter;
    CPPUNIT_ASSERT(iter != useBases.readMetadataEnd());
    CPPUNIT_ASSERT_EQUAL(iter->readNumber_, 2UL);
    CPPUNIT_ASSERT_EQUAL(iter->readType_, bcl2fastq::common::ReadType::INDEX);
    CPPUNIT_ASSERT_EQUAL(iter->firstCycle_, 59UL);
    CPPUNIT_ASSERT_EQUAL(iter->lastCycle_, 64UL);

    ++iter;
    CPPUNIT_ASSERT(iter != useBases.readMetadataEnd());
    CPPUNIT_ASSERT_EQUAL(iter->readNumber_, 2UL);
    CPPUNIT_ASSERT_EQUAL(iter->readType_, bcl2fastq::common::ReadType::DATA);
    CPPUNIT_ASSERT_EQUAL(iter->firstCycle_, 66UL);
    CPPUNIT_ASSERT_EQUAL(iter->lastCycle_, 115UL);

    ++iter;
    CPPUNIT_ASSERT(iter == useBases.readMetadataEnd());
}

void TestUseBasesMask::testMaskingR1()
{
    bcl2fastq::layout::UseBasesMask useBases("y50nn,y10n",
                                             readMetadata_[1].begin(),
                                             readMetadata_[1].end());

    bcl2fastq::layout::ReadMetadataContainer::const_iterator iter = useBases.readMetadataBegin();

    CPPUNIT_ASSERT(iter != useBases.readMetadataEnd());
    CPPUNIT_ASSERT_EQUAL(iter->readNumber_, 1UL);
    CPPUNIT_ASSERT_EQUAL(iter->readType_, bcl2fastq::common::ReadType::DATA);
    CPPUNIT_ASSERT_EQUAL(iter->firstCycle_, 1UL);
    CPPUNIT_ASSERT_EQUAL(iter->lastCycle_, 50UL);

    ++iter;
    CPPUNIT_ASSERT(iter != useBases.readMetadataEnd());
    CPPUNIT_ASSERT_EQUAL(iter->readNumber_, 2UL);
    CPPUNIT_ASSERT_EQUAL(iter->readType_, bcl2fastq::common::ReadType::DATA);
    CPPUNIT_ASSERT_EQUAL(iter->firstCycle_, 53UL);
    CPPUNIT_ASSERT_EQUAL(iter->lastCycle_, 62UL);

    ++iter;
    CPPUNIT_ASSERT(iter == useBases.readMetadataEnd());
}

void TestUseBasesMask::testMaskingR2()
{
    bcl2fastq::layout::UseBasesMask useBases("y12,y5nn",
                                             readMetadata_[2].begin(),
                                             readMetadata_[2].end());

    bcl2fastq::layout::ReadMetadataContainer::const_iterator iter = useBases.readMetadataBegin();

    CPPUNIT_ASSERT(iter != useBases.readMetadataEnd());
    CPPUNIT_ASSERT_EQUAL(iter->readNumber_, 1UL);
    CPPUNIT_ASSERT_EQUAL(iter->readType_, bcl2fastq::common::ReadType::DATA);
    CPPUNIT_ASSERT_EQUAL(iter->firstCycle_, 1UL);
    CPPUNIT_ASSERT_EQUAL(iter->lastCycle_, 12UL);

    ++iter;
    CPPUNIT_ASSERT(iter != useBases.readMetadataEnd());
    CPPUNIT_ASSERT_EQUAL(iter->readNumber_, 2UL);
    CPPUNIT_ASSERT_EQUAL(iter->readType_, bcl2fastq::common::ReadType::DATA);
    CPPUNIT_ASSERT_EQUAL(iter->firstCycle_, 13UL);
    CPPUNIT_ASSERT_EQUAL(iter->lastCycle_, 17UL);

    ++iter;
    CPPUNIT_ASSERT(iter == useBases.readMetadataEnd());
}

void TestUseBasesMask::testMixedDataIndex()
{
    bcl2fastq::layout::UseBasesMask useBases("y10i10n,i10y11",
                                             readMetadata_[3].begin(),
                                             readMetadata_[3].end());

    bcl2fastq::layout::ReadMetadataContainer::const_iterator iter = useBases.readMetadataBegin();

    CPPUNIT_ASSERT(iter != useBases.readMetadataEnd());
    CPPUNIT_ASSERT_EQUAL(iter->readNumber_, 1UL);
    CPPUNIT_ASSERT_EQUAL(iter->readType_, bcl2fastq::common::ReadType::DATA);
    CPPUNIT_ASSERT_EQUAL(iter->firstCycle_, 1UL);
    CPPUNIT_ASSERT_EQUAL(iter->lastCycle_, 10UL);

    ++iter;
    CPPUNIT_ASSERT(iter != useBases.readMetadataEnd());
    CPPUNIT_ASSERT_EQUAL(iter->readNumber_, 1UL);
    CPPUNIT_ASSERT_EQUAL(iter->readType_, bcl2fastq::common::ReadType::INDEX);
    CPPUNIT_ASSERT_EQUAL(iter->firstCycle_, 11UL);
    CPPUNIT_ASSERT_EQUAL(iter->lastCycle_, 20UL);

    ++iter;
    CPPUNIT_ASSERT(iter != useBases.readMetadataEnd());
    CPPUNIT_ASSERT_EQUAL(iter->readNumber_, 2UL);
    CPPUNIT_ASSERT_EQUAL(iter->readType_, bcl2fastq::common::ReadType::INDEX);
    CPPUNIT_ASSERT_EQUAL(iter->firstCycle_, 22UL);
    CPPUNIT_ASSERT_EQUAL(iter->lastCycle_, 31UL);

    ++iter;
    CPPUNIT_ASSERT(iter != useBases.readMetadataEnd());
    CPPUNIT_ASSERT_EQUAL(iter->readNumber_, 2UL);
    CPPUNIT_ASSERT_EQUAL(iter->readType_, bcl2fastq::common::ReadType::DATA);
    CPPUNIT_ASSERT_EQUAL(iter->firstCycle_, 32UL);
    CPPUNIT_ASSERT_EQUAL(iter->lastCycle_, 42UL);

    ++iter;
    CPPUNIT_ASSERT(iter == useBases.readMetadataEnd());
}

void TestUseBasesMask::testExpansion()
{
    bcl2fastq::layout::UseBasesMask useBases("y10i*n,i10y*",
                                             readMetadata_[3].begin(),
                                             readMetadata_[3].end());

    bcl2fastq::layout::ReadMetadataContainer::const_iterator iter = useBases.readMetadataBegin();

    CPPUNIT_ASSERT(iter != useBases.readMetadataEnd());
    CPPUNIT_ASSERT_EQUAL(iter->readNumber_, 1UL);
    CPPUNIT_ASSERT_EQUAL(iter->readType_, bcl2fastq::common::ReadType::DATA);
    CPPUNIT_ASSERT_EQUAL(iter->firstCycle_, 1UL);
    CPPUNIT_ASSERT_EQUAL(iter->lastCycle_, 10UL);

    ++iter;
    CPPUNIT_ASSERT(iter != useBases.readMetadataEnd());
    CPPUNIT_ASSERT_EQUAL(iter->readNumber_, 1UL);
    CPPUNIT_ASSERT_EQUAL(iter->readType_, bcl2fastq::common::ReadType::INDEX);
    CPPUNIT_ASSERT_EQUAL(iter->firstCycle_, 11UL);
    CPPUNIT_ASSERT_EQUAL(iter->lastCycle_, 20UL);

    ++iter;
    CPPUNIT_ASSERT(iter != useBases.readMetadataEnd());
    CPPUNIT_ASSERT_EQUAL(iter->readNumber_, 2UL);
    CPPUNIT_ASSERT_EQUAL(iter->readType_, bcl2fastq::common::ReadType::INDEX);
    CPPUNIT_ASSERT_EQUAL(iter->firstCycle_, 22UL);
    CPPUNIT_ASSERT_EQUAL(iter->lastCycle_, 31UL);

    ++iter;
    CPPUNIT_ASSERT(iter != useBases.readMetadataEnd());
    CPPUNIT_ASSERT_EQUAL(iter->readNumber_, 2UL);
    CPPUNIT_ASSERT_EQUAL(iter->readType_, bcl2fastq::common::ReadType::DATA);
    CPPUNIT_ASSERT_EQUAL(iter->firstCycle_, 32UL);
    CPPUNIT_ASSERT_EQUAL(iter->lastCycle_, 42UL);

    ++iter;
    CPPUNIT_ASSERT(iter == useBases.readMetadataEnd());

    // Test * for 1 base
    bcl2fastq::layout::UseBasesMask useBases2("i20n*,y*",
                                             readMetadata_[3].begin(),
                                             readMetadata_[3].end());

    bcl2fastq::layout::ReadMetadataContainer::const_iterator iter2 = useBases2.readMetadataBegin();

    CPPUNIT_ASSERT(iter2 != useBases2.readMetadataEnd());
    CPPUNIT_ASSERT_EQUAL(iter2->readNumber_, 1UL);
    CPPUNIT_ASSERT_EQUAL(iter2->readType_, bcl2fastq::common::ReadType::INDEX);
    CPPUNIT_ASSERT_EQUAL(iter2->firstCycle_, 1UL);
    CPPUNIT_ASSERT_EQUAL(iter2->lastCycle_, 20UL);

    ++iter2;
    CPPUNIT_ASSERT(iter2 != useBases2.readMetadataEnd());
    CPPUNIT_ASSERT_EQUAL(iter2->readNumber_, 1UL);
    CPPUNIT_ASSERT_EQUAL(iter2->readType_, bcl2fastq::common::ReadType::DATA);
    CPPUNIT_ASSERT_EQUAL(iter2->firstCycle_, 22UL);
    CPPUNIT_ASSERT_EQUAL(iter2->lastCycle_, 42UL);

    ++iter2;
    CPPUNIT_ASSERT(iter2 == useBases2.readMetadataEnd());
}

void TestUseBasesMask::testMixedDataIndexSingleRead()
{
    bcl2fastq::layout::UseBasesMask useBases("y40i10y40i10",
                                             readMetadata_[4].begin(),
                                             readMetadata_[4].end());

    bcl2fastq::layout::ReadMetadataContainer::const_iterator iter = useBases.readMetadataBegin();

    CPPUNIT_ASSERT(iter != useBases.readMetadataEnd());
    CPPUNIT_ASSERT_EQUAL(iter->readNumber_, 1UL);
    CPPUNIT_ASSERT_EQUAL(iter->readType_, bcl2fastq::common::ReadType::DATA);
    CPPUNIT_ASSERT_EQUAL(iter->firstCycle_, 1UL);
    CPPUNIT_ASSERT_EQUAL(iter->lastCycle_, 40UL);

    ++iter;
    CPPUNIT_ASSERT(iter != useBases.readMetadataEnd());
    CPPUNIT_ASSERT_EQUAL(iter->readNumber_, 1UL);
    CPPUNIT_ASSERT_EQUAL(iter->readType_, bcl2fastq::common::ReadType::INDEX);
    CPPUNIT_ASSERT_EQUAL(iter->firstCycle_, 41UL);
    CPPUNIT_ASSERT_EQUAL(iter->lastCycle_, 50UL);

    ++iter;
    CPPUNIT_ASSERT(iter != useBases.readMetadataEnd());
    CPPUNIT_ASSERT_EQUAL(iter->readNumber_, 2UL);
    CPPUNIT_ASSERT_EQUAL(iter->readType_, bcl2fastq::common::ReadType::DATA);
    CPPUNIT_ASSERT_EQUAL(iter->firstCycle_, 51UL);
    CPPUNIT_ASSERT_EQUAL(iter->lastCycle_, 90UL);

    ++iter;
    CPPUNIT_ASSERT(iter != useBases.readMetadataEnd());
    CPPUNIT_ASSERT_EQUAL(iter->readNumber_, 2UL);
    CPPUNIT_ASSERT_EQUAL(iter->readType_, bcl2fastq::common::ReadType::INDEX);
    CPPUNIT_ASSERT_EQUAL(iter->firstCycle_, 91UL);
    CPPUNIT_ASSERT_EQUAL(iter->lastCycle_, 100UL);

    ++iter;
    CPPUNIT_ASSERT(iter == useBases.readMetadataEnd());
}

void TestUseBasesMask::testIsomorphic()
{
    CPPUNIT_ASSERT_THROW(
        bcl2fastq::layout::UseBasesMask useBases("y50n,y10", readMetadata_[1].begin(), readMetadata_[1].end()),
        bcl2fastq::layout::UseBasesMaskFormatError
    );
}

void TestUseBasesMask::testMismatchingCyclesInRead()
{
    // The number of cycles in each read do not match, but the total cycles match so it is ok.
    bcl2fastq::layout::UseBasesMask useBases("y50n,y12", readMetadata_[1].begin(), readMetadata_[1].end());

    bcl2fastq::layout::ReadMetadataContainer::const_iterator iter = useBases.readMetadataBegin();

    CPPUNIT_ASSERT(iter != useBases.readMetadataEnd());
    CPPUNIT_ASSERT_EQUAL(iter->readNumber_, 1UL);
    CPPUNIT_ASSERT_EQUAL(iter->readType_, bcl2fastq::common::ReadType::DATA);
    CPPUNIT_ASSERT_EQUAL(iter->firstCycle_, 1UL);
    CPPUNIT_ASSERT_EQUAL(iter->lastCycle_, 50UL);

    ++iter;

    CPPUNIT_ASSERT(iter != useBases.readMetadataEnd());
    CPPUNIT_ASSERT_EQUAL(iter->readNumber_, 2UL);
    CPPUNIT_ASSERT_EQUAL(iter->readType_, bcl2fastq::common::ReadType::DATA);
    CPPUNIT_ASSERT_EQUAL(iter->firstCycle_, 52UL);
    CPPUNIT_ASSERT_EQUAL(iter->lastCycle_, 63UL);

    ++iter;
    CPPUNIT_ASSERT(iter == useBases.readMetadataEnd());
}

void TestUseBasesMask::testZero()
{
    CPPUNIT_ASSERT_THROW(
        bcl2fastq::layout::UseBasesMask useBases("y0n,y10", readMetadata_[5].begin(), readMetadata_[5].end()),
        bcl2fastq::layout::UseBasesMaskFormatError
    );
}

void TestUseBasesMask::testMissingRead()
{
    bcl2fastq::layout::UseBasesMask useBases("Y*,I5nn,n*,Y*",
                                             readMetadata_[0].begin(),
                                             readMetadata_[0].end());

    bcl2fastq::layout::ReadMetadataContainer::const_iterator iter = useBases.readMetadataBegin();

    CPPUNIT_ASSERT(iter != useBases.readMetadataEnd());
    CPPUNIT_ASSERT_EQUAL(iter->readNumber_, 1UL);
    CPPUNIT_ASSERT_EQUAL(iter->readType_, bcl2fastq::common::ReadType::DATA);
    CPPUNIT_ASSERT_EQUAL(iter->firstCycle_, 1UL);
    CPPUNIT_ASSERT_EQUAL(iter->lastCycle_, 51UL);

    ++iter;
    CPPUNIT_ASSERT(iter != useBases.readMetadataEnd());
    CPPUNIT_ASSERT_EQUAL(iter->readNumber_, 1UL);
    CPPUNIT_ASSERT_EQUAL(iter->readType_, bcl2fastq::common::ReadType::INDEX);
    CPPUNIT_ASSERT_EQUAL(iter->firstCycle_, 52UL);
    CPPUNIT_ASSERT_EQUAL(iter->lastCycle_, 56UL);

    ++iter;
    CPPUNIT_ASSERT(iter != useBases.readMetadataEnd());
    CPPUNIT_ASSERT_EQUAL(iter->readNumber_, 2UL);
    CPPUNIT_ASSERT_EQUAL(iter->readType_, bcl2fastq::common::ReadType::DATA);
    CPPUNIT_ASSERT_EQUAL(iter->firstCycle_, 66UL);
    CPPUNIT_ASSERT_EQUAL(iter->lastCycle_, 116UL);

    ++iter;
    CPPUNIT_ASSERT(iter == useBases.readMetadataEnd());
}

void TestUseBasesMask::testMissingAsteriskNone()
{
    // The "n*" matches 0 bases
    bcl2fastq::layout::UseBasesMask useBases("Y*,I5nn,i7n*,Y*",
                                             readMetadata_[0].begin(),
                                             readMetadata_[0].end());

    bcl2fastq::layout::ReadMetadataContainer::const_iterator iter = useBases.readMetadataBegin();

    CPPUNIT_ASSERT(iter != useBases.readMetadataEnd());
    CPPUNIT_ASSERT_EQUAL(iter->readNumber_, 1UL);
    CPPUNIT_ASSERT_EQUAL(iter->readType_, bcl2fastq::common::ReadType::DATA);
    CPPUNIT_ASSERT_EQUAL(iter->firstCycle_, 1UL);
    CPPUNIT_ASSERT_EQUAL(iter->lastCycle_, 51UL);

    ++iter;
    CPPUNIT_ASSERT(iter != useBases.readMetadataEnd());
    CPPUNIT_ASSERT_EQUAL(iter->readNumber_, 1UL);
    CPPUNIT_ASSERT_EQUAL(iter->readType_, bcl2fastq::common::ReadType::INDEX);
    CPPUNIT_ASSERT_EQUAL(iter->firstCycle_, 52UL);
    CPPUNIT_ASSERT_EQUAL(iter->lastCycle_, 56UL);

    ++iter;
    CPPUNIT_ASSERT(iter != useBases.readMetadataEnd());
    CPPUNIT_ASSERT_EQUAL(iter->readNumber_, 2UL);
    CPPUNIT_ASSERT_EQUAL(iter->readType_, bcl2fastq::common::ReadType::INDEX);
    CPPUNIT_ASSERT_EQUAL(iter->firstCycle_, 59UL);
    CPPUNIT_ASSERT_EQUAL(iter->lastCycle_, 65UL);

    ++iter;
    CPPUNIT_ASSERT(iter != useBases.readMetadataEnd());
    CPPUNIT_ASSERT_EQUAL(iter->readNumber_, 2UL);
    CPPUNIT_ASSERT_EQUAL(iter->readType_, bcl2fastq::common::ReadType::DATA);
    CPPUNIT_ASSERT_EQUAL(iter->firstCycle_, 66UL);
    CPPUNIT_ASSERT_EQUAL(iter->lastCycle_, 116UL);

    ++iter;
    CPPUNIT_ASSERT(iter == useBases.readMetadataEnd());
}
