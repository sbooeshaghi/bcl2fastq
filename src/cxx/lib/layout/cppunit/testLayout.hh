/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file testLayout.hh
 *
 * \brief Layout cppunit test declarations.
 *
 * \author Aaron Day
 */


#ifndef BCL2FASTQ_LAYOUT_TEST_LAYOUT_HH
#define BCL2FASTQ_LAYOUT_TEST_LAYOUT_HH


#include <string>
#include <memory>

#include "config/SampleSheetCsv.hh"
#include "layout/RunInfoXml.hh"

#include <cppunit/extensions/HelperMacros.h>


/// \brief Test suite for Layout.
class TestLayout : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE(TestLayout);
    CPPUNIT_TEST(testDetectReadLayoutDefaultCycles);
    CPPUNIT_TEST(testUmiCyclesAtBeginning);
    CPPUNIT_TEST(testUmiCyclesAtEnd);
    CPPUNIT_TEST(testUmiCyclesInMiddle);
    CPPUNIT_TEST(testUmiCyclesAtBeginningOfMaskedRead);
    CPPUNIT_TEST(testUmiCyclesAtEndOfMaskedRead);
    
    CPPUNIT_TEST(testDetectReadLayoutDefaultCyclesTrim);
    CPPUNIT_TEST(testUmiCyclesAtBeginningTrim);
    CPPUNIT_TEST(testUmiCyclesAtEndTrim);
    CPPUNIT_TEST(testUmiCyclesInMiddleTrim);
    CPPUNIT_TEST(testUmiCyclesAtBeginningOfMaskedReadTrim);
    CPPUNIT_TEST(testUmiCyclesAtEndOfMaskedReadTrim);
    CPPUNIT_TEST_SUITE_END();

private:
    struct ExpectedCycleInfosForRead
    {
        ExpectedCycleInfosForRead(size_t numCycles,
                                  bcl2fastq::common::ReadNumber readNumber,
                                  bcl2fastq::common::CycleNumber firstCycle)
        : numCycles_(numCycles)
        , readNumber_(readNumber)
        , firstCycle_(firstCycle)
        {
        }

        size_t numCycles_;
        bcl2fastq::common::ReadNumber readNumber_;
        bcl2fastq::common::CycleNumber firstCycle_;
    };

    void getLaneInfos(bcl2fastq::layout::LaneInfosContainer& laneInfos);
    void assertExpectedCycles(const bcl2fastq::layout::LaneInfo& laneInfo,
                              const std::vector<ExpectedCycleInfosForRead>& expectedCycles,
                              const std::vector<ExpectedCycleInfosForRead>& expectedUnmaskedCycles,
                              const std::vector<ExpectedCycleInfosForRead>& expectedCyclesToLoad,
                              const std::vector<bcl2fastq::common::CycleRange>& umiCycles,
                              const std::vector<bcl2fastq::common::CycleNumber>& bclBufferOffsets);

    template<typename T>
    void verifyCycles(bcl2fastq::common::CycleNumber startCycle,
                      const T& cycleInfos);

    std::vector<std::vector<size_t>> barcodeLengthsForLane_;
    std::shared_ptr<bcl2fastq::layout::RunInfoXml> runInfoXml_;
    std::shared_ptr<bcl2fastq::config::SampleSheetCsv> sampleSheet_;

public:
    TestLayout() : barcodeLengthsForLane_(), runInfoXml_(), sampleSheet_() { }
    void setUp();
    void tearDown();

    void testDetectReadLayoutDefaultCycles();
    void testUmiCyclesAtBeginning();
    void testUmiCyclesAtEnd();
    void testUmiCyclesInMiddle();
    void testUmiCyclesAtBeginningOfMaskedRead();
    void testUmiCyclesAtEndOfMaskedRead();

    void testDetectReadLayoutDefaultCyclesTrim();
    void testUmiCyclesAtBeginningTrim();
    void testUmiCyclesAtEndTrim();
    void testUmiCyclesInMiddleTrim();
    void testUmiCyclesAtBeginningOfMaskedReadTrim();
    void testUmiCyclesAtEndOfMaskedReadTrim();
};

#endif // BCL2FASTQ_LAYOUT_TEST_LAYOUT_HH


