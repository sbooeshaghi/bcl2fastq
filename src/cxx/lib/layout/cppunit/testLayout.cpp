/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file testLayout.cpp
 *
 * \brief Layout cppunit test declarations.
 *
 * \author Aaron Day
 */


#include <string>

#include "RegistryName.hh"
#include "testLayout.hh"

#include "io/Xml.hh"
#include "layout/Layout.hh"

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(TestLayout, registryName("Layout"));

void TestLayout::setUp()
{
    // Turn off logging
    bcl2fastq::common::BCL2FASTQ_MIN_LOG_LEVEL = bcl2fastq::common::LogLevel::NONE;

    barcodeLengthsForLane_.push_back(std::vector<size_t>(1, 8));

    const std::string runInfoXmlData(
        "<?xml version=\"1.0\"?>\n"
        "<RunInfo xmlns:xsd=\"http://www.w3.org/2001/XMLSchema\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" Version=\"2\">\n"
        "  <Run Id=\"TestRun\" Number=\"42\">\n"
        "    <Flowcell>TestFlowcell</Flowcell>\n"
        "    <Instrument>CSSIM</Instrument>\n"
        "    <Date>8/20/2013 3:29:17 PM</Date>\n"
        "    <Reads>\n"
        "      <Read Number=\"1\" NumCycles=\"151\" IsIndexedRead=\"N\" />\n"
        "      <Read Number=\"2\" NumCycles=\"8\" IsIndexedRead=\"Y\" />\n"
        "      <Read Number=\"3\" NumCycles=\"151\" IsIndexedRead=\"N\" />\n"
        "    </Reads>\n"
        "    <FlowcellLayout LaneCount=\"4\" SurfaceCount=\"2\" SwathCount=\"3\" TileCount=\"11\" />\n"
        "  </Run>\n"
        "</RunInfo>\n"
    );

    runInfoXml_ = std::shared_ptr<bcl2fastq::layout::RunInfoXml>(
        new bcl2fastq::layout::RunInfoXml(bcl2fastq::io::parseXmlData(
            runInfoXmlData.begin(),
            runInfoXmlData.end())));

    const std::string sampleSheetData(
        "[Settings]\n"
        "trimUmi,1\n"
        "read1UmiLength,2\n"
        "read2UmiLength,151\n");

    sampleSheet_ = std::shared_ptr<bcl2fastq::config::SampleSheetCsv>(
        new bcl2fastq::config::SampleSheetCsv(
            bcl2fastq::common::parseCsvData(
                sampleSheetData.begin(),
                sampleSheetData.end()),
            boost::filesystem::path()));
}

void TestLayout::tearDown()
{
}

void TestLayout::getLaneInfos(bcl2fastq::layout::LaneInfosContainer& laneInfos)
{
    laneInfos.push_back(bcl2fastq::layout::LaneInfo(1));
    laneInfos.back().addSample(bcl2fastq::layout::SampleInfo());

    std::vector<std::string> barcodeComponents;
    barcodeComponents.push_back("AAAAGGGG");
    bcl2fastq::layout::Barcode barcode(barcodeComponents.begin(), barcodeComponents.end());
    bcl2fastq::layout::SampleInfo sampleInfo(1, "MySampleId", "MySampleName", "MyProject");
    sampleInfo.addBarcode(barcode);
    laneInfos.back().addSample(sampleInfo);
}

// expectedCycles: Cycles to put in the FASTQ sequence
// expectedUnmaskedCycles: All cycles in RunInfo.xml
// expectedCyclesToLoad: Cycles we need to load from BCL (including UMIs, which might not be present in the FASTQ sequence)
void TestLayout::assertExpectedCycles(const bcl2fastq::layout::LaneInfo& laneInfo,
                                      const std::vector<ExpectedCycleInfosForRead>& expectedCycles,
                                      const std::vector<ExpectedCycleInfosForRead>& expectedUnmaskedCycles,
                                      const std::vector<ExpectedCycleInfosForRead>& expectedCyclesToLoad,
                                      const std::vector<bcl2fastq::common::CycleRange>& umiCycles,
                                      const std::vector<bcl2fastq::common::CycleNumber>& bclBufferOffsets)
{
    for (const auto& expectedCyclesForRead : expectedCycles)
    {
        CPPUNIT_ASSERT_EQUAL(laneInfo.getReadInfos().at(expectedCyclesForRead.readNumber_).cycleInfos().size(), expectedCyclesForRead.numCycles_);

        verifyCycles(expectedCyclesForRead.firstCycle_,
                     laneInfo.getReadInfos().at(expectedCyclesForRead.readNumber_).cycleInfos());
    }

    for (const auto& expectedCyclesForRead : expectedUnmaskedCycles)
    {
        CPPUNIT_ASSERT_EQUAL(laneInfo.getReadInfos().at(expectedCyclesForRead.readNumber_).unmaskedCycleInfos().size(), expectedCyclesForRead.numCycles_);

        verifyCycles(expectedCyclesForRead.firstCycle_,
                     laneInfo.getReadInfos().at(expectedCyclesForRead.readNumber_).unmaskedCycleInfos());
    }

    for (const auto& expectedCyclesForRead : expectedCyclesToLoad)
    {
        CPPUNIT_ASSERT_EQUAL(laneInfo.getReadInfos().at(expectedCyclesForRead.readNumber_).cyclesToLoad().size(), expectedCyclesForRead.numCycles_);

        verifyCycles(expectedCyclesForRead.firstCycle_,
                     laneInfo.getReadInfos().at(expectedCyclesForRead.readNumber_).cyclesToLoad());
    }

    for (size_t i = 0; i < laneInfo.getReadInfos().size(); ++i)
    {
        CPPUNIT_ASSERT_EQUAL(laneInfo.getReadInfos().at(i).getUmiCycles().first, umiCycles[i].first);
        CPPUNIT_ASSERT_EQUAL(laneInfo.getReadInfos().at(i).getUmiCycles().second, umiCycles[i].second);
        CPPUNIT_ASSERT_EQUAL(laneInfo.getReadInfos().at(i).getBclBufferOffset(), bclBufferOffsets.at(i));
    }
}

template<typename T>
void TestLayout::verifyCycles(bcl2fastq::common::CycleNumber startCycle,
                              const T& cycleInfos)
{
    for (const auto& cycle : cycleInfos)
    {
        CPPUNIT_ASSERT_EQUAL(cycle.getNumber(), startCycle++);
    }
}

void TestLayout::testDetectReadLayoutDefaultCycles()
{
    bcl2fastq::config::BarcodeMismatchCountsContainer barcodeMismatches;

    bcl2fastq::layout::LaneInfosContainer laneInfos;
    getLaneInfos(laneInfos);

    bcl2fastq::layout::detectReadLayout(
        *runInfoXml_,
        *sampleSheet_, // used for adapters only
        std::vector<std::string>(), // use bases masks
        35,
        0, // read1StartCycle
        0, // read2StartCycle
        0, // read1EndCycle
        0, // read2EndCycle
        0, // read1UmiLength
        0, // read2UmiLength
        1, // read1UmiStartFromCycle
        1, // read2UmiStartFromCycle
        false, // trimUmi
        false,
        barcodeMismatches,
        laneInfos,
        barcodeLengthsForLane_);

    CPPUNIT_ASSERT_EQUAL(laneInfos.back().getReadInfos().size(),3UL);
    assertExpectedCycles(laneInfos.back(),
                         std::vector<ExpectedCycleInfosForRead>{ExpectedCycleInfosForRead(151, 0, 1),
                                                                ExpectedCycleInfosForRead(8, 1, 152),
                                                                ExpectedCycleInfosForRead(151, 2, 160)},
                         std::vector<ExpectedCycleInfosForRead>{ExpectedCycleInfosForRead(151, 0, 1),
                                                                ExpectedCycleInfosForRead(8, 1, 152),
                                                                ExpectedCycleInfosForRead(151, 2, 160)},
                         std::vector<ExpectedCycleInfosForRead>{ExpectedCycleInfosForRead(151, 0, 1),
                                                                ExpectedCycleInfosForRead(8, 1, 152),
                                                                ExpectedCycleInfosForRead(151, 2, 160)},
                         std::vector<bcl2fastq::common::CycleRange>{bcl2fastq::common::CycleRange(0, 0),
                                                                    bcl2fastq::common::CycleRange(0, 0),
                                                                    bcl2fastq::common::CycleRange(0, 0)},
                         std::vector<bcl2fastq::common::CycleNumber>{0, 0, 0});
}

void TestLayout::testUmiCyclesAtBeginning()
{
    bcl2fastq::config::BarcodeMismatchCountsContainer barcodeMismatches;

    bcl2fastq::layout::LaneInfosContainer laneInfos;
    getLaneInfos(laneInfos);

    bcl2fastq::layout::detectReadLayout(
        *runInfoXml_,
        *sampleSheet_, // used for adapters only
        std::vector<std::string>(), // use bases masks
        35,
        0, // read1StartCycle
        0, // read2StartCycle
        0, // read1EndCycle
        0, // read2EndCycle
        10, // read1UmiLength
        12, // read2UmiLength
        1, // read1UmiStartFromCycle
        1, // read2UmiStartFromCycle
        false, // trimUmi
        false,
        barcodeMismatches,
        laneInfos,
        barcodeLengthsForLane_);

    CPPUNIT_ASSERT_EQUAL(laneInfos.back().getReadInfos().size(),3UL);
    assertExpectedCycles(laneInfos.back(),
                         std::vector<ExpectedCycleInfosForRead>{ExpectedCycleInfosForRead(151, 0, 1),
                                                                ExpectedCycleInfosForRead(8, 1, 152),
                                                                ExpectedCycleInfosForRead(151, 2, 160)},
                         std::vector<ExpectedCycleInfosForRead>{ExpectedCycleInfosForRead(151, 0, 1),
                                                                ExpectedCycleInfosForRead(8, 1, 152),
                                                                ExpectedCycleInfosForRead(151, 2, 160)},
                         std::vector<ExpectedCycleInfosForRead>{ExpectedCycleInfosForRead(151, 0, 1),
                                                                ExpectedCycleInfosForRead(8, 1, 152),
                                                                ExpectedCycleInfosForRead(151, 2, 160)},
                         std::vector<bcl2fastq::common::CycleRange>{bcl2fastq::common::CycleRange(0, 10),
                                                                    bcl2fastq::common::CycleRange(0, 0),
                                                                    bcl2fastq::common::CycleRange(0, 12)},
                         std::vector<bcl2fastq::common::CycleNumber>{0, 0, 0});
}

void TestLayout::testUmiCyclesAtEnd()
{
    bcl2fastq::config::BarcodeMismatchCountsContainer barcodeMismatches;

    bcl2fastq::layout::LaneInfosContainer laneInfos;
    getLaneInfos(laneInfos);

    bcl2fastq::layout::detectReadLayout(
        *runInfoXml_,
        *sampleSheet_, // used for adapters only
        std::vector<std::string>(), // use bases masks
        35,
        0, // read1StartCycle
        0, // read2StartCycle
        0, // read1EndCycle
        0, // read2EndCycle
        10, // read1UmiLength
        12, // read2UmiLength
        142, // read1UmiStartFromCycle
        140, // read2UmiStartFromCycle
        false, // trimUmi
        false,
        barcodeMismatches,
        laneInfos,
        barcodeLengthsForLane_);

    CPPUNIT_ASSERT_EQUAL(laneInfos.back().getReadInfos().size(),3UL);
    assertExpectedCycles(laneInfos.back(),
                         std::vector<ExpectedCycleInfosForRead>{ExpectedCycleInfosForRead(151, 0, 1),
                                                                ExpectedCycleInfosForRead(8, 1, 152),
                                                                ExpectedCycleInfosForRead(151, 2, 160)},
                         std::vector<ExpectedCycleInfosForRead>{ExpectedCycleInfosForRead(151, 0, 1),
                                                                ExpectedCycleInfosForRead(8, 1, 152),
                                                                ExpectedCycleInfosForRead(151, 2, 160)},
                         std::vector<ExpectedCycleInfosForRead>{ExpectedCycleInfosForRead(151, 0, 1),
                                                                ExpectedCycleInfosForRead(8, 1, 152),
                                                                ExpectedCycleInfosForRead(151, 2, 160)},
                         std::vector<bcl2fastq::common::CycleRange>{bcl2fastq::common::CycleRange(141, 151),
                                                                    bcl2fastq::common::CycleRange(0, 0),
                                                                    bcl2fastq::common::CycleRange(139, 151)},
                         std::vector<bcl2fastq::common::CycleNumber>{0, 0, 0});
}

void TestLayout::testUmiCyclesInMiddle()
{
    bcl2fastq::config::BarcodeMismatchCountsContainer barcodeMismatches;

    bcl2fastq::layout::LaneInfosContainer laneInfos;
    getLaneInfos(laneInfos);

    CPPUNIT_ASSERT_THROW(
        bcl2fastq::layout::detectReadLayout(
            *runInfoXml_,
            *sampleSheet_, // used for adapters only
            std::vector<std::string>(), // use bases masks
            35,
            0, // read1StartCycle
            0, // read2StartCycle
            0, // read1EndCycle
            0, // read2EndCycle
            8, // read1UmiLength
            12, // read2UmiLength
            142, // read1UmiStartFromCycle
            140, // read2UmiStartFromCycle
            false, // trimUmi
            false,
            barcodeMismatches,
            laneInfos,
            barcodeLengthsForLane_),
        bcl2fastq::common::InputDataError);
}

void TestLayout::testUmiCyclesAtBeginningOfMaskedRead()
{
    bcl2fastq::config::BarcodeMismatchCountsContainer barcodeMismatches;

    bcl2fastq::layout::LaneInfosContainer laneInfos;
    getLaneInfos(laneInfos);

    bcl2fastq::layout::detectReadLayout(
        *runInfoXml_,
        *sampleSheet_, // used for adapters only
        std::vector<std::string>(), // use bases masks
        35,
        50, // read1StartCycle
        4, // read2StartCycle
        0, // read1EndCycle
        0, // read2EndCycle
        2, // read1UmiLength
        12, // read2UmiLength
        2, // read1UmiStartFromCycle
        3, // read2UmiStartFromCycle
        false, // trimUmi
        false,
        barcodeMismatches,
        laneInfos,
        barcodeLengthsForLane_);

    CPPUNIT_ASSERT_EQUAL(laneInfos.back().getReadInfos().size(),3UL);
    assertExpectedCycles(laneInfos.back(),
                         std::vector<ExpectedCycleInfosForRead>{ExpectedCycleInfosForRead(102, 0, 50),
                                                                ExpectedCycleInfosForRead(8, 1, 152),
                                                                ExpectedCycleInfosForRead(148, 2, 163)},
                         std::vector<ExpectedCycleInfosForRead>{ExpectedCycleInfosForRead(151, 0, 1),
                                                                ExpectedCycleInfosForRead(8, 1, 152),
                                                                ExpectedCycleInfosForRead(151, 2, 160)},
                         std::vector<ExpectedCycleInfosForRead>{},//ExpectedCycleInfosForRead(104, 0, 2),
                                                                //ExpectedCycleInfosForRead(8, 1, 152),
                                                                //ExpectedCycleInfosForRead(149, 2, 162)},
                         std::vector<bcl2fastq::common::CycleRange>{bcl2fastq::common::CycleRange(0, 2),
                                                                    bcl2fastq::common::CycleRange(0, 0),
                                                                    bcl2fastq::common::CycleRange(0, 12)},
                         std::vector<bcl2fastq::common::CycleNumber>{2, 0, 1});

    CPPUNIT_ASSERT_EQUAL(laneInfos.back().getReadInfos().at(0).cyclesToLoad().size(), 104UL);
}

void TestLayout::testUmiCyclesAtEndOfMaskedRead()
{
    bcl2fastq::config::BarcodeMismatchCountsContainer barcodeMismatches;

    bcl2fastq::layout::LaneInfosContainer laneInfos;
    getLaneInfos(laneInfos);

    bcl2fastq::layout::detectReadLayout(
        *runInfoXml_,
        *sampleSheet_, // used for adapters only
        std::vector<std::string>(), // use bases masks
        35,
        2, // read1StartCycle
        3, // read2StartCycle
        100, // read1EndCycle
        99, // read2EndCycle
        10, // read1UmiLength
        12, // read2UmiLength
        142, // read1UmiStartFromCycle
        140, // read2UmiStartFromCycle
        false, // trimUmi
        false,
        barcodeMismatches,
        laneInfos,
        barcodeLengthsForLane_);

    CPPUNIT_ASSERT_EQUAL(laneInfos.back().getReadInfos().size(),3UL);
    assertExpectedCycles(laneInfos.back(),
                         std::vector<ExpectedCycleInfosForRead>{ExpectedCycleInfosForRead(99, 0, 2),
                                                                ExpectedCycleInfosForRead(8, 1, 152),
                                                                ExpectedCycleInfosForRead(97, 2, 162)},
                         std::vector<ExpectedCycleInfosForRead>{ExpectedCycleInfosForRead(151, 0, 1),
                                                                ExpectedCycleInfosForRead(8, 1, 152),
                                                                ExpectedCycleInfosForRead(151, 2, 160)},
                         std::vector<ExpectedCycleInfosForRead>{},//ExpectedCycleInfosForRead(109, 0, 2),
                                                                //ExpectedCycleInfosForRead(8, 1, 152),
                                                                //ExpectedCycleInfosForRead(149, 2, 162)},
                         std::vector<bcl2fastq::common::CycleRange>{bcl2fastq::common::CycleRange(99, 109),
                                                                    bcl2fastq::common::CycleRange(0, 0),
                                                                    bcl2fastq::common::CycleRange(97, 109)},
                         std::vector<bcl2fastq::common::CycleNumber>{0, 0, 0});
}



//
// Repeate above but this time trim the UMI
//


void TestLayout::testDetectReadLayoutDefaultCyclesTrim()
{
    bcl2fastq::config::BarcodeMismatchCountsContainer barcodeMismatches;

    bcl2fastq::layout::LaneInfosContainer laneInfos;
    getLaneInfos(laneInfos);

    bcl2fastq::layout::detectReadLayout(
        *runInfoXml_,
        *sampleSheet_, // used for adapters only
        std::vector<std::string>(), // use bases masks
        35,
        0, // read1StartCycle
        0, // read2StartCycle
        0, // read1EndCycle
        0, // read2EndCycle
        0, // read1UmiLength
        0, // read2UmiLength
        1, // read1UmiStartFromCycle
        1, // read2UmiStartFromCycle
        true, // trimUmi
        false,
        barcodeMismatches,
        laneInfos,
        barcodeLengthsForLane_);

    CPPUNIT_ASSERT_EQUAL(laneInfos.back().getReadInfos().size(),3UL);
    assertExpectedCycles(laneInfos.back(),
                         std::vector<ExpectedCycleInfosForRead>{ExpectedCycleInfosForRead(151, 0, 1),
                                                                ExpectedCycleInfosForRead(8, 1, 152),
                                                                ExpectedCycleInfosForRead(151, 2, 160)},
                         std::vector<ExpectedCycleInfosForRead>{ExpectedCycleInfosForRead(151, 0, 1),
                                                                ExpectedCycleInfosForRead(8, 1, 152),
                                                                ExpectedCycleInfosForRead(151, 2, 160)},
                         std::vector<ExpectedCycleInfosForRead>{ExpectedCycleInfosForRead(151, 0, 1),
                                                                ExpectedCycleInfosForRead(8, 1, 152),
                                                                ExpectedCycleInfosForRead(151, 2, 160)},
                         std::vector<bcl2fastq::common::CycleRange>{bcl2fastq::common::CycleRange(0, 0),
                                                                    bcl2fastq::common::CycleRange(0, 0),
                                                                    bcl2fastq::common::CycleRange(0, 0)},
                         std::vector<bcl2fastq::common::CycleNumber>{0, 0, 0});
}

void TestLayout::testUmiCyclesAtBeginningTrim()
{
    bcl2fastq::config::BarcodeMismatchCountsContainer barcodeMismatches;

    bcl2fastq::layout::LaneInfosContainer laneInfos;
    getLaneInfos(laneInfos);

    bcl2fastq::layout::detectReadLayout(
        *runInfoXml_,
        *sampleSheet_, // used for adapters only
        std::vector<std::string>(), // use bases masks
        35,
        0, // read1StartCycle
        0, // read2StartCycle
        0, // read1EndCycle
        0, // read2EndCycle
        10, // read1UmiLength
        12, // read2UmiLength
        1, // read1UmiStartFromCycle
        1, // read2UmiStartFromCycle
        true, // trimUmi
        false,
        barcodeMismatches,
        laneInfos,
        barcodeLengthsForLane_);

    CPPUNIT_ASSERT_EQUAL(laneInfos.back().getReadInfos().size(),3UL);
    assertExpectedCycles(laneInfos.back(),
                         std::vector<ExpectedCycleInfosForRead>{ExpectedCycleInfosForRead(141, 0, 11),
                                                                ExpectedCycleInfosForRead(8, 1, 152),
                                                                ExpectedCycleInfosForRead(139, 2, 172)},
                         std::vector<ExpectedCycleInfosForRead>{ExpectedCycleInfosForRead(151, 0, 1),
                                                                ExpectedCycleInfosForRead(8, 1, 152),
                                                                ExpectedCycleInfosForRead(151, 2, 160)},
                         std::vector<ExpectedCycleInfosForRead>{ExpectedCycleInfosForRead(151, 0, 1),
                                                                ExpectedCycleInfosForRead(8, 1, 152),
                                                                ExpectedCycleInfosForRead(151, 2, 160)},
                         std::vector<bcl2fastq::common::CycleRange>{bcl2fastq::common::CycleRange(0, 10),
                                                                    bcl2fastq::common::CycleRange(0, 0),
                                                                    bcl2fastq::common::CycleRange(0, 12)},
                         std::vector<bcl2fastq::common::CycleNumber>{10, 0, 12});
}

void TestLayout::testUmiCyclesAtEndTrim()
{
    bcl2fastq::config::BarcodeMismatchCountsContainer barcodeMismatches;

    bcl2fastq::layout::LaneInfosContainer laneInfos;
    getLaneInfos(laneInfos);

    bcl2fastq::layout::detectReadLayout(
        *runInfoXml_,
        *sampleSheet_, // used for adapters only
        std::vector<std::string>(), // use bases masks
        35,
        0, // read1StartCycle
        0, // read2StartCycle
        0, // read1EndCycle
        0, // read2EndCycle
        10, // read1UmiLength
        12, // read2UmiLength
        142, // read1UmiStartFromCycle
        140, // read2UmiStartFromCycle
        true, // trimUmi
        false,
        barcodeMismatches,
        laneInfos,
        barcodeLengthsForLane_);

    CPPUNIT_ASSERT_EQUAL(laneInfos.back().getReadInfos().size(),3UL);
    assertExpectedCycles(laneInfos.back(),
                         std::vector<ExpectedCycleInfosForRead>{ExpectedCycleInfosForRead(141, 0, 1),
                                                                ExpectedCycleInfosForRead(8, 1, 152),
                                                                ExpectedCycleInfosForRead(139, 2, 160)},
                         std::vector<ExpectedCycleInfosForRead>{ExpectedCycleInfosForRead(151, 0, 1),
                                                                ExpectedCycleInfosForRead(8, 1, 152),
                                                                ExpectedCycleInfosForRead(151, 2, 160)},
                         std::vector<ExpectedCycleInfosForRead>{ExpectedCycleInfosForRead(151, 0, 1),
                                                                ExpectedCycleInfosForRead(8, 1, 152),
                                                                ExpectedCycleInfosForRead(151, 2, 160)},
                         std::vector<bcl2fastq::common::CycleRange>{bcl2fastq::common::CycleRange(141, 151),
                                                                    bcl2fastq::common::CycleRange(0, 0),
                                                                    bcl2fastq::common::CycleRange(139, 151)},
                         std::vector<bcl2fastq::common::CycleNumber>{0, 0, 0});
}

void TestLayout::testUmiCyclesInMiddleTrim()
{
    bcl2fastq::config::BarcodeMismatchCountsContainer barcodeMismatches;

    bcl2fastq::layout::LaneInfosContainer laneInfos;
    getLaneInfos(laneInfos);

    CPPUNIT_ASSERT_THROW(
        bcl2fastq::layout::detectReadLayout(
            *runInfoXml_,
            *sampleSheet_, // used for adapters only
            std::vector<std::string>(), // use bases masks
            35,
            0, // read1StartCycle
            0, // read2StartCycle
            0, // read1EndCycle
            0, // read2EndCycle
            8, // read1UmiLength
            12, // read2UmiLength
            142, // read1UmiStartFromCycle
            140, // read2UmiStartFromCycle
            true, // trimUmi
            false,
            barcodeMismatches,
            laneInfos,
            barcodeLengthsForLane_),
        bcl2fastq::common::InputDataError);
}




void TestLayout::testUmiCyclesAtBeginningOfMaskedReadTrim()
{
    bcl2fastq::config::BarcodeMismatchCountsContainer barcodeMismatches;

    bcl2fastq::layout::LaneInfosContainer laneInfos;
    getLaneInfos(laneInfos);

    bcl2fastq::layout::detectReadLayout(
        *runInfoXml_,
        *sampleSheet_, // used for adapters only
        std::vector<std::string>(), // use bases masks
        35,
        50, // read1StartCycle
        4, // read2StartCycle
        0, // read1EndCycle
        0, // read2EndCycle
        2, // read1UmiLength
        12, // read2UmiLength
        2, // read1UmiStartFromCycle
        3, // read2UmiStartFromCycle
        true, // trimUmi
        false,
        barcodeMismatches,
        laneInfos,
        barcodeLengthsForLane_);

    CPPUNIT_ASSERT_EQUAL(laneInfos.back().getReadInfos().size(),3UL);
    assertExpectedCycles(laneInfos.back(),
                         std::vector<ExpectedCycleInfosForRead>{ExpectedCycleInfosForRead(102, 0, 50),
                                                                ExpectedCycleInfosForRead(8, 1, 152),
                                                                ExpectedCycleInfosForRead(137, 2, 174)},
                         std::vector<ExpectedCycleInfosForRead>{ExpectedCycleInfosForRead(151, 0, 1),
                                                                ExpectedCycleInfosForRead(8, 1, 152),
                                                                ExpectedCycleInfosForRead(151, 2, 160)},
                         std::vector<ExpectedCycleInfosForRead>{},//ExpectedCycleInfosForRead(104, 0, 2),
                                                                //ExpectedCycleInfosForRead(8, 1, 152),
                                                                //ExpectedCycleInfosForRead(149, 2, 162)},
                         std::vector<bcl2fastq::common::CycleRange>{bcl2fastq::common::CycleRange(0, 2),
                                                                    bcl2fastq::common::CycleRange(0, 0),
                                                                    bcl2fastq::common::CycleRange(0, 12)},
                         std::vector<bcl2fastq::common::CycleNumber>{2, 0, 12});

    CPPUNIT_ASSERT_EQUAL(laneInfos.back().getReadInfos().at(0).cyclesToLoad().size(), 104UL);
}

void TestLayout::testUmiCyclesAtEndOfMaskedReadTrim()
{
    bcl2fastq::config::BarcodeMismatchCountsContainer barcodeMismatches;

    bcl2fastq::layout::LaneInfosContainer laneInfos;
    getLaneInfos(laneInfos);

    bcl2fastq::layout::detectReadLayout(
        *runInfoXml_,
        *sampleSheet_, // used for adapters only
        std::vector<std::string>(), // use bases masks
        35,
        2, // read1StartCycle
        3, // read2StartCycle
        100, // read1EndCycle
        99, // read2EndCycle
        10, // read1UmiLength
        12, // read2UmiLength
        140, // read1UmiStartFromCycle
        138, // read2UmiStartFromCycle
        true, // trimUmi
        false,
        barcodeMismatches,
        laneInfos,
        barcodeLengthsForLane_);

    CPPUNIT_ASSERT_EQUAL(laneInfos.back().getReadInfos().size(),3UL);
    assertExpectedCycles(laneInfos.back(),
                         std::vector<ExpectedCycleInfosForRead>{ExpectedCycleInfosForRead(99, 0, 2),
                                                                ExpectedCycleInfosForRead(8, 1, 152),
                                                                ExpectedCycleInfosForRead(97, 2, 162)},
                         std::vector<ExpectedCycleInfosForRead>{ExpectedCycleInfosForRead(151, 0, 1),
                                                                ExpectedCycleInfosForRead(8, 1, 152),
                                                                ExpectedCycleInfosForRead(151, 2, 160)},
                         std::vector<ExpectedCycleInfosForRead>{},//ExpectedCycleInfosForRead(109, 0, 2),
                                                                //ExpectedCycleInfosForRead(8, 1, 152),
                                                                //ExpectedCycleInfosForRead(149, 2, 162)},
                         std::vector<bcl2fastq::common::CycleRange>{bcl2fastq::common::CycleRange(99, 109),
                                                                    bcl2fastq::common::CycleRange(0, 0),
                                                                    bcl2fastq::common::CycleRange(97, 109)},
                         std::vector<bcl2fastq::common::CycleNumber>{0, 0, 0});
}

