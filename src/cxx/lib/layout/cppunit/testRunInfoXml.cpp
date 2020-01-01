/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file testRunInfoXml.cpp
 *
 * \brief Sample sheet helper cppunit test declarations.
 *
 * \author Marek Balint
 */


#include <string>

#include "RegistryName.hh"
#include "testRunInfoXml.hh"

#include "io/Xml.hh"
#include "layout/RunInfoXml.hh"


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(TestRunInfoXml, registryName("RunInfoXml"));


void TestRunInfoXml::setUp()
{
}

void TestRunInfoXml::tearDown()
{
}

void TestRunInfoXml::testNova()
{
    const std::string runInfoXmlData(
        "<?xml version=\"1.0\"?>\n"
        "<RunInfo xmlns:xsd=\"http://www.w3.org/2001/XMLSchema\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" Version=\"2\">\n"
        "  <Run Id=\"TestRun\" Number=\"42\">\n"
        "    <Flowcell>TestFlowcell</Flowcell>\n"
        "    <Instrument>CSSIM</Instrument>\n"
        "    <Date>8/20/2013 3:29:17 PM</Date>\n"
        "    <Reads>\n"
        "      <Read FirstCycle=\"1\" LastCycle=\"91\" />\n"
        "      <Read NumCycles=\"8\" IsIndexedRead=\"Y\" />\n"
        "      <Read FirstCycle=\"100\" LastCycle=\"107\" >\n"
        "        <Index />\n"
        "      </Read>\n"
        "      <Read FirstCycle=\"108\" LastCycle=\"200\" IsIndexedRead=\"N\" />\n"
        "    </Reads>\n"
        "    <FlowcellLayout LaneCount=\"4\" SurfaceCount=\"2\" SwathCount=\"3\" TileCount=\"11\" SectionPerLane=\"3\" LanePerSection=\"2\" />\n"
        "  </Run>\n"
        "</RunInfo>\n"
    );

    bcl2fastq::layout::RunInfoXml runInfoXml(bcl2fastq::io::parseXmlData(
        runInfoXmlData.begin(),
        runInfoXmlData.end()
    ));

    CPPUNIT_ASSERT_EQUAL(runInfoXml.getInstrument(), std::string("CSSIM"));
    CPPUNIT_ASSERT_EQUAL(runInfoXml.getRunNumber(), std::string("42"));
    CPPUNIT_ASSERT_EQUAL(runInfoXml.getFlowcellId(), std::string("TestFlowcell"));
    CPPUNIT_ASSERT_EQUAL(runInfoXml.getLanesCount(), bcl2fastq::common::LaneNumber(4));

    bcl2fastq::layout::ReadMetadataContainer::const_iterator readMetadataIter = runInfoXml.readMetadataBegin();
    const bcl2fastq::layout::ReadMetadataContainer::const_iterator readMetadataEnd = runInfoXml.readMetadataEnd();
    CPPUNIT_ASSERT(readMetadataIter != readMetadataEnd);
    CPPUNIT_ASSERT_EQUAL(readMetadataIter->readNumber_, bcl2fastq::common::ReadNumber(1));
    CPPUNIT_ASSERT_EQUAL(readMetadataIter->readType_, bcl2fastq::common::ReadType::DATA);
    CPPUNIT_ASSERT_EQUAL(readMetadataIter->firstCycle_, bcl2fastq::common::CycleNumber(1));
    CPPUNIT_ASSERT_EQUAL(readMetadataIter->lastCycle_, bcl2fastq::common::CycleNumber(91));
    ++readMetadataIter;
    CPPUNIT_ASSERT(readMetadataIter != readMetadataEnd);
    CPPUNIT_ASSERT_EQUAL(readMetadataIter->readNumber_, bcl2fastq::common::ReadNumber(1));
    CPPUNIT_ASSERT_EQUAL(readMetadataIter->readType_, bcl2fastq::common::ReadType::INDEX);
    CPPUNIT_ASSERT_EQUAL(readMetadataIter->firstCycle_, bcl2fastq::common::CycleNumber(92));
    CPPUNIT_ASSERT_EQUAL(readMetadataIter->lastCycle_, bcl2fastq::common::CycleNumber(99));
    ++readMetadataIter;
    CPPUNIT_ASSERT(readMetadataIter != readMetadataEnd);
    CPPUNIT_ASSERT_EQUAL(readMetadataIter->readNumber_, bcl2fastq::common::ReadNumber(2));
    CPPUNIT_ASSERT_EQUAL(readMetadataIter->readType_, bcl2fastq::common::ReadType::INDEX);
    CPPUNIT_ASSERT_EQUAL(readMetadataIter->firstCycle_, bcl2fastq::common::CycleNumber(100));
    CPPUNIT_ASSERT_EQUAL(readMetadataIter->lastCycle_, bcl2fastq::common::CycleNumber(107));
    ++readMetadataIter;
    CPPUNIT_ASSERT(readMetadataIter != readMetadataEnd);
    CPPUNIT_ASSERT_EQUAL(readMetadataIter->readNumber_, bcl2fastq::common::ReadNumber(2));
    CPPUNIT_ASSERT_EQUAL(readMetadataIter->readType_, bcl2fastq::common::ReadType::DATA);
    CPPUNIT_ASSERT_EQUAL(readMetadataIter->firstCycle_, bcl2fastq::common::CycleNumber(108));
    CPPUNIT_ASSERT_EQUAL(readMetadataIter->lastCycle_, bcl2fastq::common::CycleNumber(200));
    ++readMetadataIter;
    CPPUNIT_ASSERT(readMetadataIter == readMetadataEnd);
}


