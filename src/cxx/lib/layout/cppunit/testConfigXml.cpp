/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file testConfigXml.cpp
 *
 * \brief Sample sheet helper cppunit test declarations.
 *
 * \author Marek Balint
 */


#include <string>

#include "RegistryName.hh"
#include "testConfigXml.hh"

#include "io/Xml.hh"
#include "layout/ConfigXml.hh"


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(TestConfigXml, registryName("ConfigXml"));


void TestConfigXml::setUp()
{
}

void TestConfigXml::tearDown()
{
}

void TestConfigXml::testHcs()
{
    const std::string configXmlData(
        "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n"
        "<BaseCallAnalysis xmlns:xsd=\"http://www.w3.org/2001/XMLSchema\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">\n"
        "  <Run Name=\"BaseCalls\">\n"
        "    <TileSelection>\n"
        "      <Lane Index=\"6\">\n"
        "        <Sample>s</Sample>\n"
        "        <Tile>6001</Tile>\n"
        "        <Tile>6002</Tile>\n"
        "      </Lane>\n"
        "      <Lane Index=\"5\">\n"
        "        <Sample>s</Sample>\n"
        "        <Tile>5001</Tile>\n"
        "        <Tile>5002</Tile>\n"
        "      </Lane>\n"
        "      <Lane Index=\"8\">\n"
        "        <Sample>s</Sample>\n"
        "        <Tile>8001</Tile>\n"
        "        <Tile>8002</Tile>\n"
        "      </Lane>\n"
        "      <Lane Index=\"7\">\n"
        "        <Sample>s</Sample>\n"
        "        <Tile>7001</Tile>\n"
        "        <Tile>7002</Tile>\n"
        "      </Lane>\n"
        //"      <Lane Index=\"2\">\n"
        //"        <Sample>s</Sample>\n"
        //"        <Tile>2001</Tile>\n"
        //"        <Tile>2002</Tile>\n"
        //"      </Lane>\n"
        "      <Lane Index=\"1\">\n"
        "        <Sample>s</Sample>\n"
        "        <Tile>1001</Tile>\n"
        "        <Tile>1002</Tile>\n"
        "      </Lane>\n"
        "      <Lane Index=\"4\">\n"
        "        <Sample>s</Sample>\n"
        //"        <Tile>4001</Tile>\n"
        //"        <Tile>4002</Tile>\n"
        "      </Lane>\n"
        "      <Lane Index=\"3\">\n"
        "        <Sample>s</Sample>\n"
        "        <Tile>3001</Tile>\n"
        "        <Tile>3002</Tile>\n"
        "      </Lane>\n"
        "    </TileSelection>\n"
        "  </Run>\n"
        "</BaseCallAnalysis>\n"
    );

    bcl2fastq::layout::ConfigXml configXml(bcl2fastq::io::parseXmlData(
        configXmlData.begin(),
        configXmlData.end()
    ));

    for (bcl2fastq::common::LaneNumber laneNumber = 1; laneNumber <= 8; ++laneNumber)
    {
        bcl2fastq::layout::ConfigXml::TileMetadataContainer::const_iterator tileMetadataIter = configXml.tileMetadataBegin(laneNumber);
        const bcl2fastq::layout::ConfigXml::TileMetadataContainer::const_iterator tileMetadataEnd = configXml.tileMetadataEnd(laneNumber);
        if ((laneNumber == 2) || (laneNumber == 4))
        {
            CPPUNIT_ASSERT(tileMetadataIter == tileMetadataEnd);
        }
        else
        {
            CPPUNIT_ASSERT(tileMetadataIter != tileMetadataEnd);
            CPPUNIT_ASSERT_EQUAL(tileMetadataIter->tileNumber_, bcl2fastq::common::TileNumber((laneNumber * 1000) + 1));
            ++tileMetadataIter;
            CPPUNIT_ASSERT(tileMetadataIter != tileMetadataEnd);
            CPPUNIT_ASSERT_EQUAL(tileMetadataIter->tileNumber_, bcl2fastq::common::TileNumber((laneNumber * 1000) + 2));
            ++tileMetadataIter;
            CPPUNIT_ASSERT(tileMetadataIter == tileMetadataEnd);
        }
    }
}


