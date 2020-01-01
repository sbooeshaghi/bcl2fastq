/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file testSampleSheetCsv.cpp
 *
 * \brief Sample sheet helper cppunit test declarations.
 *
 * \author Marek Balint
 */


#include <string>

#include "RegistryName.hh"
#include "testSampleSheetCsv.hh"

#include "common/CsvGrammar.hh"
#include "config/SampleSheetCsv.hh"

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(TestSampleSheetCsv, registryName("SampleSheetCsv"));


void TestSampleSheetCsv::setUp()
{
}

void TestSampleSheetCsv::tearDown()
{
}

void TestSampleSheetCsv::testNova()
{
    const std::string sampleSheetData(
        "[Header],,,,\n"
        "Investigator Name,Isabelle,,,\n"
        "Project Name,Nova,,,\n"
        "Experiment Name,Orbital death ray research volume LXXIV,,,\n"
        "Date,5/27/2025,,,\n"
        "Workflow,GenerateFASTQ,,,\n"
        ",,,,\n"
        "[UnknownSection1],,,,\n"
        "Key1,Value1,,,\n"
        "Key2,Value2,,,\n"
        ",,,,\n"
        ",,,,\n"
        "[Settings],,,,\n"
        "MaskAdapter,CGCGTATACGCGTATA,,,\n"
        "MaskAdapterRead2,GCGCATATGCGCATAT,,,\n"
        "TrimAdapter,TATACGCGTATACGCG,,,\n"
        "TrimAdapterRead2,ATATGCGCATATGCGC,,,\n"
        ",,,,\n"
        "[UnknownSection2],,,,\n"
        "Key1,Value1,,,\n"
        "Key2,Value2,,,\n"
        ",,,,\n"
        "[Data],,,,\n"
        "LaNe,SampleID,SampleName,project,unknown,index,index2\n"
        "1,Aa,aA,P1,+_)(*(*,AAAAAAAA,AAAAAAAA\n"
        "2,Cc,cC,P1,!@#$,CCCCCCCC,CCCCCCCC\n"
        "3+4,Gg,gG,P2,#$%,GGGGGGGG,GGGGGGGG\n"
        "5,Tt,tT,P2,#2s^,TTTTTTTT,TTTTTTTT\n"
        "7,Xx,xX,P3,~!@,ACGTACGT,ACGTACGT\n"
        ",,,,\n"
        "[UnknownSection3],,,,\n"
        "Key1,Value1,,,\n"
        "Key2,Value2,,,\n"
        ",,,,\n"
    );

    bcl2fastq::config::SampleSheetCsv sampleSheetCsv(bcl2fastq::common::parseCsvData(
        sampleSheetData.begin(),
        sampleSheetData.end()
    ),
                                                     boost::filesystem::path());

    bcl2fastq::config::SampleSheetCsv::AdaptersContainer::const_iterator maskAdaptersIter = sampleSheetCsv.maskAdaptersBegin();
    const bcl2fastq::config::SampleSheetCsv::AdaptersContainer::const_iterator maskAdaptersEnd = sampleSheetCsv.maskAdaptersEnd();
    CPPUNIT_ASSERT(maskAdaptersIter != maskAdaptersEnd);
    CPPUNIT_ASSERT_EQUAL(maskAdaptersIter->first, bcl2fastq::common::ReadNumber(0));
    CPPUNIT_ASSERT_EQUAL(maskAdaptersIter->second, std::string("CGCGTATACGCGTATA"));
    ++maskAdaptersIter;
    CPPUNIT_ASSERT_EQUAL(maskAdaptersIter->first, bcl2fastq::common::ReadNumber(2));
    CPPUNIT_ASSERT_EQUAL(maskAdaptersIter->second, std::string("GCGCATATGCGCATAT"));
    ++maskAdaptersIter;
    CPPUNIT_ASSERT(maskAdaptersIter == maskAdaptersEnd);

    bcl2fastq::config::SampleSheetCsv::AdaptersContainer::const_iterator trimAdaptersIter = sampleSheetCsv.trimAdaptersBegin();
    const bcl2fastq::config::SampleSheetCsv::AdaptersContainer::const_iterator trimAdaptersEnd = sampleSheetCsv.trimAdaptersEnd();
    CPPUNIT_ASSERT(trimAdaptersIter != trimAdaptersEnd);
    CPPUNIT_ASSERT_EQUAL(trimAdaptersIter->first, bcl2fastq::common::ReadNumber(0));
    CPPUNIT_ASSERT_EQUAL(trimAdaptersIter->second, std::string("TATACGCGTATACGCG"));
    ++trimAdaptersIter;
    CPPUNIT_ASSERT_EQUAL(trimAdaptersIter->first, bcl2fastq::common::ReadNumber(2));
    CPPUNIT_ASSERT_EQUAL(trimAdaptersIter->second, std::string("ATATGCGCATATGCGC"));
    ++trimAdaptersIter;
    CPPUNIT_ASSERT(trimAdaptersIter == trimAdaptersEnd);

    bcl2fastq::common::SampleMetadataContainer::const_iterator sampleMetadataIter = sampleSheetCsv.getSampleMetadata().begin();
    const bcl2fastq::common::SampleMetadataContainer::const_iterator sampleMetadataEnd = sampleSheetCsv.getSampleMetadata().end();
    CPPUNIT_ASSERT(sampleMetadataIter != sampleMetadataEnd);
    CPPUNIT_ASSERT_EQUAL(sampleMetadataIter->id_, std::string("Aa"));
    CPPUNIT_ASSERT_EQUAL(sampleMetadataIter->name_, std::string("aA"));
    CPPUNIT_ASSERT_EQUAL(sampleMetadataIter->project_, std::string("P1"));
    CPPUNIT_ASSERT_EQUAL(sampleMetadataIter->lanes_.size(), 1UL);
    CPPUNIT_ASSERT_EQUAL(sampleMetadataIter->lanes_[0], 1UL);
    {
        bcl2fastq::common::SampleMetadata::BarcodesContainer::const_iterator barcodesIter = sampleMetadataIter->barcodes_.begin();
        const bcl2fastq::common::SampleMetadata::BarcodesContainer::const_iterator barcodesEnd = sampleMetadataIter->barcodes_.end();
        CPPUNIT_ASSERT(barcodesIter != barcodesEnd);
        {
            bcl2fastq::common::SampleMetadata::BarcodesContainer::value_type::const_iterator componentsIter = barcodesIter->begin();
            const bcl2fastq::common::SampleMetadata::BarcodesContainer::value_type::const_iterator componentsEnd = barcodesIter->end();
            CPPUNIT_ASSERT(componentsIter != componentsEnd);
            CPPUNIT_ASSERT_EQUAL(*componentsIter, std::string("AAAAAAAA"));
            ++componentsIter;
            CPPUNIT_ASSERT(componentsIter != componentsEnd);
            CPPUNIT_ASSERT_EQUAL(*componentsIter, std::string("AAAAAAAA"));
            ++componentsIter;
            CPPUNIT_ASSERT(componentsIter == componentsEnd);
        }
    }
    ++sampleMetadataIter;
    CPPUNIT_ASSERT(sampleMetadataIter != sampleMetadataEnd);
    CPPUNIT_ASSERT_EQUAL(sampleMetadataIter->id_, std::string("Cc"));
    CPPUNIT_ASSERT_EQUAL(sampleMetadataIter->name_, std::string("cC"));
    CPPUNIT_ASSERT_EQUAL(sampleMetadataIter->project_, std::string("P1"));
    CPPUNIT_ASSERT_EQUAL(sampleMetadataIter->lanes_.size(), 1UL);
    CPPUNIT_ASSERT_EQUAL(sampleMetadataIter->lanes_[0], 2UL);
    {
        bcl2fastq::common::SampleMetadata::BarcodesContainer::const_iterator barcodesIter = sampleMetadataIter->barcodes_.begin();
        const bcl2fastq::common::SampleMetadata::BarcodesContainer::const_iterator barcodesEnd = sampleMetadataIter->barcodes_.end();
        CPPUNIT_ASSERT(barcodesIter != barcodesEnd);
        {
            bcl2fastq::common::SampleMetadata::BarcodesContainer::value_type::const_iterator componentsIter = barcodesIter->begin();
            const bcl2fastq::common::SampleMetadata::BarcodesContainer::value_type::const_iterator componentsEnd = barcodesIter->end();
            CPPUNIT_ASSERT(componentsIter != componentsEnd);
            CPPUNIT_ASSERT_EQUAL(*componentsIter, std::string("CCCCCCCC"));
            ++componentsIter;
            CPPUNIT_ASSERT(componentsIter != componentsEnd);
            CPPUNIT_ASSERT_EQUAL(*componentsIter, std::string("CCCCCCCC"));
            ++componentsIter;
            CPPUNIT_ASSERT(componentsIter == componentsEnd);
        }
    }
    ++sampleMetadataIter;
    CPPUNIT_ASSERT(sampleMetadataIter != sampleMetadataEnd);
    CPPUNIT_ASSERT_EQUAL(sampleMetadataIter->id_, std::string("Gg"));
    CPPUNIT_ASSERT_EQUAL(sampleMetadataIter->name_, std::string("gG"));
    CPPUNIT_ASSERT_EQUAL(sampleMetadataIter->project_, std::string("P2"));
    CPPUNIT_ASSERT_EQUAL(sampleMetadataIter->lanes_.size(), 2UL);
    CPPUNIT_ASSERT_EQUAL(sampleMetadataIter->lanes_[0], 3UL);
    CPPUNIT_ASSERT_EQUAL(sampleMetadataIter->lanes_[1], 4UL);
    {
        bcl2fastq::common::SampleMetadata::BarcodesContainer::const_iterator barcodesIter = sampleMetadataIter->barcodes_.begin();
        const bcl2fastq::common::SampleMetadata::BarcodesContainer::const_iterator barcodesEnd = sampleMetadataIter->barcodes_.end();
        CPPUNIT_ASSERT(barcodesIter != barcodesEnd);
        {
            bcl2fastq::common::SampleMetadata::BarcodesContainer::value_type::const_iterator componentsIter = barcodesIter->begin();
            const bcl2fastq::common::SampleMetadata::BarcodesContainer::value_type::const_iterator componentsEnd = barcodesIter->end();
            CPPUNIT_ASSERT(componentsIter != componentsEnd);
            CPPUNIT_ASSERT_EQUAL(*componentsIter, std::string("GGGGGGGG"));
            ++componentsIter;
            CPPUNIT_ASSERT(componentsIter != componentsEnd);
            CPPUNIT_ASSERT_EQUAL(*componentsIter, std::string("GGGGGGGG"));
            ++componentsIter;
            CPPUNIT_ASSERT(componentsIter == componentsEnd);
        }
    }
    ++sampleMetadataIter;
    CPPUNIT_ASSERT(sampleMetadataIter != sampleMetadataEnd);
    CPPUNIT_ASSERT_EQUAL(sampleMetadataIter->id_, std::string("Tt"));
    CPPUNIT_ASSERT_EQUAL(sampleMetadataIter->name_, std::string("tT"));
    CPPUNIT_ASSERT_EQUAL(sampleMetadataIter->project_, std::string("P2"));
    CPPUNIT_ASSERT_EQUAL(sampleMetadataIter->lanes_.size(), 1UL);
    CPPUNIT_ASSERT_EQUAL(sampleMetadataIter->lanes_[0], 5UL);
    {
        bcl2fastq::common::SampleMetadata::BarcodesContainer::const_iterator barcodesIter = sampleMetadataIter->barcodes_.begin();
        const bcl2fastq::common::SampleMetadata::BarcodesContainer::const_iterator barcodesEnd = sampleMetadataIter->barcodes_.end();
        CPPUNIT_ASSERT(barcodesIter != barcodesEnd);
        {
            bcl2fastq::common::SampleMetadata::BarcodesContainer::value_type::const_iterator componentsIter = barcodesIter->begin();
            const bcl2fastq::common::SampleMetadata::BarcodesContainer::value_type::const_iterator componentsEnd = barcodesIter->end();
            CPPUNIT_ASSERT(componentsIter != componentsEnd);
            CPPUNIT_ASSERT_EQUAL(*componentsIter, std::string("TTTTTTTT"));
            ++componentsIter;
            CPPUNIT_ASSERT(componentsIter != componentsEnd);
            CPPUNIT_ASSERT_EQUAL(*componentsIter, std::string("TTTTTTTT"));
            ++componentsIter;
            CPPUNIT_ASSERT(componentsIter == componentsEnd);
        }
    }
    ++sampleMetadataIter;
    CPPUNIT_ASSERT(sampleMetadataIter != sampleMetadataEnd);
    CPPUNIT_ASSERT_EQUAL(sampleMetadataIter->id_, std::string("Xx"));
    CPPUNIT_ASSERT_EQUAL(sampleMetadataIter->name_, std::string("xX"));
    CPPUNIT_ASSERT_EQUAL(sampleMetadataIter->project_, std::string("P3"));
    CPPUNIT_ASSERT_EQUAL(sampleMetadataIter->lanes_.size(), 1UL);
    CPPUNIT_ASSERT_EQUAL(sampleMetadataIter->lanes_[0], 7UL);
    {
        bcl2fastq::common::SampleMetadata::BarcodesContainer::const_iterator barcodesIter = sampleMetadataIter->barcodes_.begin();
        const bcl2fastq::common::SampleMetadata::BarcodesContainer::const_iterator barcodesEnd = sampleMetadataIter->barcodes_.end();
        CPPUNIT_ASSERT(barcodesIter != barcodesEnd);
        {
            bcl2fastq::common::SampleMetadata::BarcodesContainer::value_type::const_iterator componentsIter = barcodesIter->begin();
            const bcl2fastq::common::SampleMetadata::BarcodesContainer::value_type::const_iterator componentsEnd = barcodesIter->end();
            CPPUNIT_ASSERT(componentsIter != componentsEnd);
            CPPUNIT_ASSERT_EQUAL(*componentsIter, std::string("ACGTACGT"));
            ++componentsIter;
            CPPUNIT_ASSERT(componentsIter != componentsEnd);
            CPPUNIT_ASSERT_EQUAL(*componentsIter, std::string("ACGTACGT"));
            ++componentsIter;
            CPPUNIT_ASSERT(componentsIter == componentsEnd);
        }
    }
    ++sampleMetadataIter;
    CPPUNIT_ASSERT(sampleMetadataIter == sampleMetadataEnd);
}

void TestSampleSheetCsv::testInvalidCharacters()
{
    const std::string invalidSampleIdData(
        "[Data],,,,\n"
        "SampleID,SampleName,project,unknown,index,index2\n"
        "Aa@,aA,P1,+_)(*(*,AAAAAAAA,AAAAAAAA\n"
        ",,,,\n"
    );

    CPPUNIT_ASSERT_THROW(bcl2fastq::config::SampleSheetCsv sampleSheetCsv(bcl2fastq::common::parseCsvData(
                         invalidSampleIdData.begin(),
                         invalidSampleIdData.end()),
                                                                          boost::filesystem::path()),
                         bcl2fastq::common::InputDataError);

    const std::string invalidSampleNameData(
        "[Data],,,,\n"
        "SampleID,SampleName,project,unknown,index,index2\n"
        "Aa,aA&,P1,+_)(*(*,AAAAAAAA,AAAAAAAA\n"
        ",,,,\n"
    );

    CPPUNIT_ASSERT_THROW(bcl2fastq::config::SampleSheetCsv sampleSheetCsv(bcl2fastq::common::parseCsvData(
                         invalidSampleNameData.begin(),
                         invalidSampleNameData.end()),
                                                                          boost::filesystem::path()),
                         bcl2fastq::common::InputDataError);

    const std::string invalidProjectData(
        "[Data],,,,\n"
        "SampleID,SampleName,project,unknown,index,index2\n"
        "Aa,aA,P1),+_)(*(*,AAAAAAAA,AAAAAAAA\n"
        ",,,,\n"
    );

    CPPUNIT_ASSERT_THROW(bcl2fastq::config::SampleSheetCsv sampleSheetCsv(bcl2fastq::common::parseCsvData(
                         invalidProjectData.begin(),
                         invalidProjectData.end()),
                                                                          boost::filesystem::path()),
                         bcl2fastq::common::InputDataError);

    const std::string invalidLaneData(
        "[Data],,,,\n"
        "Lane,SampleID,SampleName,project,unknown,index,index2\n"
        "1;2;3,Aa,aA,P1,+_)(*(*,AAAAAAAA,AAAAAAAA\n"
        ",,,,\n"
    );

    CPPUNIT_ASSERT_THROW(bcl2fastq::config::SampleSheetCsv sampleSheetCsv(bcl2fastq::common::parseCsvData(
                         invalidLaneData.begin(),
                         invalidLaneData.end()),
                                                                          boost::filesystem::path()),
                         bcl2fastq::common::InputDataError);
}

void TestSampleSheetCsv::testDuplicateBarcode()
{
    const std::string duplicateBarcodeAltLaneData(
        "[Data],,,,\n"
        "Lane,SampleID,SampleName,project,unknown,index,index2\n"
        "1+2+3,Aa,aA,P1,+_)(*(*,AAAAAAAA,AAAAAAAA\n"
        "4,Aa,aA,P1,+_)(*(*,AAAAAAAA,AAAAAAAA\n"
        ",,,,\n"
    );

    CPPUNIT_ASSERT_NO_THROW(bcl2fastq::config::SampleSheetCsv sampleSheetCsv(
        bcl2fastq::common::parseCsvData(
            duplicateBarcodeAltLaneData.begin(),
            duplicateBarcodeAltLaneData.end()),
        boost::filesystem::path()));

    const std::string duplicateBarcodeSameLaneData(
        "[Data],,,,\n"
        "Lane,SampleID,SampleName,project,unknown,index,index2\n"
        "1+2+3,Aa,aA,P1,+_)(*(*,AAAAAAAA,AAAAAAAA\n"
        "1,Bb,bB,P2,+_)(*(*,AAAAAAAA,AAAAAAAA\n"
        ",,,,\n"
    );

    CPPUNIT_ASSERT_THROW(
        bcl2fastq::config::SampleSheetCsv sampleSheetCsv(
            bcl2fastq::common::parseCsvData(duplicateBarcodeSameLaneData.begin(),
                                            duplicateBarcodeSameLaneData.end()),
            boost::filesystem::path()),
        bcl2fastq::common::InputDataError);

}
