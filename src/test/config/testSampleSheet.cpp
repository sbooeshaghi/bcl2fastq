// testSampleSheet

#include <vector>
#include <string>
#include <boost/format.hpp>
#include <iostream>
#include <ios>

#include "common/Debug.hh"
#include "common/SystemCompatibility.hh"
#include "common/CsvGrammar.hh"
#include "config/SampleSheetCsv.hh"
#include "io/SyncFile.hh"

#include "data/InteropFile.hh"
#include "gtest/gtest.h"

using namespace bcl2fastq;

// TODO: Move this to helper class - elsewhere
namespace {

    config::SampleSheetCsvPtr createSampleSheet(
        std::string& sampleSheetData) {

        return std::make_shared<config::SampleSheetCsv>(
            common::parseCsvData(
                sampleSheetData.begin(), sampleSheetData.end()),
                boost::filesystem::path());
    }

    config::SampleSheetCsvPtr createSampleSheet(
        std::vector< std::pair< std::string, std::string > >& settings) {

        std::string sampleSheetData = "[Settings]\n";
        for( auto keyValue : settings) {
            sampleSheetData += keyValue.first + "," + keyValue.second + "\n";
        }

        return std::make_shared<config::SampleSheetCsv>(
            common::parseCsvData(
                sampleSheetData.begin(), sampleSheetData.end()),
                boost::filesystem::path());
    }

    config::SampleSheetCsvPtr createAdapterSampleSheet(
        const std::string &trimAdapter1,
        const std::string &trimAdapter2) {

        std::vector< std::pair< std::string, std::string > > settings;
        settings.push_back( std::make_pair("TrimAdapter", trimAdapter1) );
        settings.push_back( std::make_pair("TrimAdapterRead2", trimAdapter2) );

        return createSampleSheet(settings);
    }

    config::SampleSheetCsvPtr createMaskAdapterSampleSheet(
        const std::string &maskAdapter1,
        const std::string &maskAdapter2) {

        std::vector< std::pair< std::string, std::string > > settings;
        settings.push_back( std::make_pair("MaskAdapter", maskAdapter1) );
        settings.push_back( std::make_pair("MaskAdapterRead2", maskAdapter2) );

        return createSampleSheet(settings);
    }

} // end unnamed namespace


//=========================================
//
// TESTING Empty Sample Sheet
//
//=========================================

TEST(SampleSheetTest, EmptySampleSheet) {

    // An empty string
    std::string sampleSheetData = "";
    config::SampleSheetCsvPtr sampleSheetCsv =
        std::make_shared<config::SampleSheetCsv>(
            common::parseCsvData(
                sampleSheetData.begin(), sampleSheetData.end()),
                boost::filesystem::path());

    EXPECT_EQ(config::SampleSheetCsv::INDETERMINATE, sampleSheetCsv->findAdaptersWithIndels());

    // "[Settings]," with comma only
    sampleSheetData = "[Settings],\n";
    sampleSheetCsv =
        std::make_shared<config::SampleSheetCsv>(
            common::parseCsvData(
                sampleSheetData.begin(), sampleSheetData.end()),
                boost::filesystem::path());

    EXPECT_EQ(config::SampleSheetCsv::INDETERMINATE, sampleSheetCsv->findAdaptersWithIndels());

// THIS THROW EXCEPTION - If the above cases are valid, this should be too or visa versa
//    // "[Settings]" without comma only
//    sampleSheetData = "[Settings]\n";
//    sampleSheetCsv =
//        std::make_shared<config::SampleSheetCsv>(
//            common::parseCsvData(
//                sampleSheetData.begin(), sampleSheetData.end()),
//                boost::filesystem::path());
//
//    EXPECT_EQ(config::SampleSheetCsv::INDETERMINATE, sampleSheetCsv->findAdaptersWithIndels());

}

//=========================================
//
// TESTING TrimAdapter
//
//=========================================

TEST(SampleSheetTest, DefineTrimAdapter) {

    const std::string trimAdapter1 = "TATACGCGTATACGCG";
    const std::string trimAdapter2 = "ATATGCGCATATGCGC";
    const std::string empty;

    // Define both trim adapters
    {
        std::vector< std::pair< std::string, std::string > > settings;
        settings.push_back( std::make_pair("TrimAdapter", trimAdapter1) );
        settings.push_back( std::make_pair("TrimAdapterRead2", trimAdapter2) );
        auto sampleSheetCsv = createSampleSheet(settings);

        config::SampleSheetCsv::AdaptersContainer::const_iterator it = sampleSheetCsv->trimAdaptersBegin();

        EXPECT_EQ(2, sampleSheetCsv->trimAdaptersEnd() - it);
        EXPECT_EQ(0,it->first);
        EXPECT_EQ(trimAdapter1,it->second);
        EXPECT_EQ(2, (++it)->first);
        EXPECT_EQ(trimAdapter2, it->second);
    }

    // Define both trim adapters (alt definition)
    {
        std::vector< std::pair< std::string, std::string > > settings;
        settings.push_back( std::make_pair("Adapter", trimAdapter1) );
        settings.push_back( std::make_pair("AdapterRead2", trimAdapter2) );
        auto sampleSheetCsv = createSampleSheet(settings);

        config::SampleSheetCsv::AdaptersContainer::const_iterator it = sampleSheetCsv->trimAdaptersBegin();

        EXPECT_EQ(2, sampleSheetCsv->trimAdaptersEnd() -it);
        EXPECT_EQ(0,it->first);
        EXPECT_EQ(trimAdapter1,it->second);
        EXPECT_EQ(2, (++it)->first);
        EXPECT_EQ(trimAdapter2, it->second);
    }

    // Only the first trim adapter
    {
        auto sampleSheetCsv = createAdapterSampleSheet( trimAdapter1, empty);
        config::SampleSheetCsv::AdaptersContainer::const_iterator it = sampleSheetCsv->trimAdaptersBegin();

        EXPECT_EQ(1, sampleSheetCsv->trimAdaptersEnd()-it);
        EXPECT_EQ(0,it->first);
        EXPECT_EQ(trimAdapter1,it->second);
    }

    // Only the second trim adapter ---- IS THIS WRONG?
    {
        auto sampleSheetCsv = createAdapterSampleSheet( empty, trimAdapter2);
        config::SampleSheetCsv::AdaptersContainer::const_iterator it = sampleSheetCsv->trimAdaptersBegin();

        EXPECT_EQ(1, sampleSheetCsv->trimAdaptersEnd()-it);
        EXPECT_EQ(2,it->first);
        EXPECT_EQ(trimAdapter2,it->second);
    }

    // Three trim adapters ---- IS THIS WRONG?
    {
        std::vector< std::pair< std::string, std::string > > settings;
        settings.push_back( std::make_pair("TrimAdapter", trimAdapter1) );
        settings.push_back( std::make_pair("TrimAdapterRead2", trimAdapter2) );
        settings.push_back( std::make_pair("TrimAdapterRead3", "ACGTACGTACGTACGT") );
        auto sampleSheetCsv = createSampleSheet(settings);

        config::SampleSheetCsv::AdaptersContainer::const_iterator it = sampleSheetCsv->trimAdaptersBegin();

        EXPECT_EQ(3, sampleSheetCsv->trimAdaptersEnd()-it);
        it += 2;
        EXPECT_EQ(3,it->first);
        EXPECT_EQ("ACGTACGTACGTACGT",it->second);
    }
}

TEST(SampleSheetTest, TestInvalidTrimAdapter) {

    const std::string badAdapter1  = "ZSDFGHJKASDFGHHJ";
    const std::string badAdapter2  = "YUUYUHFGFCVBHFGH";
    const std::string empty;

    EXPECT_THROW(auto sampleSheetCsv = createAdapterSampleSheet( badAdapter1, empty), common::InputDataError);
    EXPECT_THROW(auto sampleSheetCsv = createAdapterSampleSheet( badAdapter1, badAdapter2), common::InputDataError);
    EXPECT_THROW(auto sampleSheetCsv = createAdapterSampleSheet( "ACTGTGCAGT", badAdapter2), common::InputDataError);
//    EXPECT_THROW(auto sampleSheetCsv = createAdapterSampleSheet( empty, "ACTGTGCAGT"), common::InputDataError);
}

//=========================================
//
// TESTING MaskAdapter
//
//=========================================

TEST(SampleSheetTest, DefineMaskAdapter) {

    const std::string maskAdapter1 = "TATACGCGTATACGCG";
    const std::string maskAdapter2 = "ATATGCGCATATGCGC";
    const std::string empty;

    // Define both mask adapters
    {
        std::vector< std::pair< std::string, std::string > > settings;
        settings.push_back( std::make_pair("MaskAdapter", maskAdapter1) );
        settings.push_back( std::make_pair("MaskAdapterRead2", maskAdapter2) );
        auto sampleSheetCsv = createSampleSheet(settings);

        config::SampleSheetCsv::AdaptersContainer::const_iterator it = sampleSheetCsv->maskAdaptersBegin();

        EXPECT_EQ(2, sampleSheetCsv->maskAdaptersEnd() - it);
        EXPECT_EQ(0,it->first);
        EXPECT_EQ(maskAdapter1,it->second);
        EXPECT_EQ(2, (++it)->first);
        EXPECT_EQ(maskAdapter2, it->second);
    }

    // Only the first mask adapter
    {
        auto sampleSheetCsv = createMaskAdapterSampleSheet( maskAdapter1, empty);
        config::SampleSheetCsv::AdaptersContainer::const_iterator it = sampleSheetCsv->maskAdaptersBegin();

        EXPECT_EQ(1, sampleSheetCsv->maskAdaptersEnd()-it);
        EXPECT_EQ(0,it->first);
        EXPECT_EQ(maskAdapter1,it->second);
    }

    // Only the second mask adapter ---- IS THIS WRONG?
    {
        auto sampleSheetCsv = createMaskAdapterSampleSheet( empty, maskAdapter2);
        config::SampleSheetCsv::AdaptersContainer::const_iterator it = sampleSheetCsv->maskAdaptersBegin();

        EXPECT_EQ(1, sampleSheetCsv->maskAdaptersEnd()-it);
        EXPECT_EQ(2,it->first);
        EXPECT_EQ(maskAdapter2,it->second);
    }

    // Three mask adapters ---- IS THIS WRONG?
    {
        std::vector< std::pair< std::string, std::string > > settings;
        settings.push_back( std::make_pair("MaskAdapter", maskAdapter1) );
        settings.push_back( std::make_pair("MaskAdapterRead2", maskAdapter2) );
        settings.push_back( std::make_pair("MaskAdapterRead3", "ACGTACGTACGTACGT") );
        auto sampleSheetCsv = createSampleSheet(settings);

        config::SampleSheetCsv::AdaptersContainer::const_iterator it = sampleSheetCsv->maskAdaptersBegin();

        EXPECT_EQ(3, sampleSheetCsv->maskAdaptersEnd()-it);
        it += 2;
        EXPECT_EQ(3,it->first);
        EXPECT_EQ("ACGTACGTACGTACGT",it->second);
    }
}

TEST(SampleSheetTest, TestInvalidMaskAdapter) {

    const std::string badAdapter1  = "ZSDFGHJKASDFGHHJ";
    const std::string badAdapter2  = "YUUYUHFGFCVBHFGH";
    const std::string empty;

    EXPECT_THROW(auto sampleSheetCsv = createMaskAdapterSampleSheet( badAdapter1, empty), common::InputDataError);
    EXPECT_THROW(auto sampleSheetCsv = createMaskAdapterSampleSheet( badAdapter1, badAdapter2), common::InputDataError);
    EXPECT_THROW(auto sampleSheetCsv = createMaskAdapterSampleSheet( "ACTGTGCAGT", badAdapter2), common::InputDataError);
//    EXPECT_THROW(auto sampleSheetCsv = createMaskAdapterSampleSheet( empty, "ACTGTGCAGT"), common::InputDataError);
}

//=========================================
//
// TESTING findAdaptersWithIndels
// TESTING findAdapterWithIndels
//
//=========================================

namespace {

    config::SampleSheetCsvPtr createFindAdaptersWithIndelsData(
        const std::string &trueFalseValue,
        const std::string &key = "findAdaptersWithIndels") {

        std::vector< std::pair< std::string, std::string > > settings;
        settings.push_back( std::make_pair(key, trueFalseValue) );

        return createSampleSheet(settings);
    }

} // end unnamed namespace

// TEST FindAdapterWithIndels
TEST(SampleSheetTest, FindAdapterWithIndels) {

    // Undefined
    EXPECT_EQ(config::SampleSheetCsv::INDETERMINATE, createFindAdaptersWithIndelsData("ACGT", "MaskAdapter2")->findAdaptersWithIndels());

    // True
    EXPECT_EQ(config::SampleSheetCsv::TRUE, createFindAdaptersWithIndelsData("y")->findAdaptersWithIndels());
    EXPECT_EQ(config::SampleSheetCsv::TRUE, createFindAdaptersWithIndelsData("yEs")->findAdaptersWithIndels());
    EXPECT_EQ(config::SampleSheetCsv::TRUE, createFindAdaptersWithIndelsData("true")->findAdaptersWithIndels());
    EXPECT_EQ(config::SampleSheetCsv::TRUE, createFindAdaptersWithIndelsData("t")->findAdaptersWithIndels());
    EXPECT_EQ(config::SampleSheetCsv::TRUE, createFindAdaptersWithIndelsData("1")->findAdaptersWithIndels());

    // deprecated spelling
    EXPECT_EQ(config::SampleSheetCsv::TRUE, createFindAdaptersWithIndelsData("y", "findadapterwithindels")->findAdaptersWithIndels());

    // False
    EXPECT_EQ(config::SampleSheetCsv::FALSE, createFindAdaptersWithIndelsData("n")->findAdaptersWithIndels());
    EXPECT_EQ(config::SampleSheetCsv::FALSE, createFindAdaptersWithIndelsData("no")->findAdaptersWithIndels());
    EXPECT_EQ(config::SampleSheetCsv::FALSE, createFindAdaptersWithIndelsData("FalSE")->findAdaptersWithIndels());
    EXPECT_EQ(config::SampleSheetCsv::FALSE, createFindAdaptersWithIndelsData("f")->findAdaptersWithIndels());
    EXPECT_EQ(config::SampleSheetCsv::FALSE, createFindAdaptersWithIndelsData("0")->findAdaptersWithIndels());
}

// TEST FindAdapterWithIndels
TEST(SampleSheetTest, BadFindAdapterWithIndels) {

    // Undefined
    EXPECT_THROW(createFindAdaptersWithIndelsData("duh")->findAdaptersWithIndels(), common::InputDataError);
    EXPECT_THROW(createFindAdaptersWithIndelsData("ye")->findAdaptersWithIndels(), common::InputDataError);
    EXPECT_THROW(createFindAdaptersWithIndelsData("fals")->findAdaptersWithIndels(), common::InputDataError);

}

//=========================================
//
// TESTING Data/Lane
//
//=========================================

TEST(SampleSheetTest, DataContent) {

    // Test sample table header names ....
    std::string sampleSheetData =
        "[Data]\n"
        "SampleID,SampleName,project,Lane,Index\n"
        "SampleID_TEST1,Sample2,Project,2,CCCCCCCC\n"
        "1,Sample1,Project,1,AAAAAAAA\n";

    config::SampleSheetCsvPtr sampleSheetCsv = createSampleSheet(sampleSheetData);
    common::SampleMetadataContainer::const_iterator sampleIt = sampleSheetCsv->getSampleMetadata().begin();

    EXPECT_EQ(2, sampleSheetCsv->getSampleMetadata().end()-sampleIt);
    EXPECT_EQ("SampleID_TEST1", sampleIt->id_);
    EXPECT_EQ("Sample2", sampleIt->name_);
    EXPECT_EQ("Project", sampleIt->project_);
    EXPECT_EQ(2, sampleIt->lanes_[0]);
    EXPECT_EQ("CCCCCCCC", sampleIt->barcodes_[0][0]);

    // Test Alternate header names ....
    sampleSheetData =
        "[Data]\n"
        "Sample_ID,Sample_Name,sample_project,IGNORED,LaneS,Index2\n"
        "ID_TEST2,Sample2_TEST,Project2,thing1,2,CCCCCCCC\n"
        "1,Sample1,Project,thing2,1,AAAAAAAA\n";

    sampleSheetCsv = createSampleSheet(sampleSheetData);
    sampleIt = sampleSheetCsv->getSampleMetadata().begin();

    EXPECT_EQ(2, sampleSheetCsv->getSampleMetadata().end()-sampleIt);
    EXPECT_EQ("ID_TEST2", sampleIt->id_);
    EXPECT_EQ("Sample2_TEST", sampleIt->name_);
    EXPECT_EQ("Project2", sampleIt->project_);
    EXPECT_EQ(2, sampleIt->lanes_[0]);
    EXPECT_EQ(1, sampleIt->barcodes_.size());
    EXPECT_EQ("CCCCCCCC", sampleIt->barcodes_[0][1]);

    // Test white space removed
    const std::string id1 = " ID1 ";
    sampleSheetData =
        "[Data]\n"
        "Sample_ID,Sample_Name,sample_project,IGNORED,LaneS,Index2\n"
        + id1 + ",Sample2_TEST,Project2,thing1,2,CCCCCCCC\n"
        "1,Sample1,Project,thing2,1,AAAAAAAA\n";

    sampleSheetCsv = createSampleSheet(sampleSheetData);
    sampleIt = sampleSheetCsv->getSampleMetadata().begin();

    EXPECT_EQ(2, sampleSheetCsv->getSampleMetadata().end()-sampleIt);
    EXPECT_EQ("ID1", sampleIt->id_);

//    // Isolated \r is treated as \n -- NASTY!
//    const std::string id3 = "\rID3\r";
//    sampleSheetData =
//        "[Data]\n"
//        "Sample_ID,Sample_Name,sample_project,IGNORED,LaneS,Index2\n"
//        + id3 + ",Sample2_TEST,Project2,thing1,2,CCCCCCCC\n"
//        "1,Sample1,Project,thing2,1,AAAAAAAA\n";
//
//    sampleSheetCsv = createSampleSheet(sampleSheetData);
//    sampleIt = sampleSheetCsv->getSampleMetadata().begin();
//
//    EXPECT_EQ(2, sampleSheetCsv->getSampleMetadata().end()-sampleIt);
//    EXPECT_EQ("ID3", sampleIt->id_);

}

TEST(SampleSheetTest, InvalidDataContent) {

    // Valid chars [a-zA-Z0-9_-]+
    const std::string invalidText = "A%ID";
    const std::string invalidHiddenText = "\x07ID";
    const std::string invalidTabText = "\tID\t";

    // Invalid visible text
    std::string sampleSheetData =
        "[Data]\n"
        "SampleID,SampleName,project,Lanes,Index\n"
        "SampleID_TEST1,Sample2,Project,2,CCCCCCCC\n"
        + invalidText + ",Sample1,Project,1,AAAAAAAA\n";

    EXPECT_THROW(config::SampleSheetCsvPtr sampleSheetCsv = createSampleSheet(sampleSheetData), common::InputDataError);

    // Invalid hidden text (NASTY)
    sampleSheetData =
        "[Data]\n"
        "SampleID,SampleName,project,Lanes,Index\n"
        "SampleID_TEST1,Sample2,Project,2,CCCCCCCC\n"
        + invalidHiddenText + ",Sample1,Project,1,AAAAAAAA\n";

    EXPECT_THROW(config::SampleSheetCsvPtr sampleSheetCsv = createSampleSheet(sampleSheetData), common::InputDataError);

    // Invalid tab text (NASTY)
    sampleSheetData =
        "[Data]\n"
        "SampleID,SampleName,project,Lanes,Index\n"
        "SampleID_TEST1,Sample2,Project,2,CCCCCCCC\n"
        + invalidTabText + ",Sample1,Project,1,AAAAAAAA\n";

    EXPECT_THROW(config::SampleSheetCsvPtr sampleSheetCsv = createSampleSheet(sampleSheetData), common::InputDataError);
}

TEST(SampleSheetTest, ReadCycle) {

    //
    std::string sampleSheetData =
        "[Settings]\n"
        "Read1StartFromCycle,5\n"  // Read1StartWITHcyle silently doesn't startFROMcycle
        "Read1EndWithCycle,307\n"
        "\n";

    config::SampleSheetCsvPtr sampleSheetCsv = createSampleSheet(sampleSheetData);
    EXPECT_EQ(  5, sampleSheetCsv->getRead1StartCycle());
    EXPECT_EQ(307, sampleSheetCsv->getRead1EndCycle());

    //
    sampleSheetData =
        "[Settings]\n"
        "Read2StartFromCycle,6\n"  // Read1StartWITHcyle silently doesn't startFROMcycle
        "Read2EndWithCycle,306\n"
        "\n";

    sampleSheetCsv = createSampleSheet(sampleSheetData);
    EXPECT_EQ(  6, sampleSheetCsv->getRead2StartCycle());
    EXPECT_EQ(306, sampleSheetCsv->getRead2EndCycle());

}

//=========================================
//
// TESTING ReverseCompliment Setting
//
//=========================================

TEST(SampleSheetTest, ReverseComplement) {

    //
    std::string sampleSheetData =
        "[Settings]\n"
        "Things,1\n"
        "\n";

    config::SampleSheetCsvPtr sampleSheetCsv = createSampleSheet(sampleSheetData);
    EXPECT_EQ(config::SampleSheetCsv::INDETERMINATE, sampleSheetCsv->generateReverseComplementFastqs());

    sampleSheetData =
        "[Settings]\n"
        "ReverseComplement,1\n"
        "\n";

    sampleSheetCsv = createSampleSheet(sampleSheetData);
    EXPECT_EQ(config::SampleSheetCsv::TRUE, sampleSheetCsv->generateReverseComplementFastqs());

    sampleSheetData =
        "[Settings]\n"
        "ReverseComplement,0\n"
        "\n";

    sampleSheetCsv = createSampleSheet(sampleSheetData);
    EXPECT_EQ(config::SampleSheetCsv::FALSE, sampleSheetCsv->generateReverseComplementFastqs());
}
