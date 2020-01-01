// data/testLayout.cpp

#include "layout/Layout.hh"
#include "gtest/gtest.h"
#include "../bgzfData.hh"
#include "../helper/Bcl2FastqTestHelper.hh"

using namespace bcl2fastq;

namespace
{

    static const std::string inputPathString("testData/cbcl/btof1110/Data/Intensities/BaseCalls/");

    class TileLayoutDetectorAccessor : public layout::TileLayoutDetector
    {
        // Provides access to protected members of class.
    public:
        static void access_readTilesFromCbclHeaders(
            const boost::filesystem::path&   inputDir,
            const config::SampleSheetCsv&    sampleSheet,
            const std::vector<boost::regex>& tileRegexps,
            const std::string&               flowcellId,
            std::vector<layout::LaneInfo>&   lanes,
            common::TileFileMap&             tileFileMap,
            common::NumBasesPerByte&         numBasesPerByte,
            common::CycleNumber              read1StartCycle,
            common::CycleNumber              read1EndCycle,
            bool                             ignoreMissingBcls)
        {
            layout::TileLayoutDetector::readTilesFromCbclHeaders(
                inputDir,
                sampleSheet,
                tileRegexps,
                flowcellId,
                lanes,
                tileFileMap,
                numBasesPerByte,
                read1StartCycle,
                read1EndCycle,
                ignoreMissingBcls);
        }
    }; // end class TileLayoutDetectorAccessor

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

    config::SampleSheetCsvPtr createMaskAdapterSampleSheet(
        const std::string &maskAdapter1,
        const std::string &maskAdapter2) {

        std::vector< std::pair< std::string, std::string > > settings;
        settings.push_back( std::make_pair("MaskAdapter", maskAdapter1) );
        settings.push_back( std::make_pair("MaskAdapterRead2", maskAdapter2) );

        return createSampleSheet(settings);
    }

    // TODO: Add helper::SampleSheet-orator
    // Whatever just need a sample sheet without ExcludeTiles
    config::SampleSheetCsvPtr getASampleSheet() {
        const std::string trimAdapter1 = "TATACGCGTATACGCG";
        const std::string trimAdapter2 = "ATATGCGCATATGCGC";
        std::vector< std::pair< std::string, std::string > > settings;
        settings.push_back( std::make_pair("TrimAdapter", trimAdapter1) );
        settings.push_back( std::make_pair("TrimAdapterRead2", trimAdapter2) );
        return createSampleSheet(settings);
    }

} // end namespace

//=========================================
//
// TESTING A
//
//=========================================

// Three bad cycles throws exception on initialization.
TEST(TileLayoutDetectorTest, testRead3BadCycles) {

    helper::Bcl2FastqOptionsOratorPtr options = helper::Bcl2FastqOptionsOrator::createOptionsWithRunFolder("testData/cbcl/btof1110/testRead3BadCycles");

    const boost::filesystem::path   inputDir(options->options_.getInputDir());
    const std::vector<boost::regex> tileRegexps;
    const std::string               flowcellId;
    std::vector<layout::LaneInfo>   lanes;
    common::TileFileMap             tileFileMap;
    common::NumBasesPerByte         numBasesPerByte = common::NumBasesPerByte::TWO;
    bool                            ignoreMissingBcls = true;

    config::SampleSheetCsvPtr sampleSheetPtr = getASampleSheet();

    // Two lanes
    lanes.push_back(layout::LaneInfo(1));
    lanes.push_back(layout::LaneInfo(2));

    EXPECT_THROW(
        TileLayoutDetectorAccessor::access_readTilesFromCbclHeaders(
            inputDir,
            *sampleSheetPtr,
            tileRegexps,
            flowcellId,
            lanes,
            tileFileMap,
            numBasesPerByte,
            1,94,
            ignoreMissingBcls),
        common::InputDataError);
}

// Three bad cycles throws exception on initialization.
TEST(TileLayoutDetectorTest, testRead3BadCycles_startOn2) {

    helper::Bcl2FastqOptionsOratorPtr options = helper::Bcl2FastqOptionsOrator::createOptionsWithRunFolder("testData/cbcl/btof1110/testRead3BadCycles");

    const boost::filesystem::path   inputDir(options->options_.getInputDir());
    const std::vector<boost::regex> tileRegexps;
    const std::string               flowcellId;
    std::vector<layout::LaneInfo>   lanes;
    common::TileFileMap             tileFileMap;
    common::NumBasesPerByte         numBasesPerByte = common::NumBasesPerByte::TWO;
    bool                            ignoreMissingBcls = true;

    config::SampleSheetCsvPtr sampleSheetPtr = getASampleSheet();

    // Two lanes
    lanes.push_back(layout::LaneInfo(1));
    //lanes.push_back(layout::LaneInfo(2));

    EXPECT_THROW(
        TileLayoutDetectorAccessor::access_readTilesFromCbclHeaders(
            inputDir,
            *sampleSheetPtr,
            tileRegexps,
            flowcellId,
            lanes,
            tileFileMap,
            numBasesPerByte,
            2,94,
            ignoreMissingBcls),
        common::InputDataError);
}

// Cycle 1 on initialization uses cycle 2.
TEST(TileLayoutDetectorTest, testRead1BadCycle) {

    helper::Bcl2FastqOptionsOratorPtr options = helper::Bcl2FastqOptionsOrator::createOptionsWithRunFolder("testData/cbcl/btof1110/testRead1BadCycle");

    const boost::filesystem::path   inputDir(options->options_.getInputDir());
    const std::vector<boost::regex> tileRegexps;
    const std::string               flowcellId;
    std::vector<layout::LaneInfo>   lanes;
    common::TileFileMap             tileFileMap;
    common::NumBasesPerByte         numBasesPerByte = common::NumBasesPerByte::TWO;
    bool                            ignoreMissingBcls = true;

    config::SampleSheetCsvPtr sampleSheetPtr = getASampleSheet();

    // Two lanes
    lanes.push_back(layout::LaneInfo(1));
    lanes.push_back(layout::LaneInfo(2));

    TileLayoutDetectorAccessor::access_readTilesFromCbclHeaders(
        inputDir,
        *sampleSheetPtr,
        tileRegexps,
        flowcellId,
        lanes,
        tileFileMap,
        numBasesPerByte,
        1,94,
        ignoreMissingBcls);
}

// Start on cycle 2 
TEST(TileLayoutDetectorTest, testRead1BadCycles) {

    helper::Bcl2FastqOptionsOratorPtr options = helper::Bcl2FastqOptionsOrator::createOptionsWithRunFolder("testData/cbcl/btof1110/testRead1BadCycle");

    const boost::filesystem::path   inputDir(options->options_.getInputDir());
    const std::vector<boost::regex> tileRegexps;
    const std::string               flowcellId;
    std::vector<layout::LaneInfo>   lanes;
    common::TileFileMap             tileFileMap;
    common::NumBasesPerByte         numBasesPerByte = common::NumBasesPerByte::TWO;
    bool                            ignoreMissingBcls = true;

    config::SampleSheetCsvPtr sampleSheetPtr = getASampleSheet();

    // Two lanes
    lanes.push_back(layout::LaneInfo(1));
    //lanes.push_back(layout::LaneInfo(2));

    TileLayoutDetectorAccessor::access_readTilesFromCbclHeaders(
        inputDir,
        *sampleSheetPtr,
        tileRegexps,
        flowcellId,
        lanes,
        tileFileMap,
        numBasesPerByte,
        2,94,
        ignoreMissingBcls);
}


