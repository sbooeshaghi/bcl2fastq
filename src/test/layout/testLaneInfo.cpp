// data/testLaneInfo.cpp

#include <vector>
#include <string>
#include <boost/format.hpp>
#include <iostream>
#include <ios>

#include "common/Debug.hh"
#include "common/SystemCompatibility.hh"
#include "common/Options.hh"
#include "config/Bcl2FastqOptions.hh"
#include "io/SyncFile.hh"
#include "common/Logger.hh"
#include "layout/LaneInfo.hh"
#include "common/Types.hh"
//#include "data/CycleBCIFile.hh"
//#include "data/AggregatedBclFileReader.hh"
//#include "data/RawBclBuffer.hh"
#include "gtest/gtest.h"
#include "../bgzfData.hh"

using namespace bcl2fastq;

//=========================================
//
// TESTING A
//
//=========================================

TEST(TestLaneInfo, testDefault) {

    common::LaneNumber laneNumber = 1;

    // LaneInfo
    layout::LaneInfo laneInfo(laneNumber);
    EXPECT_EQ(1, laneInfo.getNumber());
    EXPECT_EQ("L001", laneInfo.getDirName().string());
    EXPECT_EQ(laneInfo.sampleInfosEnd(), laneInfo.sampleInfosBegin());
    EXPECT_EQ(0, laneInfo.getSampleInfos().size());
    EXPECT_EQ(0, laneInfo.getTileInfos().size());
    EXPECT_EQ(0, laneInfo.getReadInfos().size());
    EXPECT_EQ(0, laneInfo.getMinimumTrimmedReadLength());
    EXPECT_EQ(0, laneInfo.getNumCyclesToLoad());
    EXPECT_EQ(0, laneInfo.getClustersCount());
    EXPECT_TRUE(laneInfo.haveClustersCount());
}

//TEST(TestLaneInfo, testReads) {
//
//    layout::LaneInfosContainer lanes;
//    lanes.push_back(laneInfo);
//    //TODO: Setup reads?
//
//    // Define tiles
//    int tileIndex = 0;
//    int skippedTiles = 0;
//    for (int i = 0; i < tileCount; ++i) {
//        laneInfo.addTile(layout::TileInfo(testTiles[tileIndex], tileIndex, skippedTiles));
//        ++tileIndex;
//    }
//}

