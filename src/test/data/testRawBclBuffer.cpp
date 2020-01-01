// data/testCycleBCIFile.cpp

//#include <vector>
//#include <string>
//#include <boost/format.hpp>
//#include <iostream>
//#include <ios>
//
//#include "common/Debug.hh"
//#include "common/SystemCompatibility.hh"
//#include "common/Options.hh"
//#include "config/Bcl2FastqOptions.hh"
//#include "io/SyncFile.hh"
//#include "common/Logger.hh"
#include "data/RawBclBuffer.hh"
#include "gtest/gtest.h"

using namespace bcl2fastq;

TEST(TestRawDataBuffer,RawDataBuffer) {

    common::RawDataBuffer rawDataBuffer(10);

    EXPECT_EQ(10,rawDataBuffer.size());
}


// Clusters for a single cycle and tile.
TEST(TestPerCycle,PerCycleClass) {

    data::PerCycleData perCycleData;

    EXPECT_EQ(0, perCycleData.bcls_.size());
    EXPECT_TRUE(perCycleData.includeNonPf_);
    EXPECT_EQ(6, perCycleData.numBitsPerQscore_);
    EXPECT_EQ(0, perCycleData.remappedQscores_.size());

    perCycleData.bcls_.resize(32);
    EXPECT_EQ(1, perCycleData.getNumClustersPerByte());
    EXPECT_EQ(32, perCycleData.getNumClusters());

    perCycleData.numBitsPerQscore_ = 2;
    EXPECT_EQ(2, perCycleData.numBitsPerQscore_);
    EXPECT_EQ(2, perCycleData.getNumClustersPerByte());
    EXPECT_EQ(64, perCycleData.getNumClusters());
}

// All cycles and clusters for a single tile.
TEST(TestRawBclBuffer, RawBclBuffer) {

    data::RawBclBuffer rawBclBuffer;

    EXPECT_EQ(0, rawBclBuffer.getGroupNumber());
    rawBclBuffer.setGroupNumber(2);
    EXPECT_EQ(2, rawBclBuffer.getGroupNumber());

    EXPECT_FALSE(rawBclBuffer.hasCycleData());

    EXPECT_EQ( 0, rawBclBuffer.cycleData_.size());
    rawBclBuffer.cycleData_.resize(10);
    EXPECT_EQ(10, rawBclBuffer.cycleData_.size());
    EXPECT_EQ( 0, rawBclBuffer.cycleData_[4].bcls_.size());

    for (data::PerCycleData &cycleData : rawBclBuffer.cycleData_) {
        cycleData.bcls_.resize(16);
    }
    EXPECT_EQ(16, rawBclBuffer.cycleData_[2].bcls_.size());

    EXPECT_TRUE(rawBclBuffer.hasCycleData());

// gzipDecompressors_, filters_, positions_, tilOffsets, tileInfo_, uncompressedBclOffset_
}

// A thread safe sub-group of tiles.
TEST(TestRawBclBufferGroup, RawBclBufferGroup) {

    ;
}


