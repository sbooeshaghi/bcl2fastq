// data/testCycleBCIFile.cpp

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
#include "data/CycleBCIFile.hh"
#include "gtest/gtest.h"
#include "../bgzfData.hh"

using namespace bcl2fastq;

namespace {

    template<class T>
    bool readBinary(std::ifstream& in, T& buffer)
    {
        in.read(reinterpret_cast<char *>(&buffer), sizeof(T));
        return (in) ? true : false;
    }

} // end unnamed namespace


//=========================================
//
// TESTING A
//
//=========================================

TEST(testCycleBciFile, Read) {

    std::string filename = bgzfData::inputPathString + "L001/0001.bcl.bgzf.bci";
    std::ifstream in(filename.c_str(), std::ios_base::in | std::ios_base::binary);

    int counter = 0;
    uint32_t version;
    uint32_t tileCount;

    readBinary(in, version);
    readBinary(in, tileCount);
    //std::cerr << " Version: " << version << " # tiles: " << tileCount << std::endl;

    data::CycleBCIFile::Record record;
    while(in) {
        readBinary(in, record);
        //if (in) {
        //    std::cerr << ++counter << " uncOffset: " << record.uncompressedOffset << " zOffset: " << record.compressedOffset << std::endl;
        //}
    }
}

TEST(TestCycleBCIFile, getAndVerify) {

    const boost::filesystem::path inputPath(bgzfData::inputPathString);
    common::LaneNumber  laneNumber = 1;
    common::CycleNumber cycleNumber= 3;
    const bool ignoreMissingDataFile = true;
    boost::filesystem::path filename;

    bool result = data::CycleBCIFile::getAndVerifyFileName(
        inputPath, laneNumber, cycleNumber, ignoreMissingDataFile, filename);
    EXPECT_TRUE(result);

    std::string fullName = bgzfData::inputPathString + "L001/0003.bcl.bgzf.bci";
    EXPECT_EQ(fullName, filename.string());

//    data::CycleBCIFile cycleBciFile(ignoreMissingDataFile);
}

TEST(TestCycleBCIFile, openFile) {

    const boost::filesystem::path inputPath(bgzfData::inputPathString);
    common::LaneNumber  laneNumber = 1;
    common::CycleNumber cycleNumber= 3;
    const bool ignoreMissingDataFile = true;

    data::CycleBCIFile cycleBciFile(ignoreMissingDataFile);

    bool opened = cycleBciFile.openFile(inputPath, laneNumber, cycleNumber);
    EXPECT_TRUE(opened);

    std::string fullName = bgzfData::inputPathString + "L001/0003.bcl.bgzf.bci";
    std::string filename = cycleBciFile.getPath().string();
    EXPECT_EQ(fullName, filename);

    EXPECT_EQ(216, cycleBciFile.getTilesCount());
}

TEST(TestCycleBCIFile, read) {

    const boost::filesystem::path inputPath(bgzfData::inputPathString);
    common::LaneNumber  laneNumber = 1;
    common::CycleNumber cycleNumber= 1;
    const bool ignoreMissingDataFile = true;

    data::CycleBCIFile cycleBciFile(ignoreMissingDataFile);
    bool opened = cycleBciFile.openFile(inputPath, laneNumber, cycleNumber);
    EXPECT_TRUE(opened);
}

TEST(TestCycleBCIFile, FileSize) {

    const boost::filesystem::path inputPath(bgzfData::inputPathString);
    common::LaneNumber  laneNumber = 1;
    common::CycleNumber cycleNumber= 1;
    const bool ignoreMissingDataFile = true;

    data::CycleBCIFile cycleBciFile(ignoreMissingDataFile);
    bool opened = cycleBciFile.openFile(inputPath, laneNumber, cycleNumber);
    EXPECT_TRUE(opened);
}

TEST(TestCycleBCIFile, GetRecord) {

    const boost::filesystem::path inputPath(bgzfData::inputPathString);
    common::LaneNumber  laneNumber = 1;
    common::CycleNumber cycleNumber= 1;
    const bool ignoreMissingDataFile = true;

    data::CycleBCIFile cycleBciFile(ignoreMissingDataFile);
    bool opened = cycleBciFile.openFile(inputPath, laneNumber, cycleNumber);

    for (int tileIndex = 0; tileIndex < 216; ++tileIndex) {
        data::CycleBCIFile::Record record = cycleBciFile.getRecord(tileIndex);
        //std::cerr << record.compressedOffset << std::endl;
    } 
    EXPECT_TRUE(opened);
}

TEST(TestCycleBCIFile, failedRead) {

    const boost::filesystem::path inputPath(bgzfData::inputPathString);
    common::LaneNumber  laneNumber = 1;
    common::CycleNumber cycleNumber= 3; // 1;
    const bool ignoreMissingDataFile = true;

    data::CycleBCIFile cycleBciFile(ignoreMissingDataFile);
    bool opened = cycleBciFile.openFile(inputPath, laneNumber, cycleNumber);
    EXPECT_TRUE(opened);
}


