// data/testCycleBCIFile.cpp

#include "gtest/gtest.h"

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
#include "data/AggregatedBclFileReader.hh"
#include "data/RawBclBuffer.hh"
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
// TESTING
//
//=========================================

TEST(TestBgzfFile, readBgzfFileWithoutBcl2fastq) {

    std::string filename = bgzfData::inputPathString + "L001/0001.bcl.bgzf";
    std::ifstream in(filename.c_str(), std::ios_base::in | std::ios_base::binary);

    uint32_t pos = 0;

    while (in) {

        uint8_t  id1  ;  readBinary(in, id1  );
        uint8_t  id2  ;  readBinary(in, id2  );
        uint8_t  cm   ;  readBinary(in, cm   );
        uint8_t  flg  ;  readBinary(in, flg  );
        uint32_t mtime;  readBinary(in, mtime);
        uint8_t  xflg ;  readBinary(in, xflg );
        uint8_t  os   ;  readBinary(in, os   );
        uint16_t xlen ;  readBinary(in, xlen );
        uint8_t  si1  ;  readBinary(in, si1  );
        uint8_t  si2  ;  readBinary(in, si2  );
        uint16_t slen ;  readBinary(in, slen );
        uint16_t bsize;  readBinary(in, bsize);
        uint8_t  cdata;  readBinary(in, cdata);
        uint32_t crc32;  readBinary(in, crc32);
        uint32_t isize;  readBinary(in, isize);

        pos += bsize + 1;

        in.seekg((std::streampos)(pos));
    }
}

TEST(TestBgzfFile, readBgzfWithBcl2Fastq) {

    // Setup 
    const boost::filesystem::path inputPath(bgzfData::inputPathString);
    common::LaneNumber  laneNumber = 1;
    common::CycleNumber cycleNumber= 1;
    const bool ignoreMissingDataFile = true;
    int tileCount = 144;
    int cycleCount = 52;

    // DEFINE layout::LaneInfo stuff
    layout::LaneInfosContainer lanes;
    layout::LaneInfo laneInfo(laneNumber);
    lanes.push_back(laneInfo);

    //TODO: Setup reads?

    // Define tiles
    int tileIndex = 0;
    int skippedTiles = 0;
    for (int i = 0; i < tileCount; ++i) {
        laneInfo.addTile(layout::TileInfo(bgzfData::testTiles[tileIndex], tileIndex, skippedTiles));
        ++tileIndex;
    }
    std::cerr << "done adding tiles\n" << std::endl;

    // Cycle BCI file
    std::shared_ptr<data::CycleBCIFile> cycleBciFile(std::make_shared<data::CycleBCIFile>(ignoreMissingDataFile));
    bool openedBci = cycleBciFile->openFile(inputPath, laneNumber, cycleNumber);
    EXPECT_TRUE(openedBci);

    // 
    std::shared_ptr<io::SyncFile> bgzfFile(std::make_shared<io::SyncFile>(std::ios_base::in | std::ios_base::binary));

    // Read bgzf
    data::AggregatedBclFileReader bgzf(
        inputPath, laneInfo, cycleNumber, 0, ignoreMissingDataFile, bgzfFile, cycleBciFile);

    data::RawBclBufferGroup bclData(tileCount,cycleCount); // (144,52);

    data::RawBclBufferVec::iterator it;
    const layout::LaneInfo::TileInfosContainer& tileInfos = laneInfo.getTileInfos();

    std::cerr << "Associate tiles " << std::endl;
    layout::LaneInfo::TileInfosContainer::const_iterator tileInfo = tileInfos.begin();
    for ( it = bclData.begin(); it < bclData.end(); ++it) {

        it->tileInfo_ = tileInfo;
        ++tileInfo;
    }

    std::cerr << "Set buffers " << std::endl;
    for (auto& tileBuffer : bclData) {

        tileBuffer.cycleData_.resize(cycleCount);
        tileBuffer.gzipDecompressors_.resize(tileBuffer.cycleData_.size());
    }

    // todo: define RawBclBuffer content
    // todo: define cycles 
    // todo: associate tileInfoIter

    common::LogLevel::value_type oldLevel = common::BCL2FASTQ_MIN_LOG_LEVEL;
    common::BCL2FASTQ_MIN_LOG_LEVEL = common::LogLevel::TRACE;

    bgzf.read(bclData);

    common::BCL2FASTQ_MIN_LOG_LEVEL = oldLevel;
}

