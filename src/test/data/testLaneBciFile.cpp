// data/testCycleBCIFile.cpp

#include <vector>
#include <string>
#include <boost/format.hpp>
#include <iostream>
#include <fstream>
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
// TESTING
//
//=========================================

TEST(TestLaneBci,LaneBciRead) {

    std::string filename = bgzfData::inputPathString + "L001/s_1.bci";
    std::ifstream in(filename.c_str(), std::ios_base::in | std::ios_base::binary);

    int counter = 0;
    uint32_t tileNumber;
    uint32_t clusterCount;

    while(in) {
        readBinary(in, tileNumber);
        readBinary(in, clusterCount);
        //if (in) {
        //    std::cerr << ++counter << " Tile: " << tileNumber << " Clusters: " << clusterCount << std::endl;
        //}
    }
}

