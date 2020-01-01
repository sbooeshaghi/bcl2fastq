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
#include "data/PositionsFile.hh"
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

// DEBUGGING TEST === 

//TEST(TestPositionsFile, readLocs) {
//
//    std::string filename = "/illumina/scratch/STE_DataAnalysis/test_input/Bcl2Fastq/23_TDI861_miniSeq_si_pe_4_p11/Data/Intensities/L001/s_1.locs";
//    size_t fileByteCount = 260171052;
//    size_t readResult = 80000;
//    //std::string filename = "/illumina/scratch/basespace-developers/dberard/data/genfna388/Data/Intensities/L001/s_1.locs";
//    //size_t fileByteCount = 129549748;
//    //size_t readResult = 80000;
//
//    common::RawDataBuffer inputData;
//    //io::SyncFile positionsFile(std::ios_base::in | std::ios_base::binary);
//    //positionsFile.readEntireFile(inputData);
//    io::UnprocessedFile positionsFile(filename, false);
//    positionsFile.readEntireFile(inputData);
//
//    EXPECT_EQ(fileByteCount, inputData.size());
//
//    size_t targetSize = fileByteCount / 8;
//    bool ignoreErrors = false;
//    bool skipHeader = true;
//    std::vector<data::PositionsFile::Record> targetBuffer;
//
//    std::streamsize byteCount = data::PositionsFileFactory::read(
//        inputData, ignoreErrors, skipHeader, targetBuffer, targetSize);
//
//    //EXPECT_EQ(8000, byteCount);
//    //EXPECT_EQ(1000, targetBuffer.size());
//   
//    for (data::PositionsFile::Record record : targetBuffer) {
//
//        std::cerr << "Rec: " << record.x_ << " , " << record.y_ << std::endl;
//    }
//
//}

