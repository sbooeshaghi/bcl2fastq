
#include <vector>
#include <string>
#include <boost/format.hpp>
#include <iostream>
#include <ios>

#include "common/Debug.hh"
//#include "common/Logger.hh"
#include "common/SystemCompatibility.hh"
#include "io/SyncFile.hh"

#include "data/InteropFile.hh"
#include "gtest/gtest.h"

using namespace bcl2fastq;

namespace {

    template<class T>
    bool readBinary(std::ifstream& in, T& buffer)
    {
        in.read(reinterpret_cast<char *>(&buffer), sizeof(T));
        return (in) ? true : false;
    }

} // end unnamed namespace

TEST(InteropFactoryTest, CreateInteropRecordV1) {

    EXPECT_EQ((int)1, (int)(data::InteropFactory::V1));

    data::IInteropRecordWriterPtr interopRecordWriterV1 = data::InteropFactory::createInteropRecordWriter(data::InteropFactory::V1);

    data::InteropFile fileV1(data::InteropFactory::V1);
    EXPECT_EQ(data::InteropFactory::V1, fileV1.getVersion());
}

TEST(InteropFactoryTest, CreateInteropRecordV2) {

    EXPECT_EQ((int)2, (int)(data::InteropFactory::V2));

    data::IInteropRecordWriterPtr interopRecordWriterV2 = data::InteropFactory::createInteropRecordWriter(data::InteropFactory::V2);

    data::InteropFile fileV2(data::InteropFactory::V2);
    EXPECT_EQ(data::InteropFactory::V2, fileV2.getVersion());
}

TEST(InteropFileTest, CreateInteropRecord) {

    data::InteropRecord record1(1,1103,2,0x0a23,"AACCGGTT-ACGTACGT","Sample ID 423", "Project NameIt");
    EXPECT_EQ(   1, record1.laneNumber);
    EXPECT_EQ(1103, record1.tileNumber);
    EXPECT_EQ(   2, record1.readNumber);
    EXPECT_EQ(2595, record1.clustersCount);
    EXPECT_EQ("AACCGGTT-ACGTACGT", record1.indexName);
    EXPECT_EQ("Sample ID 423", record1.sampleId);
    EXPECT_EQ("Project NameIt", record1.projectName);

//
    bcl2fastq::data::InteropRecord fiveDigitRecord(1,11102,2,0x0a23,"AACCGGTT-ACGTACGT","Sample ID 423", "Project NameIt");
    EXPECT_EQ(    1, fiveDigitRecord.laneNumber);
    EXPECT_EQ(11102, fiveDigitRecord.tileNumber);
    EXPECT_EQ(    2, fiveDigitRecord.readNumber);
    EXPECT_EQ( 2595, fiveDigitRecord.clustersCount);
    EXPECT_EQ("AACCGGTT-ACGTACGT", fiveDigitRecord.indexName);
    EXPECT_EQ("Sample ID 423", fiveDigitRecord.sampleId);
    EXPECT_EQ("Project NameIt", fiveDigitRecord.projectName);

}

TEST(InteropFileTest, WriterInteropRecordV1) {

    data::InteropRecord record1(1,1103,2,0x0a23,"AACCGGTT-ACGTACGT","Sample ID 423", "Project NameIt");
    data::IInteropRecordWriterPtr interopRecordWriterV1 = data::InteropFactory::createInteropRecordWriter(data::InteropFactory::V1);

    bcl2fastq::data::CharVector byteVector;
    size_t recordSize = interopRecordWriterV1->writeInteropRecord(byteVector, record1);
    size_t expected = 16 + record1.indexName.size()
            + record1.sampleId.size()
            + record1.projectName.size();
    EXPECT_EQ(expected, byteVector.size());
    EXPECT_EQ(expected, recordSize);

//
    bcl2fastq::data::InteropRecord fiveDigitRecord(1,11102,2,0x0a23,"AACCGGTT-ACGTACGT","Sample ID 423", "Project NameIt");

    recordSize = interopRecordWriterV1->writeInteropRecord(byteVector, fiveDigitRecord);
    expected = 16 + fiveDigitRecord.indexName.size()
            + fiveDigitRecord.sampleId.size()
            + fiveDigitRecord.projectName.size();
    EXPECT_EQ(expected, byteVector.size());
    EXPECT_EQ(expected, recordSize);
}

TEST(InteropFileTest, WriterInteropRecordV2) {

    data::InteropRecord record1(1,1103,2,0x0a23,"AACCGGTT-ACGTACGT","Sample ID 423", "Project NameIt");
    data::IInteropRecordWriterPtr interopRecordWriterV2 = data::InteropFactory::createInteropRecordWriter(data::InteropFactory::V2);

    bcl2fastq::data::CharVector byteVector;
    size_t recordSize = interopRecordWriterV2->writeInteropRecord(byteVector, record1);
    size_t expected = 22 + record1.indexName.size()
            + record1.sampleId.size()
            + record1.projectName.size();
    EXPECT_EQ(expected, byteVector.size());
    EXPECT_EQ(expected, recordSize);

//
    bcl2fastq::data::InteropRecord fiveDigitRecord(1,11102,2,0x0a23,"AACCGGTT-ACGTACGT","Sample ID 423", "Project NameIt");

    recordSize = interopRecordWriterV2->writeInteropRecord(byteVector, fiveDigitRecord);
    expected = 22 + fiveDigitRecord.indexName.size()
            + fiveDigitRecord.sampleId.size()
            + fiveDigitRecord.projectName.size();
    EXPECT_EQ(expected, byteVector.size());
    EXPECT_EQ(expected, recordSize);
}

TEST(InteropFileTest, WriteInteropRecordToFileV1) {

    const boost::filesystem::path path("IndexMetricsOut.bin");
    {
        bcl2fastq::data::InteropFile interopFile(data::InteropFactory::V1);
        interopFile.open(path);

        bcl2fastq::data::InteropRecord r1(1,1103,2,0x0a23,"AACCGGTT-ACGTACGT","Sample ID 423", "Project NameIt");
        bcl2fastq::data::InteropRecord r2(1,1103,2,0x0d23,"AACCGGTT-ACGTACGT","Sample ID 523", "Project Big Wigg");
        interopFile.writeInteropRecord(r1);
        interopFile.writeInteropRecord(r2);
    }

    std::ifstream in("IndexMetricsOut.bin", std::ios::ate | std::ios_base::binary);

    std::streamsize fileSize = in.tellg();
    EXPECT_EQ(123, fileSize);

    in.seekg(0);
    char version;
    readBinary(in, version);
    EXPECT_EQ(1, (int)version);
}

TEST(InteropFileTest, WriteInteropRecordToFileV2) {

    const boost::filesystem::path path("IndexMetricsOutV2.bin");
    {
        bcl2fastq::data::InteropFile interopFile(data::InteropFactory::V2);
        interopFile.open(path);

        bcl2fastq::data::InteropRecord r1(1,1103,2,0x0a23,"AACCGGTT-ACGTACGT","Sample ID 423", "Project NameIt");
        bcl2fastq::data::InteropRecord r2(1,1103,2,0x0d23,"AACCGGTT-ACGTACGT","Sample ID 523", "Project Big Wigg");
        interopFile.writeInteropRecord(r1);
        interopFile.writeInteropRecord(r2);
    }

    std::ifstream in("IndexMetricsOutV2.bin", std::ios::ate | std::ios_base::binary);

    std::streamsize fileSize = in.tellg();
    EXPECT_EQ(135, fileSize);

    in.seekg(0);
    char version;
    readBinary(in, version);
    EXPECT_EQ((char)(2), version);
}

//TEST(InteropFileTest, ViewInterop) {
//
//    std::string filename = "/illumina/scratch/Voyager/170417_VP2-04_0258_A05P4DNXVY/InterOp/IndexMetricsOut.bin";
//    std::ifstream in(filename.c_str(), std::ios_base::in | std::ios_base::binary);
//
//    char version;
//    char buffer[256];
//
//    // Version
//    in.read(&version, sizeof(char));
//
//    for (int i = 0; i < 100; ++i) {
//        uint16_t  laneNumber;    readBinary(in, laneNumber);    std::cerr << " Lane: " << laneNumber;
//        uint16_t  tileNumber;    readBinary(in, tileNumber);    std::cerr << " Tile: " << tileNumber;
//        uint16_t  readNumber;    readBinary(in, readNumber);    std::cerr << " Read: " << readNumber;
//
//        uint16_t  length;        readBinary(in, length);
//        in.read(reinterpret_cast<char *>(&buffer), length);
//        std::string barcode(buffer); barcode.resize(length);
//        std::cerr << " Index: " << barcode;
//
//        uint32_t  clustersCnt;   readBinary(in, clustersCnt);   std::cerr << " Clust:" << clustersCnt;
//
//        readBinary(in, length);
//        in.read(reinterpret_cast<char *>(&buffer), length);
//        std::string sampleId(buffer); sampleId.resize(length);
//        std::cerr << " Sample: " << sampleId;
//
//        readBinary(in, length);
//        in.read(reinterpret_cast<char *>(&buffer), length);
//        std::string project(buffer); project.resize(length);
//        std::cerr << " project: " << project << std::endl;
//    }
//}

