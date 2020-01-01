// testReadInfo

#include "common/Types.hh"
#include "layout/ReadInfo.hh"
#include "gtest/gtest.h"

using namespace bcl2fastq;

namespace {

    const std::string runInfoXmlData(
        "<?xml version=\"1.0\"?>\n"
        "<RunInfo xmlns:xsd=\"http://www.w3.org/2001/XMLSchema\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" Version=\"2\">\n"
        "  <Run Id=\"TestRun\" Number=\"42\">\n"
        "    <Flowcell>TestFlowcell</Flowcell>\n"
        "    <Instrument>CSSIM</Instrument>\n"
        "    <Date>8/20/2013 3:29:17 PM</Date>\n"
        "    <Reads>\n"
        "      <Read Number=\"1\" NumCycles=\"151\" IsIndexedRead=\"N\" />\n"
        "      <Read Number=\"2\" NumCycles=\"8\" IsIndexedRead=\"Y\" />\n"
        "      <Read Number=\"3\" NumCycles=\"151\" IsIndexedRead=\"N\" />\n"
        "    </Reads>\n"
        "    <FlowcellLayout LaneCount=\"4\" SurfaceCount=\"2\" SwathCount=\"3\" TileCount=\"11\" />\n"
        "  </Run>\n"
        "</RunInfo>\n"
    );
 

} // end unnamed namespace

TEST(TestReadInfo, ReadInfoCtor) {

    common::ReadNumber readNumber  = 1;
    common::ReadType   readType    = common::DATA; // INDEX, UMI
    common::CycleRange cyclesToUse = std::make_pair(1,151);
    common::CycleRange unmasked    = std::make_pair(1,151);
    common::CycleRange umi         = std::make_pair(0,0);

    layout::ReadInfo readInfo(readNumber, readType, cyclesToUse, unmasked, umi);

    EXPECT_EQ(1, readInfo.getNumber());
    EXPECT_EQ(common::DATA, readInfo.getReadType());

    EXPECT_EQ(151, readInfo.cycleInfos().size());
    EXPECT_EQ(151, readInfo.unmaskedCycleInfos().size());
    EXPECT_EQ(0, readInfo.getUmiCycles().first);
}

