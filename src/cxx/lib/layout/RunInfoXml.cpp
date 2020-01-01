/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file RunInfoXml.cpp
 *
 * \brief Implementation of RunInfo.xml helper.
 *
 * \author Marek Balint
 */


#include <fstream>

#include <boost/property_tree/xml_parser.hpp>
#include <boost/optional.hpp>
#include <boost/format.hpp>

#include "common/Logger.hh"
#include "common/Exceptions.hh"
#include "io/Xml.hh"
#include "layout/RunInfoXml.hh"


namespace bcl2fastq {


namespace layout {


namespace detail {


static const std::string runInfoXmlFileName("RunInfo.xml");


} // namespace detail


RunInfoXml::RunInfoXml(boost::property_tree::ptree runInfoData)
: ptree_(runInfoData)
, reads_()
{
    this->initReadMetadata();
    this->initTileNumbers();
}

std::string RunInfoXml::getInstrument() const
{
    boost::optional<std::string> ret = ptree_.get_optional<std::string>("RunInfo.Run.Instrument");
    return ret ? ret.get() : "";
}

std::string RunInfoXml::getRunNumber() const
{
    boost::optional<std::string> ret = ptree_.get_optional<std::string>("RunInfo.Run.<xmlattr>.Number");
    return ret ? ret.get() : "";
}

std::string RunInfoXml::getRunId() const
{
    boost::optional<std::string> ret = ptree_.get_optional<std::string>("RunInfo.Run.<xmlattr>.Id");
    return ret ? ret.get() : "";
}

std::string RunInfoXml::getFlowcellId() const
{
    boost::optional<std::string> ret = ptree_.get_optional<std::string>("RunInfo.Run.Flowcell");
    return ret ? ret.get() : "";
}

common::LaneNumber RunInfoXml::getLanesCount() const
{
    boost::optional<common::LaneNumber> ret = ptree_.get_optional<common::LaneNumber>("RunInfo.Run.FlowcellLayout.<xmlattr>.LaneCount");
    if(!ret) {
        BOOST_THROW_EXCEPTION(common::InputDataError("Required attribute 'LaneCount' of element 'RunInfo.Run.FlowcellLayout' missing in RunInfo.xml file"));
    }
    return ret.get();
}

ReadMetadataContainer::const_iterator RunInfoXml::readMetadataBegin() const
{
    return reads_.begin();
}

ReadMetadataContainer::const_iterator RunInfoXml::readMetadataEnd() const
{
    return reads_.end();
}

RunInfoXml::TileNumbersContainer::const_iterator RunInfoXml::tileNumbersBegin( common::LaneNumber lane ) const
{
    BCL2FASTQ_ASSERT_MSG( 0 < lane && lane <= getLanesCount(), "Lane '" << lane << "' is out of range" );
    return tiles_[lane-1].begin();
}

RunInfoXml::TileNumbersContainer::const_iterator RunInfoXml::tileNumbersEnd( common::LaneNumber lane ) const
{
    BCL2FASTQ_ASSERT_MSG( 0 < lane && lane <= getLanesCount(), "Lane '" << lane << "' is out of range" );
    return tiles_[lane-1].end();
}

void RunInfoXml::initReadMetadata()
{
    boost::optional<boost::property_tree::ptree &> readsElement = ptree_.get_child_optional("RunInfo.Run.Reads");
    if (!readsElement) {
        BOOST_THROW_EXCEPTION(common::InputDataError("Required element 'RunInfo.Run.Reads' missing in RunInfo.xml file"));
    }

    // magic number 1: read number is 1-based
    common::ReadNumber dataReadsCounter = 1;
    common::ReadNumber indexReadsCounter = 1;
    // magic number 1: cycle number is 1-based
    common::CycleNumber cyclesCounter = 1;

    BOOST_FOREACH (const boost::property_tree::ptree::value_type &readElement, readsElement.get())
    {
        static const std::string readsElementName("Read");
        if (readElement.first != readsElementName)
        {
            continue;
        }

        boost::optional<const boost::property_tree::ptree &> indexElement = readElement.second.get_child_optional("Index");
        boost::optional<std::string> indexAttribute = readElement.second.get_optional<std::string>("<xmlattr>.IsIndexedRead");
        common::ReadNumber readNumber;
        bool indexReadFlag;
        if (indexElement || (indexAttribute && indexAttribute.get() == "Y"))
        {
            readNumber = indexReadsCounter++;
            indexReadFlag = true;
        }
        else
        {
            readNumber = dataReadsCounter++;
            indexReadFlag = false;
        }

        boost::optional<common::CycleNumber> firstCycle = readElement.second.get_optional<common::CycleNumber>("<xmlattr>.FirstCycle");
        boost::optional<common::CycleNumber> lastCycle = readElement.second.get_optional<common::CycleNumber>("<xmlattr>.LastCycle");
        if (!firstCycle)
        {
            firstCycle = cyclesCounter;
        }
        if (!lastCycle)
        {
            boost::optional<common::CycleNumber> numCycles = readElement.second.get_optional<common::CycleNumber>("<xmlattr>.NumCycles");
            if(!numCycles) {
                BOOST_THROW_EXCEPTION(common::InputDataError("Either 'FirstCycle' and 'LastCycle' or 'NumCycles' attribute in element 'RunInfo.Run.Reads.Read' is required in RunInfo.xml file"));
            }
            // magic number -1: lastCycle needs to reference last existing cycle as opposed to one-past-last cycle
            lastCycle = firstCycle.get() + numCycles.get() - 1;
        }
        // magic nubmer +1: lastCycle refers to last existing cycle in this read, next read needs to start with the next cycle
        cyclesCounter = lastCycle.get() + 1;

        reads_.push_back(ReadMetadata(ReadMetadata::Range(firstCycle.get(),
                                                          lastCycle.get()),
                                                          readNumber,
                                                          indexReadFlag ? common::ReadType::INDEX : common::ReadType::DATA));
    }
}

void RunInfoXml::initTileNumbers()
{
    boost::optional<boost::property_tree::ptree &> flowcellLayout = ptree_.get_child_optional("RunInfo.Run.FlowcellLayout");
    if (!flowcellLayout) {
        BOOST_THROW_EXCEPTION(common::InputDataError("Required element 'RunInfo.Run.FlowcellLayout' missing in RunInfo.xml file"));
    }
    const unsigned laneCount = flowcellLayout.get().get<unsigned>("<xmlattr>.LaneCount");
    const unsigned surfaceCount = flowcellLayout.get().get<unsigned>("<xmlattr>.SurfaceCount");
    const unsigned swathCount = flowcellLayout.get().get<unsigned>("<xmlattr>.SwathCount");
    const unsigned tileCount = flowcellLayout.get().get<unsigned>("<xmlattr>.TileCount");
    const unsigned sectionPerLane = flowcellLayout.get().get<unsigned>("<xmlattr>.SectionPerLane", 0);
    const unsigned lanePerSection = flowcellLayout.get().get<unsigned>("<xmlattr>.LanePerSection", 0);

    for (unsigned lane = 0; lane < laneCount; ++lane)
    {
        TileNumbersContainer tmp;
        for (unsigned surface = 1; surface <= surfaceCount; ++surface)
        {
            for (unsigned swath = 1; swath <= swathCount; ++swath)
            {
                for (unsigned tile = 1; tile <= tileCount; ++tile)
                {
                    if (sectionPerLane > 1 || lanePerSection > 1)
                    {
                        BCL2FASTQ_ASSERT_MSG(lanePerSection, "Unexpected Zero lanePerSection");
                        const unsigned firstSection = 1 + sectionPerLane * (lane/lanePerSection);
                        for (unsigned section = firstSection; section <= firstSection + sectionPerLane; ++section)
                        {
                            tmp.push_back( common::TileNumber( tile + section * 100 + swath * 1000 + surface * 10000 ));
                        }
                    }
                    else
                    {
                        tmp.push_back( common::TileNumber( tile + swath * 100 + surface * 1000 ));
                    }
                }
            }
        }
        tiles_.push_back( tmp );
    }
}

RunInfoXml createRunInfoXml(boost::filesystem::path runfolderDir)
{
    boost::filesystem::path file(runfolderDir / boost::filesystem::path(detail::runInfoXmlFileName));
    BCL2FASTQ_LOG(common::LogLevel::INFO) << "RunInfo.xml: '" << file << "'" << std::endl;
    return RunInfoXml(io::parseXmlFile(file));
}


} // namespace layout


} // namespace bcl2fastq


