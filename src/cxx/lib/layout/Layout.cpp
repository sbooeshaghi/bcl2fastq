/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file Layout.cpp
 *
 * \brief Implementation of data layout.
 *
 * \author Marek Balint
 * \author Mauricio Varea
 */


#include <utility>
#include <numeric>
#include <sstream>

#include <boost/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>

#include "common/Logger.hh"
#include "common/Types.hh"
#include "layout/Layout.hh"
#include "layout/UseBasesMask.hh"
#include "layout/RunInfoXml.hh"
#include "config/SampleSheetCsv.hh"
#include "layout/ConfigXml.hh"
#include "layout/BCIndex.hh"
#include "common/DirectoryValidator.hh"
#include "layout/BarcodeCollisionDetector.hh"
#include "layout/FileExistenceVerifier.hh"
#include "data/CbclFile.hh"
#include "io/SyncFile.hh"


namespace bcl2fastq {


namespace layout {

typedef std::map<std::string, common::SampleNumber> SampleIdNumberMap;

const FlowcellInfo & Layout::getFlowcellInfo() const
{
    return flowcellInfo_;
}

LaneInfosContainer::const_iterator Layout::laneInfosBegin() const
{
    return laneInfos_.begin();
}

LaneInfosContainer::const_iterator Layout::laneInfosEnd() const
{
    return laneInfos_.end();
}

bool includeTile(const std::string&               flowcell,
                 common::LaneNumber               laneNumber,
                 common::TileNumber               tileNumber,
                 const std::vector<boost::regex>& tileRegexps,
                 const config::SampleSheetCsv&    sampleSheet)
{
    if (tileRegexps.empty())
    {
        return !sampleSheet.isTileExcluded(flowcell, laneNumber, tileNumber);
    }
    else
    {
        std::string tileString((boost::format("s_%d_%04d") % laneNumber % tileNumber).str());
        BOOST_FOREACH(const boost::regex &re, tileRegexps)
        {
            if (boost::regex_search(tileString, re))
            {
                return true;
            }
        }
    }
    return false;
}

bool checkExistence(const boost::filesystem::path& laneDir,
                    const std::string& fileExt,
                    std::vector<boost::filesystem::path>& files,
                    bool findOne = true)
{
    BCL2FASTQ_LOG(common::LogLevel::DEBUG) << "Checking if there are any " << fileExt << " files in: " << laneDir << std::endl;

    boost::regex pattern(".*\\." + fileExt + "$", boost::regex_constants::mod_x);
    if( !boost::filesystem::is_directory( laneDir ))  return false;
    for( boost::filesystem::directory_iterator it(laneDir), end;  it != end;  ++it )
    {
        boost::smatch match;
        if( boost::filesystem::is_regular_file( it->status() )
         && boost::regex_match( it->path().leaf().string(), match, pattern ))
        {
            if (findOne)
            {
                BCL2FASTQ_LOG(common::LogLevel::DEBUG) << "\t... found " << it->path() << std::endl;
                return true;
            }

            files.push_back(it->path());
        }
    }

    if (!files.empty())
    {
        return true;
    }

    BCL2FASTQ_LOG(common::LogLevel::DEBUG) << "\t... no " << fileExt << " files found!" << std::endl;
    return false;
}

bool checkBciExistence(const boost::filesystem::path& laneDir)
{
    std::vector<boost::filesystem::path> unused;
    return checkExistence(laneDir, "bci", unused);
}

bool checkCbclExistence(const boost::filesystem::path& laneDir)
{
    std::vector<boost::filesystem::path> unused;
    return checkExistence(laneDir / "C1.1", "cbcl", unused);
}

enum AdapterType
{
    TRIM,
    MASK
};

unsigned int addAdaptersToRead( config::SampleSheetCsv::AdaptersContainer::const_iterator adaptersBegin,
                                config::SampleSheetCsv::AdaptersContainer::const_iterator adaptersEnd,
                                ReadInfo &readInfo,
                                AdapterType adapterType)
{
    BCL2FASTQ_ASSERT_MSG(adapterType == TRIM || adapterType == MASK, "Unrecognized adapter type");

    unsigned int adapterCount = 0;
    bool haveNonZeroReadAdapter = false;
    BOOST_FOREACH (const config::SampleSheetCsv::AdaptersContainer::value_type &adapter, std::make_pair(adaptersBegin, adaptersEnd))
    {
        if (adapter.first != 0)
        {
        	haveNonZeroReadAdapter = true;
        }
    }
    BOOST_FOREACH (const config::SampleSheetCsv::AdaptersContainer::value_type &adapter, std::make_pair(adaptersBegin, adaptersEnd))
    {
        bool acceptAdapter = false;
        if(!haveNonZeroReadAdapter)
        {
            acceptAdapter = true;
        }
        else if (adapter.first == readInfo.getNumber())
        {
            acceptAdapter = true;
        }
        else if((adapter.first == 0) && (readInfo.getNumber() == 1))
        {
            acceptAdapter = true;
        }

        if (acceptAdapter)
        {
            std::string adapterTypeStr;
            if (adapterType == MASK)
            {
                readInfo.addMaskAdapter(adapter.second);
                adapterTypeStr = "mask";
            }
            else if (adapterType == TRIM)
            {
                readInfo.addTrimAdapter(adapter.second);
                adapterTypeStr = "trim";
            }

            ++adapterCount;
            BCL2FASTQ_LOG(common::LogLevel::INFO) << "  Adapter: " << adapter.second << (boost::format(" (%s adapter)") % adapterTypeStr).str() << std::endl;
        }
    }
    return adapterCount;
}

void validateBarcodes(const LaneInfo& laneInfo,
                      BarcodeCollisionDetector& barcodeCollisionDetector)
{
    if (laneInfo.getSampleInfos().size() == 1)
    {
        BCL2FASTQ_ASSERT_MSG(laneInfo.getSampleInfos().front().getBarcodes().empty(),
                             "Internal programming error. Barcode exists for lane with only 1 sample");

        BCL2FASTQ_LOG(common::LogLevel::INFO) << "  Sample: " << laneInfo.getSampleInfos().front() << std::endl;

        return;
    }

    for (const auto& sampleInfo : laneInfo.getSampleInfos())
    {
        BCL2FASTQ_LOG(common::LogLevel::INFO) << "  Sample: " << sampleInfo << std::endl;

        if (sampleInfo.getNumber() == 0)
        {
            // No barcodes for the default sample
            continue;
        }

        for (const auto& barcode : sampleInfo.getBarcodes())
        {
            BCL2FASTQ_LOG(common::LogLevel::INFO) << "    Barcode: '" << barcode << "'" << std::endl;
        }

        barcodeCollisionDetector.validateBarcode(sampleInfo.getBarcodes());
    }
}

void Layout::detectFlowcellInfo(const boost::filesystem::path& intensitiesDir,
                                const RunInfoXml&              runInfoXml,
                                bool                           hasBci,
                                bool                           hasCbcl)
{
    flowcellInfo_.setInstrument(runInfoXml.getInstrument());
    flowcellInfo_.setRunNumber(runInfoXml.getRunNumber());
    flowcellInfo_.setFlowcellId(runInfoXml.getFlowcellId());
    flowcellInfo_.setRunId(runInfoXml.getRunId());

    if( boost::filesystem::exists( intensitiesDir / "s.locs" ) )
    {
        flowcellInfo_.setPatternedFlowcell();
        BCL2FASTQ_LOG(common::LogLevel::INFO) << "Patterned flowcell detected" << std::endl;
        BCL2FASTQ_LOG(common::LogLevel::DEBUG) << "This may be a 'HiSeq X Ten' run" << std::endl;
    } else {
        flowcellInfo_.resetPatternedFlowcell();
        BCL2FASTQ_LOG(common::LogLevel::DEBUG) << "NO patterned flowcell detected" << std::endl;
    }

    if( hasBci )
    {
        flowcellInfo_.setAggregateTilesFlag(common::TileAggregationMode::AGGREGATED);
        BCL2FASTQ_LOG(common::LogLevel::DEBUG) << "This may be a 'NextSeq 500' run" << std::endl;
    }
    else if ( hasCbcl )
    {
        flowcellInfo_.setAggregateTilesFlag(common::TileAggregationMode::CBCL);
    }
    else
    {
        flowcellInfo_.setAggregateTilesFlag(common::TileAggregationMode::NON_AGGREGATED);
    }
}

void parseUseBasesMasks(const std::vector<std::string>&            useBasesMaskStrs,
                        std::map<common::LaneNumber, std::string>& useBasesMasksForLane,
                        std::string&                               defaultUseBasesMask)
{
    BOOST_FOREACH(const std::string& useBasesMask, useBasesMaskStrs)
    {
        std::vector<std::string> useBasesMaskAndLane;
        boost::split(useBasesMaskAndLane, useBasesMask, boost::is_any_of(":"));

        if (useBasesMaskAndLane.size() == 1)
        {
            BCL2FASTQ_ASSERT_MSG(defaultUseBasesMask.empty(), "Multiple use-bases-mask paramters entered without specifying lane number.");

            defaultUseBasesMask = useBasesMask;
        }
        else if (useBasesMaskAndLane.size() == 2)
        {
            boost::algorithm::trim<std::string>(useBasesMaskAndLane[0]);
            boost::algorithm::trim<std::string>(useBasesMaskAndLane[1]);

            useBasesMasksForLane.insert(std::make_pair(boost::lexical_cast<unsigned int>(useBasesMaskAndLane[0]),
                                        useBasesMaskAndLane[1]));
        }
        else
        {
            BCL2FASTQ_ASSERT_MSG(false, "Invalid use-bases-mask: " + useBasesMask);
        }
    }
}

void setReadMetadata(common::CycleNumber startCycle,
                     common::CycleNumber endCycle,
                     layout::ReadMetadata& readMetadata)
{
    if (endCycle >= 1)
    {
        readMetadata.lastCycle_ = readMetadata.firstCycle_ + (endCycle-1);
    }
    if (startCycle >= 1)
    {
        readMetadata.firstCycle_ += (startCycle-1);
    }
}

void setReadMetadata(common::CycleNumber   startCycle,
                     common::CycleNumber   endCycle,
                     common::CycleNumber   umiStartFromCycle,
                     common::CycleNumber   umiLength,
                     bool                  trimUmi,
                     layout::ReadMetadata& readMetadata)
{
    common::CycleNumber cycleRange = readMetadata.lastCycle_ - readMetadata.firstCycle_ + 1;
    if (startCycle == 0)
    {
        startCycle = 1;
    }
    if (endCycle == 0)
    {
        endCycle = cycleRange;
    }

    if (startCycle < 1 ||
        startCycle > endCycle ||
        startCycle > cycleRange ||
        endCycle > cycleRange )
    {
        BOOST_THROW_EXCEPTION(common::InputDataError("Specified start/end cycles in read are outside of the range found in RunInfo.xml"));
    }


    if (umiStartFromCycle > startCycle &&
        endCycle + 1 > umiStartFromCycle + umiLength)
    {
        BOOST_THROW_EXCEPTION(common::InputDataError(
                "A UMI in the center of a read is not supported."));
    }

    if (umiStartFromCycle > readMetadata.lastCycle_)
    {
        BOOST_THROW_EXCEPTION(common::InputDataError(
                "UMI start cycle is after the last cycle in the read."));
    }

    if (trimUmi)
    {
        common::CycleNumber umiEndCycle = umiStartFromCycle + umiLength - 1;
        if (umiEndCycle >= startCycle &&
            umiStartFromCycle <= endCycle)
        {
            if (umiStartFromCycle <= startCycle)
            {
                startCycle = umiEndCycle + 1;
            }
            else
            {
                endCycle = umiStartFromCycle - 1;
            }
        }
    }

    setReadMetadata(startCycle,
                    endCycle,
                    readMetadata);
}

void setIndexReadMetadata(const std::vector<size_t>& barcodeLengths,
                          layout::ReadMetadata&      readMetadata,
                          size_t                     indexNumber)
{
    if (barcodeLengths.size() > indexNumber)
    {
        if (barcodeLengths[indexNumber] > 1 + readMetadata.lastCycle_ - readMetadata.firstCycle_)
        {
            BOOST_THROW_EXCEPTION(common::InputDataError((boost::format(
                "Barcodes in sample sheet are longer than the index length found in RunInfo.xml.")).str()));
        }

        if (barcodeLengths[indexNumber] == 0)
        {
            // Mask the read
            readMetadata.firstCycle_ = 1;
            readMetadata.lastCycle_ = 0;
        }
        else
        {
            setReadMetadata(0,
                            barcodeLengths[indexNumber],
                            readMetadata);
        }
    }
}

void validateIndexLengths(const layout::ReadMetadataContainer& readMetadataContainer,
                          const std::vector<size_t>&           barcodeLengths,
                          layout::LaneInfo&                    laneInfo)
{
    std::vector<std::pair<size_t,size_t>> indexNumberLengthPairs;
    for (const auto& readMetadata : readMetadataContainer)
    {
        if (readMetadata.readType_ == common::ReadType::INDEX)
        {
            indexNumberLengthPairs.push_back(std::make_pair(readMetadata.lastCycle_ - readMetadata.firstCycle_ + 1, readMetadata.readNumber_));
        }
    }

    for (const std::pair<size_t,size_t> indexNumberLengthPair : indexNumberLengthPairs)
    {
        if (indexNumberLengthPair.second > 2)
        {
            BOOST_THROW_EXCEPTION(common::InputDataError("Greater than 2 index reads encountered in use-bases-mask. This is not allowed."));
        }

        if ((indexNumberLengthPairs.size() == 2 && barcodeLengths.at(indexNumberLengthPair.second - 1) != indexNumberLengthPair.first) ||
            (barcodeLengths[0] != indexNumberLengthPair.first && barcodeLengths[1] != indexNumberLengthPair.first))
        {
            BOOST_THROW_EXCEPTION(common::InputDataError("Barcode lengths in the sample sheet do not match those in --use-bases-mask"));
        }
    }

    if (indexNumberLengthPairs.empty())
    {
        // There are no index sequences
        if (barcodeLengths[0] != 0 && barcodeLengths[1] != 0)
        {
            laneInfo.maskBarcode(2);
            laneInfo.maskBarcode(1);
        }
        else if (barcodeLengths[0] != 0 || barcodeLengths[1] != 0)
        {
            laneInfo.maskBarcode(1);
        }

        laneInfo.removeDefaultSample();
        if (laneInfo.getSampleInfos().size() > 1)
        {
            BOOST_THROW_EXCEPTION(common::InputDataError((boost::format(
                "Multiple samples with all barcodes masked for lane: %d") % laneInfo.getNumber()).str()));
        }
    }
    else if ((indexNumberLengthPairs.size() == 1 && barcodeLengths.size() == 2) &&
             (barcodeLengths[0] != 0 && barcodeLengths[1] != 0))
    {
        laneInfo.maskBarcode(indexNumberLengthPairs[0].second == 1 ? 2 : 1);
    }
}

void getReadMetadata(const RunInfoXml&              runInfoXml,
                     const std::vector<size_t>&     barcodeLengths,
                     const std::string&             useBasesMaskStr,
                     common::CycleNumber            read1StartCycle,
                     common::CycleNumber            read2StartCycle,
                     common::CycleNumber            read1EndCycle,
                     common::CycleNumber            read2EndCycle,
                     common::CycleNumber            read1UmiLength,
                     common::CycleNumber            read2UmiLength,
                     common::CycleNumber            read1UmiStartFromCycle,
                     common::CycleNumber            read2UmiStartFromCycle,
                     bool                           trimUmi,
                     layout::LaneInfo&              laneInfo,
                     layout::ReadMetadataContainer& readMetadataContainer)
{
    readMetadataContainer =
        layout::ReadMetadataContainer(runInfoXml.readMetadataBegin(),
                                      runInfoXml.readMetadataEnd());

    size_t indexNumber = 0;
    common::ReadNumber readNumber = 1;
    for (ReadMetadata& readMetadata : readMetadataContainer)
    {
        if (readMetadata.readType_ == common::ReadType::INDEX)
        {
            setIndexReadMetadata(barcodeLengths,
                                 readMetadata,
                                 indexNumber++);
        }
        else
        {
            if (readNumber == 1)
            {
                setReadMetadata(read1StartCycle,
                                read1EndCycle,
                                read1UmiStartFromCycle,
                                read1UmiLength,
                                trimUmi,
                                readMetadata);
            }
            else if (readNumber == 2)
            {
                setReadMetadata(read2StartCycle,
                                read2EndCycle,
                                read2UmiStartFromCycle,
                                read2UmiLength,
                                trimUmi,
                                readMetadata);
            }

            ++readNumber;
        }
    }

    if (!useBasesMaskStr.empty())
    {
        if (read1UmiLength != 0 || read2UmiLength != 0)
        {
            BOOST_THROW_EXCEPTION(common::InputDataError("Use of --use-bases-mask is not supported with a UMI sequence entered in the sample sheet."));
        }

        // useBasesMask overrides the data from RunInfoXml
        layout::UseBasesMask useBasesMask(useBasesMaskStr,
                                          runInfoXml.readMetadataBegin(),
                                          runInfoXml.readMetadataEnd());

        if (!useBasesMask.isEmpty())
        {
            readMetadataContainer = layout::ReadMetadataContainer(useBasesMask.readMetadataBegin(),
                                                                  useBasesMask.readMetadataEnd());

            if (!barcodeLengths.empty())
            {
                validateIndexLengths(readMetadataContainer,
                                     barcodeLengths,
                                     laneInfo);
            }
        }
    }
}

void detectReadLayout(const RunInfoXml&               runInfoXml,
                      const config::SampleSheetCsv&   sampleSheetCsv,
                      const std::vector<std::string>& useBasesMaskStrs,
                      size_t                          minimumTrimmedReadLength,
                      common::CycleNumber             read1StartCycle,
                      common::CycleNumber             read2StartCycle,
                      common::CycleNumber             read1EndCycle,
                      common::CycleNumber             read2EndCycle,
                      common::CycleNumber             read1UmiLength,
                      common::CycleNumber             read2UmiLength,
                      common::CycleNumber             read1UmiStartFromCycle,
                      common::CycleNumber             read2UmiStartFromCycle,
                      bool                            trimUmi,
                      bool                                    autoSetToZeroMismatches,
                      config::BarcodeMismatchCountsContainer& componentMaxMismatches,
                      LaneInfosContainer&                     lanes,
                      const std::vector<std::vector<size_t>>& barcodeLengthsForLane)
{
    std::map<common::LaneNumber, std::string> useBasesMasksForLane;
    std::string defaultUseBasesMask;
    parseUseBasesMasks(useBasesMaskStrs,
                       useBasesMasksForLane,
                       defaultUseBasesMask);

    for (LaneInfo& laneInfo : lanes)
    {
        BarcodeCollisionDetector barcodeCollisionDetector(componentMaxMismatches,
                                                          autoSetToZeroMismatches);

        std::string useBasesMaskStr(defaultUseBasesMask);
        std::map<common::LaneNumber, std::string>::const_iterator pos = useBasesMasksForLane.find(laneInfo.getNumber());
        if (pos != useBasesMasksForLane.end())
        {
            useBasesMaskStr = pos->second;
        }

        BCL2FASTQ_LOG(common::LogLevel::INFO) << "Lane: " << laneInfo.getNumber() << std::endl;
        BCL2FASTQ_LOG(common::LogLevel::INFO) << "Mask: " << useBasesMaskStr << std::endl;

        unsigned int totalAdapters = 0;
        size_t minimumNonIndexReadLength = minimumTrimmedReadLength;

        layout::ReadMetadataContainer readMetadataContainer;
        getReadMetadata(runInfoXml,
                        barcodeLengthsForLane[laneInfo.getNumber()-1],
                        useBasesMaskStr,
                        read1StartCycle,
                        read2StartCycle,
                        read1EndCycle,
                        read2EndCycle,
                        read1UmiLength,
                        read2UmiLength,
                        read1UmiStartFromCycle,
                        read2UmiStartFromCycle,
                        trimUmi,
                        laneInfo,
                        readMetadataContainer);

        validateBarcodes(laneInfo, barcodeCollisionDetector);

        for (const auto& readMetadata : readMetadataContainer)
        {
            // read
            ReadInfo readInfo(readMetadata.readNumber_,
                              readMetadata.readType_,
                              common::CycleRange(readMetadata.firstCycle_, readMetadata.lastCycle_),
                              common::CycleRange(readMetadata.firstUnmaskedCycle_, readMetadata.lastUnmaskedCycle_),
                              readMetadata.readNumber_ == 1 ?
                                  common::CycleRange(read1UmiStartFromCycle, read1UmiLength) :
                                  common::CycleRange(read2UmiStartFromCycle, read2UmiLength));

            // adapters
            if (readInfo.isDataRead())
            {
                size_t readLength = readMetadata.lastCycle_ - readMetadata.firstCycle_ + 1;

                if (readLength < minimumNonIndexReadLength && readLength != 0)
                {
                    minimumNonIndexReadLength = readLength;
                }

                totalAdapters += addAdaptersToRead(sampleSheetCsv.maskAdaptersBegin(), sampleSheetCsv.maskAdaptersEnd(), readInfo, MASK);
                totalAdapters += addAdaptersToRead(sampleSheetCsv.trimAdaptersBegin(), sampleSheetCsv.trimAdaptersEnd(), readInfo, TRIM);
            }

            laneInfo.addRead(readInfo);
        }
        BCL2FASTQ_LOG(common::LogLevel::DEBUG) << "Total number of adapters: " << totalAdapters << std::endl;

        if (minimumTrimmedReadLength > minimumNonIndexReadLength)
        {
            BCL2FASTQ_LOG(common::LogLevel::WARNING) << "Option: '--minimum-trimmed-read-length' with value: "
                                                     << minimumTrimmedReadLength
                                                     << " is being overwritten by the shortest non-index read length: "
                                                     << minimumNonIndexReadLength << std::endl;

            minimumTrimmedReadLength = minimumNonIndexReadLength;
        }

        laneInfo.setMinimumTrimmedReadLength(minimumTrimmedReadLength);
    }
}

typedef std::map<std::string, SampleInfo> SampleMap;

const std::string& getSampleName(const std::string& sampleName, const std::string& sampleId)
{
    return sampleName.empty() ?
        (sampleId.empty() ? SampleInfo::defaultName : sampleId) :
         sampleName;
}

void checkForDuplicateSample(const std::string& sampleId,
                             const SampleInfo& sample1,
                             const SampleInfo& sample2)
{
    if (sample1.getProject() != sample2.getProject())
    {
        BOOST_THROW_EXCEPTION(common::InputDataError((boost::format(
           "Only 1 project is allowed for a given sampleId. SampleId: '" + sampleId + "'.")).str()));
    }
    if (sample1.getName() != sample2.getName())
    {
        BOOST_THROW_EXCEPTION(common::InputDataError((boost::format(
            "Only 1 sample name is allowed for a given sampleId. SampleId: '" + sampleId + "'.")).str()));
    }
}

void handleDuplicateSampleId(SampleMap& sampleMap,
                             SampleInfo& sample)
{
    std::pair<SampleMap::iterator, bool> pos = sampleMap.insert(std::make_pair(sample.getId(), sample));
    if (!pos.second)
    {
        // Duplicate sample id.
        sample.setNumber(pos.first->second.getNumber());

        checkForDuplicateSample(sample.getId(), sample, pos.first->second);
    }
}

void handleDefaultSample(bool hasBarcodes,
                         bool hasNoSamples,
                         size_t numSamplesInLane,
                         LaneInfo& laneInfo,
                         bool& sampleWithoutBarcodesInLane)
{
    if (hasBarcodes)
    {
        if (hasNoSamples)
        {
            SampleInfo defaultSample;
            laneInfo.addSample(defaultSample);
        }
    }
    else
    {
        sampleWithoutBarcodesInLane = true;
    }

    if (sampleWithoutBarcodesInLane && numSamplesInLane != 1)
    {
        BOOST_THROW_EXCEPTION(common::InputDataError((boost::format("Missing a barcode in lane %d") % laneInfo.getNumber()).str()));
    }
}

void addBarcodes(SampleInfo& sampleInfo,
                 const common::SampleMetadata::BarcodesContainer& barcodes,
                 std::vector<size_t>& barcodeLengths)
{
    for (const auto& barcodeComponents : barcodes)
    {
        common::SampleMetadata::BarcodesContainer::value_type nonEmptyBarcodeComponents;
        for (const auto& component : barcodeComponents)
        {
            if (!component.empty())
            {
                nonEmptyBarcodeComponents.push_back(component);
            }
        }
        Barcode barcode (nonEmptyBarcodeComponents.begin(), nonEmptyBarcodeComponents.end());
        sampleInfo.addBarcode(barcode);
    }

    if (barcodeLengths.empty() && !barcodes.empty())
    {
        for (const auto& component : barcodes.front())
        {
            barcodeLengths.push_back(component.size());
        }
    }
}

void addSamplesForLane(const common::SampleMetadataContainer& sampleMetadataContainer,
                       SampleMap& sampleMap,
                       SampleIdNumberMap& sampleIdNumberMap,
                       LaneInfo& laneInfo,
                       std::vector<size_t>& barcodeLengths)
{
    common::LaneNumber laneNumber = laneInfo.getNumber();

    size_t numSamplesInLane = 0;
    bool sampleWithoutBarcodesInLane = false;
    std::map<common::SampleNumber, size_t> laneInfoSampleMap;
    SampleMap sampleMapForLane;
    for (const auto& sampleMetadata : sampleMetadataContainer)
    {
        if (!sampleMetadata.lanes_.empty() &&
            std::find(sampleMetadata.lanes_.begin(), sampleMetadata.lanes_.end(), laneNumber) == sampleMetadata.lanes_.end())
        {
            // SampleSheet.csv specified a lane for this sample, and this lane is not it.
            continue;
        }

        ++numSamplesInLane;

        bool hasBarcodes       = !sampleMetadata.barcodes_.empty();
        std::string sampleId   = sampleMetadata.id_.empty() ? SampleInfo::defaultId : sampleMetadata.id_;
        std::string sampleName = getSampleName(sampleMetadata.name_, sampleMetadata.id_);

        handleDefaultSample(hasBarcodes,
                            laneInfo.sampleInfosBegin() == laneInfo.sampleInfosEnd(),
                            numSamplesInLane,
                            laneInfo,
                            sampleWithoutBarcodesInLane);

        SampleInfo sample(sampleIdNumberMap[sampleId],
                          sampleId,
                          sampleName,
                          sampleMetadata.project_.empty() ? SampleInfo::defaultProject : sampleMetadata.project_);

        handleDuplicateSampleId(sampleMap,
                                sample);

        auto pos = sampleMapForLane.insert(std::make_pair(sample.getId(), sample));
        addBarcodes(pos.first->second,
                    sampleMetadata.barcodes_,
                    barcodeLengths);
    }

    for (const auto& keyValue : sampleMapForLane)
    {
        laneInfo.addSample(keyValue.second);
    }

    // Sort by sample number to improve readability for output
    laneInfo.sortSampleInfos();
}

void getSampleNumbers(const common::SampleMetadataContainer& sampleMetadataContainer,
                      SampleIdNumberMap& sampleIdNumberMap)
{
    common::SampleNumber sampleNumber = 0;
    for (const auto& sampleMetadata : sampleMetadataContainer)
    {
        common::SampleNumber& sampleNumberForId = sampleIdNumberMap[sampleMetadata.id_];
        if (sampleNumberForId == 0)
        {
            sampleNumberForId = ++sampleNumber;
        }
    }
}
void detectLaneLayout(const layout::RunInfoXml&                       runInfoXml,
                      const config::SampleSheetCsv&                   sampleSheetCsv,
                      const boost::filesystem::path&                  inputDir,
                      bool                                            hasConfigXml,
                      bool                                            includeNonPfClusters,
                      bool&                                           hasBci,
                      bool&                                           hasCbcl,
                      LaneInfosContainer&                             lanes,
                      std::vector<std::vector<size_t>>&               barcodeLengthsForLane)
{
    SampleIdNumberMap sampleIdNumberMap;
    getSampleNumbers(sampleSheetCsv.getSampleMetadata(),
                     sampleIdNumberMap);

    SampleMap sampleMap;

    const common::LaneNumber lanesCount = runInfoXml.getLanesCount();
    for (common::LaneNumber laneNumber = 1; laneNumber <= lanesCount; ++laneNumber)
    {
        // lane
        LaneInfo laneInfo(laneNumber);

        barcodeLengthsForLane.resize(barcodeLengthsForLane.size() + 1);
        if (sampleSheetCsv.getSampleMetadata().empty())
        {
            // Add a default sample
            laneInfo.addSample(SampleInfo());
        }
        else
        {
            addSamplesForLane(sampleSheetCsv.getSampleMetadata(),
                              sampleMap,
                              sampleIdNumberMap,
                              laneInfo,
                              barcodeLengthsForLane.back());

        }

        if (laneInfo.sampleInfosBegin() != laneInfo.sampleInfosEnd())
        {
            hasBci = hasBci || checkBciExistence( inputDir / laneInfo.getDirName() );
            hasCbcl = hasCbcl || checkCbclExistence( inputDir / laneInfo.getDirName() );

            if (hasCbcl && includeNonPfClusters)
            {
                BOOST_THROW_EXCEPTION(common::InputDataError("Cannot specify '--with-failed-reads' with cbcl input data. Failed reads are not stored in cbcl files."));
            }

            lanes.push_back(laneInfo);
        }
    }

    BCL2FASTQ_LOG(common::LogLevel::DEBUG) << "Has BCI: " << hasBci << std::endl;
    BCL2FASTQ_LOG(common::LogLevel::DEBUG) << "Has 'config.xml': " << hasConfigXml << std::endl;
}

void readTilesFromRunInfoXmlForLane(const layout::RunInfoXml&        runInfoXml,
                                    const config::SampleSheetCsv&    sampleSheet,
                                    const std::vector<boost::regex>& tileRegexps,
                                    common::LaneNumber               laneNumber,
                                    LaneInfo&                        laneInfo)
{
    size_t skippedTiles = 0;
    size_t tileIndex = 0;
    BOOST_FOREACH (const common::TileNumber &tileNumber,
                   std::make_pair( runInfoXml.tileNumbersBegin( laneNumber ),
                                   runInfoXml.tileNumbersEnd( laneNumber )) )
    {
        if ( includeTile(runInfoXml.getFlowcellId(), laneNumber, tileNumber, tileRegexps, sampleSheet) )
        {
            TileInfo tileInfo(tileNumber, tileIndex, skippedTiles);
            laneInfo.addTile(tileInfo);
            skippedTiles = 0;
        }
        else
        {
            ++skippedTiles;
        }
        ++tileIndex;
    }
}

void readTilesFromRunInfoXml(const layout::RunInfoXml&        runInfoXml,
                             const config::SampleSheetCsv&    sampleSheet,
                             const std::vector<boost::regex>& tileRegexps,
                             std::vector<LaneInfo>&           lanes)
{
    BOOST_FOREACH (LaneInfo &laneInfo, std::make_pair(lanes.begin(), lanes.end()))
    {
        common::LaneNumber laneNumber = laneInfo.getNumber();
        BCL2FASTQ_ASSERT_MSG( 1 <= laneNumber && laneNumber <= runInfoXml.getLanesCount(), "Lane '" << laneNumber << "' is out of range");
        BCL2FASTQ_LOG(common::LogLevel::INFO) << "Lane: " << laneNumber << std::endl;

        readTilesFromRunInfoXmlForLane(runInfoXml,
                                       sampleSheet,
                                       tileRegexps,
                                       laneNumber,
                                       laneInfo);
    }
}

void readTilesFromBciFile(const boost::filesystem::path&   inputDir,
                          const config::SampleSheetCsv&    sampleSheet,
                          const std::vector<boost::regex>& tileRegexps,
                          const std::string&               flowcellId,
                          std::vector<LaneInfo>&           lanes)
{
    BOOST_FOREACH (LaneInfo &laneInfo, std::make_pair(lanes.begin(), lanes.end()))
    {
        size_t skippedTiles = 0;
        size_t skippedClusters = 0;
        size_t tileIndex = 0;
        common::LaneNumber laneNumber = laneInfo.getNumber();
        if (checkBciFile(inputDir, laneNumber))
        {
            const BCIndex bcIndex = createBCIndex(inputDir, laneInfo.getNumber());
            BOOST_FOREACH (const BCIndex::TileMetadata &tileMetadata, std::make_pair(bcIndex.tileMetadataBegin(), bcIndex.tileMetadataEnd()))
            {
                if ( includeTile(flowcellId, laneNumber, tileMetadata.tileNumber_, tileRegexps, sampleSheet) &&
                     tileMetadata.clustersCount_ > 0)
                {
                    TileInfo tileInfo(tileMetadata.tileNumber_, tileIndex, tileMetadata.clustersCount_, skippedTiles, skippedClusters);
                    laneInfo.addTile(tileInfo);
                    skippedTiles = 0;
                    skippedClusters = 0;
                }
                else
                {
                    ++skippedTiles;
                    skippedClusters += tileMetadata.clustersCount_;
                }
                ++tileIndex;
            }
        }
    }
}

void readTilesFromConfigXmlOrRunInfoXml(const boost::filesystem::path&   inputDir,
                                        const layout::RunInfoXml&        runInfoXml,
                                        const config::SampleSheetCsv&    sampleSheet,
                                        const std::vector<boost::regex>& tileRegexps,
                                        std::vector<LaneInfo>&           lanes)
{
    BCL2FASTQ_LOG(common::LogLevel::DEBUG) << "This may be a 'HiSeq 2000' / 'MiSeq' run" << std::endl;
    const ConfigXml configXml = createConfigXml(inputDir);
    BOOST_FOREACH (LaneInfo &laneInfo, std::make_pair(lanes.begin(), lanes.end()))
    {
        size_t skippedTiles = 0;
        size_t tileIndex = 0;
        common::LaneNumber laneNumber = laneInfo.getNumber();
        BCL2FASTQ_LOG(common::LogLevel::INFO) << "Lane: " << laneNumber << std::endl;

        if (configXml.containsLaneNumber(laneNumber))
        {
            BOOST_FOREACH (const ConfigXml::TileMetadata &tileMetadata,
                           std::make_pair(configXml.tileMetadataBegin(laneInfo.getNumber()),
                                          configXml.tileMetadataEnd(laneInfo.getNumber())) )
            {
                if ( includeTile(runInfoXml.getFlowcellId(), laneNumber, tileMetadata.tileNumber_, tileRegexps, sampleSheet) )
                {
                    TileInfo tileInfo(tileMetadata.tileNumber_, tileIndex, skippedTiles);
                    laneInfo.addTile(tileInfo);
                    skippedTiles = 0;
                }
                else
                {
                    ++skippedTiles;
                }
                ++tileIndex;
            }
        }
        else
        {
            readTilesFromRunInfoXmlForLane(runInfoXml,
                                           sampleSheet,
                                           tileRegexps,
                                           laneNumber,
                                           laneInfo);
        }
    }
}

void TileLayoutDetector::readTilesFromCbclHeaders(
        const boost::filesystem::path&   inputDir,
        const config::SampleSheetCsv&    sampleSheet,
        const std::vector<boost::regex>& tileRegexps,
        const std::string&               flowcellId,
        std::vector<LaneInfo>&           lanes,
        common::TileFileMap&             tileFileMap,
        common::NumBasesPerByte&         numBasesPerByte,
        common::CycleNumber              read1StartCycle,
        common::CycleNumber              read1EndCycle,
        bool                             ignoreMissingBcls)
{
    tileFileMap.resize(lanes.back().getNumber() + 1);

    // This method reads the bcl files for one cycle to get a list of tiles. All cycles contain
    // the same set of tiles. If the tiles cannot be read, bcl2fastq will abort. To circumvent
    // aborting on cycle 1 if corrupt, retry with a few subsequent cycles.
    //
    const int maxCyclesToAttempt = 3;

    uint32_t startCycle = (read1StartCycle > 0) ? read1StartCycle : 1;

    for (auto& lane : lanes)
    {
        std::vector<std::pair<common::TileNumber, boost::filesystem::path>> tileNumbers;
        auto laneNumber = lane.getNumber();

        int cycleAttempt = 0;
        int success = false;
        while (!success) {

            int cycle = startCycle+cycleAttempt;
            ++cycleAttempt;
            std::string cycleFolder = "C" + std::to_string(cycle) + ".1";

            std::vector<boost::filesystem::path> files;
            checkExistence(inputDir / lane.getDirName() / cycleFolder,
                           "cbcl",
                           files,
                           false);

            if (files.size() == 0) {

                std::stringstream ss;
                ss << "No cbcl files found for cycle " << cycle << ". "
                   << "Attempt " << cycle << " of " << maxCyclesToAttempt << ".";

                if (cycleAttempt < maxCyclesToAttempt) {
                    BCL2FASTQ_LOG(common::LogLevel::WARNING)
                        << ss.str() << " Retry with cycle " << cycle+1 << "." << std::endl;
                }
                else {
                    ss << "\nCorrupt or missing .cbcl files were found for 3 consecutive cycles, verify the source of the data.";
                    BOOST_THROW_EXCEPTION(common::InputDataError("\n" + ss.str()));
                }

                continue;
            }

            tileNumbers.resize(0);
            for (const auto& fileName : files)
            {
                try {
                    io::SyncFile cbclFile(fileName,
                                          false, // Can't ignore errors here. Need this data.
                                          std::ios_base::in | std::ios_base::binary);

                    // Synchronization isn't really necessary here, but its safer to
                    // always require it.
                    io::SyncFile::SyncFileReader cbclFileReader(cbclFile, 0);

                    data::CbclFileReader::Header cbclHeader;
                    if (!cbclHeader.load(cbclFileReader, false))
                    {
                        cbclHeader.reset();
                    }

                    // 8 bits per byte, 2 bits per basecall
                    numBasesPerByte = static_cast<common::NumBasesPerByte>(8 / (cbclHeader.getNumberOfBitsPerQscore() + 2));

                    for (auto tileNumber : cbclHeader.getTileIds())
                    {
                        tileNumbers.push_back(std::make_pair(tileNumber, fileName.leaf()));
                    }

                    success = true;
                }
                catch(...) {

                    std::stringstream ss;
                    ss << "Found corrupt cycle " << cycle << " CBCL file during initialization. "
                       << "Attempt " << cycleAttempt << " of " << maxCyclesToAttempt << " for file\n"
                       << fileName.string();

                    if (cycleAttempt < maxCyclesToAttempt) {
                        BCL2FASTQ_LOG(common::LogLevel::WARNING)
                            << ss.str() << " Retry with cycle " << cycle+1 << "." << std::endl;
                        break;
                    }
                    else {
                        ss << "\nCorrupt or missing .cbcl files were found for 3 consecutive cycles, verify the source of the data.";
                        BOOST_THROW_EXCEPTION(common::InputDataError("\n" + ss.str()));
                    }
                }

            }
        } // retry if needed.

        std::sort(tileNumbers.begin(), tileNumbers.end());

        size_t skippedTiles = 0;
        size_t tileIndex = 0;
        for (const auto& tileNumberFilePair : tileNumbers)
        {
            if (includeTile(flowcellId, laneNumber, tileNumberFilePair.first, tileRegexps, sampleSheet) )
            {
                TileInfo tileInfo(tileNumberFilePair.first, tileIndex, skippedTiles);
                lane.addTile(tileInfo);
                skippedTiles = 0;

                tileFileMap[laneNumber].insert(tileNumberFilePair);
            }
            else
            {
                ++skippedTiles;
            }

            ++tileIndex;
        }
    }
}

void TileLayoutDetector::detectTileLayout(const std::vector<std::string>& tilesFilterList,
                      const layout::RunInfoXml&       runInfoXml,
                      const config::SampleSheetCsv&   sampleSheet,
                      const boost::filesystem::path&  inputDir,
                      std::vector<LaneInfo>&          lanes,
                      common::TileFileMap&            tileFileMap,
                      common::NumBasesPerByte&        numBasesPerByte,
                      common::TileAggregationMode     flowcellAggregateTilesFlag,
                      bool                            hasConfigXml,
                      common::CycleNumber             read1StartCycle,
                      common::CycleNumber             read1EndCycle,
                      bool                            ignoreMissingBcls)
{
    std::vector<boost::regex> tileRegexps;
    BOOST_FOREACH(const std::string &tileRegex, tilesFilterList)
    {
        tileRegexps.push_back(boost::regex(tileRegex));
    }

    switch (flowcellAggregateTilesFlag)
    {
    case common::TileAggregationMode::NON_AGGREGATED:
        if( hasConfigXml )
        {
            readTilesFromConfigXmlOrRunInfoXml(inputDir,
                                               runInfoXml,
                                               sampleSheet,
                                               tileRegexps,
                                               lanes);
        } else {  // !hasConfigXml
            readTilesFromRunInfoXml(runInfoXml,
                                    sampleSheet,
                                    tileRegexps,
                                    lanes);
        }
        break;

    case common::TileAggregationMode::AGGREGATED:
        readTilesFromBciFile(inputDir,
                             sampleSheet,
                             tileRegexps,
                             runInfoXml.getFlowcellId(),
                             lanes);
        break;
    case common::TileAggregationMode::CBCL:
        readTilesFromCbclHeaders(inputDir,
                                 sampleSheet,
                                 tileRegexps,
                                 runInfoXml.getFlowcellId(),
                                 lanes,
                                 tileFileMap,
                                 numBasesPerByte,
                                 read1StartCycle,
                                 read1EndCycle,
                                 ignoreMissingBcls);
        break;
    default:
        BCL2FASTQ_ASSERT_MSG(false, "Unrecognized tile aggregation mode");
    }

    bool hasTiles = false;
    for (const auto& laneInfo : lanes)
    {
        if (!laneInfo.getTileInfos().empty())
        {
            hasTiles = true;
            break;
        }
    }

    if (!hasTiles)
    {
        bool first = true;
        std::string tileFilterStr;
        for (const auto& filter : tilesFilterList)
        {
            if (first)
            {
                first = false;
            }
            else
            {
                tileFilterStr += ", ";
            }
            tileFilterStr += "'" + filter + "'";
        }
        BOOST_THROW_EXCEPTION(common::InputDataError("No tiles were found to process." +
            (tilesFilterList.empty() ? "" : (" Tile regular expressions used: (" + tileFilterStr + ")"))));
    }
}

Layout::Layout(
    const boost::filesystem::path& intensitiesDir,
    const boost::filesystem::path& inputDir,
    const boost::filesystem::path& outputDir,
    const boost::filesystem::path& reportsDir,
    const boost::filesystem::path& statsDir,
    const config::SampleSheetCsv& sampleSheet,
    const RunInfoXml& runInfoXml,
    const std::vector<std::string> & tilesFilterList,
    const std::vector<std::string>& useBasesMasks,
    size_t minimumTrimmedReadLength,
    bool autoSetToZeroMismatches,
    config::BarcodeMismatchCountsContainer& componentMaxMismatches,
    bool ignoreMissingBcls,
    bool ignoreMissingFilters,
    bool ignoreMissingPositions,
    common::CycleNumber read1StartCycle,
    common::CycleNumber read2StartCycle,
    common::CycleNumber read1EndCycle,
    common::CycleNumber read2EndCycle,
    common::CycleNumber read1UmiLength,
    common::CycleNumber read2UmiLength,
    common::CycleNumber read1UmiStartFromCycle,
    common::CycleNumber read2UmiStartFromCycle,
    bool trimUmi,
    bool includeNonPfClusters
)
: flowcellInfo_()
, laneInfos_()
, tileFileMap_()
, numBasesPerByte_(common::NumBasesPerByte::ONE)
{
    common::DirectoryValidator::getSingleton().addUniquePath(reportsDir, "Reports");
    common::DirectoryValidator::getSingleton().addUniquePath(statsDir, "Stats");

    bool hasConfigXml = checkConfigXml(inputDir);

    bool hasBci = false;
    bool hasCbcl = false;
    std::vector<LaneInfo> lanes;
    std::vector<std::vector<size_t>> barcodeLengthsForLane;
    detectLaneLayout(runInfoXml,
                     sampleSheet,
                     inputDir,
                     hasConfigXml,
                     includeNonPfClusters,
                     hasBci,
                     hasCbcl,
                     lanes,
                     barcodeLengthsForLane);

    detectReadLayout(runInfoXml,
                     sampleSheet,
                     useBasesMasks,
                     minimumTrimmedReadLength,
                     read1StartCycle,
                     read2StartCycle,
                     read1EndCycle,
                     read2EndCycle,
                     read1UmiLength,
                     read2UmiLength,
                     read1UmiStartFromCycle,
                     read2UmiStartFromCycle,
                     trimUmi,
                     autoSetToZeroMismatches,
                     componentMaxMismatches,
                     lanes,
                     barcodeLengthsForLane);

    detectFlowcellInfo(intensitiesDir,
                       runInfoXml,
                       hasBci,
                       hasCbcl);

    TileLayoutDetector::detectTileLayout(
        tilesFilterList,
        runInfoXml,
        sampleSheet,
        inputDir,
        lanes,
        tileFileMap_,
        numBasesPerByte_,
        flowcellInfo_.getAggregateTilesMode(),
        hasConfigXml,
        read1StartCycle,
        read1EndCycle,
        ignoreMissingBcls);

    for (const auto& laneInfo : lanes)
    {
         if (laneInfo.getTileInfos().empty())
         {
             continue;
         }

         FileExistenceVerifier::verifyAllFilesExist(inputDir,
                                                    intensitiesDir,
                                                    laneInfo,
                                                    tileFileMap_,
                                                    flowcellInfo_.getAggregateTilesMode(),
                                                    flowcellInfo_.isPatternedFlowcell(),
                                                    ignoreMissingBcls,
                                                    ignoreMissingFilters,
                                                    ignoreMissingPositions);

        laneInfos_.push_back(laneInfo);
    }
}


} // namespace layout


} // namespace bcl2fastq


