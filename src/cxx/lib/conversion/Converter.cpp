/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file Converter.cpp
 *
 * \brief Implementation of BCL to FASTQ converter.
 *
 * \author Marek Balint
 * \author Mauricio Varea
 */


#include <boost/ref.hpp>
#include <boost/bind.hpp>

#include <boost/assign/list_of.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/range/adaptor/reversed.hpp>

#include "common/Logger.hh"
#include "common/Exceptions.hh"
#include "common/FileSystem.hh"
#include "stats/BarcodeStats.hh"
#include "stats/TileStats.hpp"
#include "stats/DemultiplexingStatsXml.hh"
#include "stats/ConversionStatsXml.hh"
#include "stats/DemuxReportGenerator.hh"
#include "stats/Json.hh"

#include "conversion/Converter.hh"


namespace bcl2fastq
{
namespace conversion
{


Converter::Converter(
    const bcl2fastq::config::Bcl2FastqOptions &options,
    const layout::Layout &layout,
    const layout::LaneInfo &laneInfo,
    bcl2fastq::conversion::ConverterStats &allStats
)
: layout_(layout)
, laneInfo_(laneInfo)
, statsFinalizer_(laneInfo.getSampleInfos().size()==1,
                  allStats.getConversionStats(laneInfo_),
                  allStats.getDemuxStats(laneInfo_))
, allStats_(allStats)
, inUseRawBclBufferQueue_(MAX_BCL_BUFFER_QUEUE_DEPTH)
, recycledRawBclBufferQueue_()
, inUseBclBufferQueue_(MAX_BCL_BUFFER_QUEUE_DEPTH)
, recycledBclBufferQueue_()
, inUseDemuxBufferQueue_()
, inUseFastqBufferMap_()
, recycledFastqBufferQueue_()
, processTaskQueue_()
, ioReadTaskQueue_()
, ioWriteTaskQueue_()
, terminator_()
, barcodeTranslationTable_(
    laneInfo.sampleInfosBegin(),
    laneInfo.sampleInfosEnd(),
    options.getBarcodeMismatches()
  )
, bclReader_(
    ioReadTaskQueue_,
    inUseRawBclBufferQueue_,
    recycledRawBclBufferQueue_,
    layout,
    laneInfo,
    options.getIgnoreMissingBcls(),
    options.getIgnoreMissingFilters(),
    options.getIgnoreMissingPositions(),
    options.getIgnoreMissingControls(),
    options.getInputDir(),
    options.getIntensitiesDir()
  )
, bclLoader_(
    processTaskQueue_,
    inUseRawBclBufferQueue_,
    recycledRawBclBufferQueue_,
    inUseBclBufferQueue_,
    recycledBclBufferQueue_,
    layout,
    laneInfo,
    options.getIgnoreMissingBcls(),
    options.getIgnoreMissingFilters(),
    options.getIgnoreMissingPositions(),
    options.getIgnoreMissingControls(),
    bclReader_.getNumTilesPerBuffer()
  )
, demultiplexer_(
    processTaskQueue_,
    inUseBclBufferQueue_,
    inUseDemuxBufferQueue_,
    layout,
    laneInfo,
    barcodeTranslationTable_,
    options.getIncludeNonPfClusters(),
    options.getFastqCreatorThreadsCount(),
    allStats.getDemuxStats(laneInfo)
  )
, fastqCreator_(
    createNewFastqCreator(
        processTaskQueue_,
        inUseDemuxBufferQueue_,
        recycledBclBufferQueue_,
        inUseFastqBufferMap_,
        recycledFastqBufferQueue_,
        layout,
        laneInfo,
        options,
        allStats))
, fastqWriter_(
    options.getFastqWriterThreadsCount() == 0 ? ioReadTaskQueue_ : ioWriteTaskQueue_,
    inUseFastqBufferMap_,
    recycledFastqBufferQueue_,
    layout,
    laneInfo,
    options.getCreateFastqsForIndexReads(),
    options.getOutputDir(),
    options.noLaneSplitting()
  )
, processThreadPool_(options.getFastqCreatorThreadsCount(),
                     processTaskQueue_,
                     terminator_)
, ioReadThreadPool_(options.getBclLoaderThreadsCount(),
                    ioReadTaskQueue_,
                    terminator_)
, ioWriteThreadPool_(options.getFastqWriterThreadsCount(),
                     ioWriteTaskQueue_,
                     terminator_)
{
}

std::unique_ptr<Stage> Converter::createNewFastqCreator(TaskQueue& processTaskQueue,
                                                        ThreadSafeQueue<std::shared_ptr<DemuxBuffer>>& inUseDemuxBufferQueue,
                                                        ThreadSafeQueue<std::shared_ptr<DemuxBuffer>>& recycledBclBufferQueue,
                                                        ThreadSafeMap<uint32_t, std::shared_ptr<FastqBuffer>>& inUseFastqBufferMap,
                                                        ThreadSafeQueue<std::shared_ptr<FastqBuffer>>& recycledFastqBufferQueue,
                                                        const layout::Layout& layout,
                                                        const layout::LaneInfo& laneInfo,
                                                        const bcl2fastq::config::Bcl2FastqOptions& options,
                                                        bcl2fastq::conversion::ConverterStats& allStats)
{
    // The purpose of the FastqCreator being a template on the NumBasesPerByte is to avoid
    // costly branch logic inside of FastqIterator. Surprisingly, including even a simple switch
    // statement inside of the FastqIterator slows down all of bcl2fastq by 30%. This makes some
    // sense if you consider that we can be iterating over trillions of bases. Also, the use
    // of branches inside of FastqIterator makes vector AVX instructions impossible.
    switch (layout.getNumBasesPerByte())
    {
    case common::NumBasesPerByte::ONE:
        return std::unique_ptr<FastqCreator<common::NumBasesPerByte::ONE>>(
            new FastqCreator<common::NumBasesPerByte::ONE>(
                processTaskQueue,
                inUseDemuxBufferQueue,
                recycledBclBufferQueue,
                inUseFastqBufferMap,
                recycledFastqBufferQueue,
                layout,
                laneInfo,
                options.getMaskShortAdapterReads(),
                options.getAdapterStringency(),
                options.getGenerateReverseComplementFastqs(),
                options.getIncludeNonPfClusters(),
                options.getCreateFastqsForIndexReads(),
                options.useBgzfCompression(),
                options.getFastqCompressionLevel(),
                options.findAdaptersWithSlidingWindow(),
                options.getFastqCreatorThreadsCount(),
                allStats.getConversionStats(laneInfo),
                allStats.getUnknownBarcodeStat(laneInfo)));
    case common::NumBasesPerByte::TWO:
        return std::unique_ptr<FastqCreator<common::NumBasesPerByte::TWO>>(
            new FastqCreator<common::NumBasesPerByte::TWO>(
                processTaskQueue,
                inUseDemuxBufferQueue,
                recycledBclBufferQueue,
                inUseFastqBufferMap,
                recycledFastqBufferQueue,
                layout,
                laneInfo,
                options.getMaskShortAdapterReads(),
                options.getAdapterStringency(),
                options.getGenerateReverseComplementFastqs(),
                options.getIncludeNonPfClusters(),
                options.getCreateFastqsForIndexReads(),
                options.useBgzfCompression(),
                options.getFastqCompressionLevel(),
                options.findAdaptersWithSlidingWindow(),
                options.getFastqCreatorThreadsCount(),
                allStats.getConversionStats(laneInfo),
                allStats.getUnknownBarcodeStat(laneInfo)));
    case common::NumBasesPerByte::FOUR:
        return std::unique_ptr<FastqCreator<common::NumBasesPerByte::FOUR>>(
            new FastqCreator<common::NumBasesPerByte::FOUR>(
                processTaskQueue,
                inUseDemuxBufferQueue,
                recycledBclBufferQueue,
                inUseFastqBufferMap,
                recycledFastqBufferQueue,
                layout,
                laneInfo,
                options.getMaskShortAdapterReads(),
                options.getAdapterStringency(),
                options.getGenerateReverseComplementFastqs(),
                options.getIncludeNonPfClusters(),
                options.getCreateFastqsForIndexReads(),
                options.useBgzfCompression(),
                options.getFastqCompressionLevel(),
                options.findAdaptersWithSlidingWindow(),
                options.getFastqCreatorThreadsCount(),
                allStats.getConversionStats(laneInfo),
                allStats.getUnknownBarcodeStat(laneInfo)));
    default:
        BCL2FASTQ_ASSERT_MSG(false, "Unsuported number of bases per byte. Supported values include: One, Two, or Four");
    };

    return std::unique_ptr<Stage>();
}

void Converter::run()
{
    StageRefsContainer stages;

    stages.push_back(StageRefsContainer::value_type(bclReader_));
    stages.push_back(StageRefsContainer::value_type(bclLoader_));
    stages.push_back(StageRefsContainer::value_type(demultiplexer_));
    stages.push_back(StageRefsContainer::value_type(*fastqCreator_));
    stages.push_back(StageRefsContainer::value_type(fastqWriter_));

    terminator_.add(&bclReader_);
    terminator_.add(&bclLoader_);
    terminator_.add(&demultiplexer_);
    terminator_.add(fastqCreator_.get());
    terminator_.add(&fastqWriter_);

    common::ThreadVector stageThreads(stages.size());

    stageThreads.execute(boost::bind(&Converter::stageThreadFunction, this, boost::ref(stages), _1));

    // Done adding tasks.
    processTaskQueue_.setFinished();
    ioReadTaskQueue_.setFinished();
    ioWriteTaskQueue_.setFinished();
}

Converter::StatsFinalizer::StatsFinalizer(bool isOnly1Sample,
                                          const ConversionStatsForSampleTileBarcode& conversionStats,
                                          DemultiplexTask::DemuxStats& demuxStats)
: demuxStats_(demuxStats)
, conversionStats_(conversionStats)
, isOnly1Sample_(isOnly1Sample)
{
}

Converter::StatsFinalizer::~StatsFinalizer()
{
    if (isOnly1Sample_)
    {
        // 1 barcodeStats per tile
        size_t tileNum = 0;
        for (auto& barcodeStats : demuxStats_.front().front())
        {
            barcodeStats.setBarcodeCount(conversionStats_.front().at(tileNum).front().first[stats::TileBarcodeStats::PF].clustersCount);
            ++tileNum;
        }
    }
}

void Converter::stageThreadFunction(StageRefsContainer &stages, common::ThreadVector::size_type threadNum)
{
    BCL2FASTQ_LOG(common::LogLevel::DEBUG) << "Running stage #" << threadNum << std::endl;
    stages.at(threadNum).get().run();
    BCL2FASTQ_LOG(common::LogLevel::DEBUG) << "Stage #" << threadNum << " complete" << std::endl;
}


ConverterStats::ConverterStats( const layout::Layout &layout,
                                const boost::filesystem::path &fastqSummaryPartialFileName,
                                const boost::filesystem::path &demuxSummaryPartialFileName,
                                const boost::filesystem::path &adapterTrimmingTxtFile,
                                const boost::filesystem::path &interopFile,
                                const boost::filesystem::path &demuxStatsXmlFile,
                                const boost::filesystem::path &conversionStatsXmlFile,
                                const boost::filesystem::path &statsJsonFile)
: layout_(layout)
, fastqSummaryPartialFileName_(fastqSummaryPartialFileName)
, demuxSummaryPartialFileName_(demuxSummaryPartialFileName)
, adapterTrimmingTxtFile_(adapterTrimmingTxtFile)
, interopFile_(interopFile)
, demuxStatsXmlFile_(demuxStatsXmlFile)
, conversionStatsXmlFile_(conversionStatsXmlFile)
, statsJsonFile_(statsJsonFile)
, demuxStats_()
, conversionStats_()
, topUnknownBarcodes_()
, barcodeLaneStatsMap_()
{
    boost::filesystem::path statsPath0 = interopFile.parent_path();
    boost::filesystem::path statsPath1 = demuxStatsXmlFile.parent_path();
    boost::filesystem::path statsPath2 = conversionStatsXmlFile.parent_path();
    BCL2FASTQ_ASSERT_MSG(statsPath1 == statsPath2, "All human-readable XML statistics need to live in the same directory");
    std::vector<boost::filesystem::path> createList =
                boost::assign::list_of(statsPath0)(statsPath1);
    common::createDirectories(createList);
}

void ConverterStats::deserializeMetadata( const layout::LaneInfo &laneInfo,
                                          size_t sampleIndex,
                                          size_t barcodeIndex,
                                          std::string &indexName,
                                          std::string &sampleId,
                                          std::string &projectName,
                                          std::string &sampleName,
                                          common::SampleNumber &sampleNumber) const
{
    const layout::LaneInfo::SampleInfosContainer::value_type &sampleInfo = laneInfo.getSampleInfos().at(sampleIndex);
    indexName = sampleInfo.hasBarcodes() ? sampleInfo.getBarcode(barcodeIndex).toString() : layout::Barcode::DEFAULT_BARCODE;
    sampleId = sampleInfo.getId();
    projectName = sampleInfo.getProject();
    sampleName = sampleInfo.getName();
    sampleNumber = sampleInfo.getNumber();
}

void ConverterStats::dumpInteropStats() const
{
    data::InteropFile interopFile;
    if (!interopFile.open(interopFile_)) {
        return;
    }

    data::InteropRecord interopRecord;
    std::vector<char> buffer;  ///< \brief Intermediate buffer for demultiplexing statistics output.

    const std::string flowcellId = layout_.getFlowcellInfo().getFlowcellId();

    for (const auto& stats : demuxStats_)
    {
        layout::LaneInfo laneInfo = stats.first;
        interopRecord.laneNumber = laneInfo.getNumber();
        if( stats.second.empty() )
        {
            BCL2FASTQ_LOG(common::LogLevel::WARNING) << "No InterOp stats for lane '" << interopRecord.laneNumber << "'" << std::endl;
            continue;
        }
        for (const auto& summaryBarcodeStatsForSample : stats.second)
        {
            for (const auto& summaryBarcodeStatsForBarcode : summaryBarcodeStatsForSample)
            {
                for (const auto& summaryBarcodeStats : summaryBarcodeStatsForBarcode)
                {
                    // skip sample #0 (i.e. the default sample)
                    if( 0 == summaryBarcodeStats.index.first.sampleIndex_ )
                    {
                        continue;
                    }
                    std::string sampleName;
                    common::SampleNumber sampleNumber;
                    deserializeMetadata(laneInfo,
                                        summaryBarcodeStats.index.first.sampleIndex_,
                                        summaryBarcodeStats.index.first.barcodeIndex_,
                                        interopRecord.indexName,
                                        interopRecord.sampleId,
                                        interopRecord.projectName,
                                        sampleName, sampleNumber);

                    std::replace(interopRecord.indexName.begin(), interopRecord.indexName.end(), '+', '-');

                    std::string flowcellId = layout_.getFlowcellInfo().getFlowcellId();
                    interopRecord.tileNumber = laneInfo.getTileInfos()[summaryBarcodeStats.index.second].getNumber();
                    interopRecord.clustersCount = summaryBarcodeStats.getBarcodeCount();

                    for (const auto& readInfo : laneInfo.readInfos())
                    {
                        if (!readInfo.isDataRead())
                        {
                            continue;
                        }
                        interopRecord.readNumber = readInfo.getNumber();

                        interopFile.writeInteropRecord(interopRecord);
                    }
                }
            }
        }
    }
}

void ConverterStats::dumpDemuxStats()
{
    stats::DemultiplexingStatsXml demuxStatsXml;
    const std::string flowcellId = layout_.getFlowcellInfo().getFlowcellId();
    for (const DemuxPair& stats : demuxStats_)
    {
        const layout::LaneInfo& laneInfo = stats.first;
        unsigned int laneNumber = laneInfo.getNumber();
        if( stats.second.empty() )
        {
            BCL2FASTQ_LOG(common::LogLevel::WARNING) << "No demultiplexing stats for lane '" << laneNumber << "'" << std::endl;
            continue;
        }

        std::map<std::string, std::string> uniqueSampleNameMap;
        stats::DemuxReportGenerator::makeUniqueSampleNameMap(laneInfo, uniqueSampleNameMap);

        FlowcellProjectSampleLaneBarcodeStats flowcellProjectSampleStats;
        std::vector<stats::BarcodeStats> summaryBarcodeStats(0);
        for (const std::vector<std::vector<bcl2fastq::stats::BarcodeStats> >& barcodeStatsForSample : stats.second)
        {
            for (const std::vector<bcl2fastq::stats::BarcodeStats>& barcodeStatsForBarcode : barcodeStatsForSample)
            {
                summaryBarcodeStats.push_back(barcodeStatsForBarcode.front().index);

                for (const stats::BarcodeStats& barcodeStats : barcodeStatsForBarcode)
                {
                    summaryBarcodeStats.back() += barcodeStats;
                }
            }
        }

        barcodeLaneStatsMap_.insert(std::make_pair(laneNumber, BarcodeStatsMap()));

        std::string indexName,sampleId,projectName,sampleName;
        common::SampleNumber sampleNumber;
        BOOST_FOREACH (const stats::BarcodeStats &barcodeStats, std::make_pair(summaryBarcodeStats.begin(), summaryBarcodeStats.end()))
        {
            deserializeMetadata(laneInfo,
                                barcodeStats.index.first.sampleIndex_,
                                barcodeStats.index.first.barcodeIndex_,
                                indexName,
                                sampleId,
                                projectName,
                                sampleName,
                                sampleNumber);

            flowcellProjectSampleStats[flowcellId][projectName][sampleId] += barcodeStats;
            flowcellProjectSampleStats[flowcellId][projectName]["all"] += barcodeStats;
            flowcellProjectSampleStats[flowcellId]["all"]["all"] += barcodeStats;
            demuxStatsXml.addLaneBarcode(flowcellId, projectName, uniqueSampleNameMap[sampleId], indexName, laneNumber, barcodeStats);
            barcodeLaneStatsMap_[laneNumber].insert(std::make_pair(indexName, barcodeStats));
        }
        BOOST_FOREACH(const FlowcellProjectSampleLaneBarcodeStats::value_type &flowcellStats, flowcellProjectSampleStats)
        {
            BOOST_FOREACH(const ProjectSampleLaneBarcodeStats::value_type &projectStats, flowcellStats.second)
            {
                BOOST_FOREACH(const SampleLaneBarcodeStats::value_type &sampleStats, projectStats.second)
                {
                    demuxStatsXml.addLaneBarcode(flowcellStats.first, projectStats.first, sampleStats.first == "all" ? "all" : uniqueSampleNameMap[sampleStats.first], "all", laneNumber, sampleStats.second);
                }
            }
        }
    }

    std::ofstream os(demuxStatsXmlFile_.string().c_str());
    if (!os) {
        BOOST_THROW_EXCEPTION(common::IoError(errno, "Unable to open file for writing: " + demuxStatsXmlFile_.string()));
    }
    if (!(os << demuxStatsXml)) {
        BOOST_THROW_EXCEPTION(common::IoError(errno, "Failed to store demultiplexing statistics in : " + demuxStatsXmlFile_.string()));
    }
}

void ConverterStats::dumpConversionStats()
{
    stats::JsonObject jsonStats;

    stats::DemultiplexingStatsXml demuxStatsXml;
    const std::string flowcellId = layout_.getFlowcellInfo().getFlowcellId();

    jsonStats.add("Flowcell", std::shared_ptr<stats::JsonType>(new stats::JsonString(flowcellId)));

    try
    {
        int64_t runNumber = boost::lexical_cast<int64_t>(layout_.getFlowcellInfo().getRunNumber());
        jsonStats.add("RunNumber", std::shared_ptr<stats::JsonType>(new stats::JsonNumber(runNumber)));
    }
    catch (boost::bad_lexical_cast&)
    {
        // Unable to get the RunNumber from RunInfo.xml
        jsonStats.add("RunNumber", std::shared_ptr<stats::JsonType>(new stats::JsonNumber()));
    }

    jsonStats.add("RunId", std::shared_ptr<stats::JsonType>(new stats::JsonString(layout_.getFlowcellInfo().getRunId())));

    stats::JsonArray* jsonReadInfosForLanes = NULL;
    jsonStats.add("ReadInfosForLanes", std::shared_ptr<stats::JsonType>(jsonReadInfosForLanes = new stats::JsonArray));

    stats::JsonArray* jsonConversionResults = NULL;
    jsonStats.add("ConversionResults", std::shared_ptr<stats::JsonType>(jsonConversionResults = new stats::JsonArray));

    boost::filesystem::ofstream adapterTrimmingOstr(adapterTrimmingTxtFile_);
    adapterTrimmingOstr << "Lane\tRead\tProject\tSample Id\tSample Name\tSample Number\tTrimmedBases\tPercentageOfBases" << std::endl;

    std::ostringstream trimLengthOstr;
    stats::ConversionStatsXml conversionStatsXml;
    for (const auto &stats : conversionStats_)
    {
        const layout::LaneInfo& laneInfo = stats.first;

        std::map<std::string, std::string> uniqueSampleNameMap;
        stats::DemuxReportGenerator::makeUniqueSampleNameMap(laneInfo, uniqueSampleNameMap);

        const ConversionStatsForSampleTileBarcode &conversionStats = stats.second;
        unsigned int laneNumber = laneInfo.getNumber();

        if( conversionStats.empty() )
        {
            BCL2FASTQ_LOG(common::LogLevel::WARNING) << "No conversion stats for lane '" << laneNumber << "'" << std::endl;
        }

        stats::JsonObject* jsonReadInfosForLane = NULL;
        jsonReadInfosForLanes->add(std::shared_ptr<stats::JsonType>(jsonReadInfosForLane = new stats::JsonObject));
        jsonReadInfosForLane->add("LaneNumber", std::shared_ptr<stats::JsonType>(new stats::JsonNumber(laneNumber)));
        stats::JsonArray* jsonReadInfos = NULL;
        jsonReadInfosForLane->add("ReadInfos", std::shared_ptr<stats::JsonType>(jsonReadInfos = new stats::JsonArray));

        stats::JsonObject* jsonConversionResultsForLane = NULL;
        jsonConversionResults->add(std::shared_ptr<stats::JsonType>(jsonConversionResultsForLane = new stats::JsonObject));
        jsonConversionResultsForLane->add("LaneNumber", std::shared_ptr<stats::JsonType>(new stats::JsonNumber(laneNumber)));

        stats::JsonNumber* totalClustersRaw = NULL;
        stats::JsonNumber* totalClustersPF = NULL;
        stats::JsonNumber* totalYield = NULL;
        jsonConversionResultsForLane->add("TotalClustersRaw", std::shared_ptr<stats::JsonType>(totalClustersRaw = new stats::JsonNumber(0)));
        jsonConversionResultsForLane->add("TotalClustersPF", std::shared_ptr<stats::JsonType>(totalClustersPF = new stats::JsonNumber(0)));
        jsonConversionResultsForLane->add("Yield", std::shared_ptr<stats::JsonType>(totalYield = new stats::JsonNumber(0)));

        for (const auto& readInfo : laneInfo.getReadInfos())
        {
            stats::JsonObject* jsonReadInfo = NULL;
            jsonReadInfos->add(std::shared_ptr<stats::JsonType>(jsonReadInfo = new stats::JsonObject));
            jsonReadInfo->add("Number",  std::shared_ptr<stats::JsonType>(new stats::JsonNumber(readInfo.getNumber())));
            jsonReadInfo->add("NumCycles",  std::shared_ptr<stats::JsonType>(new stats::JsonNumber(readInfo.cycleInfos().size())));
            jsonReadInfo->add("IsIndexedRead",  std::shared_ptr<stats::JsonType>(new stats::JsonBool(readInfo.isIndexRead())));
        }

        FlowcellProjectSampleBarcodeTileReadBarcodeStats flowcellProjectSampleBarcodeStats;
        ProjectSampleReadBarcodeStats projectSampleBarcodeStats;

        std::string indexName,sampleId,projectName,sampleName;
        common::SampleNumber sampleNumber;

        size_t sampleIndex = 0;
        for (const auto& sampleStats : conversionStats)
        {
            size_t tileIndex = 0;
            for (const auto& tileStats : sampleStats)
            {
                size_t barcodeIndex = 0;
                for (const auto& barcodeStats : tileStats)
                {
                    deserializeMetadata(laneInfo,
                                        sampleIndex,
                                        barcodeIndex,
                                        indexName,
                                        sampleId,
                                        projectName,
                                        sampleName,
                                        sampleNumber);

                    layout::LaneInfo::TileInfosContainer::value_type tileInfo = laneInfo.getTileInfos()[tileIndex];
                    common::TileNumber tileNumber = tileInfo.getNumber();

                    projectSampleBarcodeStats[projectName][sampleId].first = sampleIndex;
                    projectSampleBarcodeStats[projectName][sampleId].second += barcodeStats.first;

                    flowcellProjectSampleBarcodeStats[flowcellId][projectName][sampleId][indexName][tileNumber].first += barcodeStats.first;
                    flowcellProjectSampleBarcodeStats[flowcellId][projectName][sampleId]["all"][tileNumber].first += barcodeStats.first;
                    flowcellProjectSampleBarcodeStats[flowcellId][projectName]["all"]["all"][tileNumber].first += barcodeStats.first;
                    flowcellProjectSampleBarcodeStats[flowcellId]["all"]["all"]["all"][tileNumber].first += barcodeStats.first;

                    for (const auto& readStats : barcodeStats.second)
                    {
                        common::ReadNumber readNumber = readStats.read-1;

                        flowcellProjectSampleBarcodeStats[flowcellId][projectName][sampleId][indexName][tileNumber].second[readNumber] += readStats;
                        flowcellProjectSampleBarcodeStats[flowcellId][projectName][sampleId]["all"][tileNumber].second[readNumber] += readStats;
                        flowcellProjectSampleBarcodeStats[flowcellId][projectName]["all"]["all"][tileNumber].second[readNumber] += readStats;
                        flowcellProjectSampleBarcodeStats[flowcellId]["all"]["all"]["all"][tileNumber].second[readNumber] += readStats;
                    }

                    ++barcodeIndex;
                }

                ++tileIndex;
            }

            ++sampleIndex;
        }

        std::map<std::pair<std::string, std::string>, std::vector<size_t>> trimmedCountMap;
        writeAdapterStats(laneInfo,
                          projectSampleBarcodeStats,
                          adapterTrimmingOstr,
                          trimLengthOstr,
                          trimmedCountMap);

        writeSummaryStats(flowcellProjectSampleBarcodeStats[flowcellId],
                          laneInfo);

        stats::JsonArray* jsonDemuxResults = NULL;
        jsonConversionResultsForLane->add("DemuxResults", std::shared_ptr<stats::JsonArray>(jsonDemuxResults = new stats::JsonArray()));

        BOOST_FOREACH( const FlowcellProjectSampleBarcodeTileReadBarcodeStats::value_type &flowcellStats, flowcellProjectSampleBarcodeStats)
        {
            BOOST_FOREACH( const ProjectSampleBarcodeTileReadBarcodeStats::value_type &projectStats, flowcellStats.second)
            {
                BOOST_FOREACH( const SampleBarcodeTileReadBarcodeStats::value_type &sampleStats, projectStats.second)
                {
                    stats::JsonObject* jsonDemuxSample = NULL;

                    stats::JsonNumber* numReads = NULL;
                    stats::JsonNumber* yield = NULL;
                    int64_t totalYieldValue = 0;

                    std::map<common::ReadNumber, stats::JsonNumber*> readYieldMap;
                    std::map<common::ReadNumber, stats::JsonNumber*> readYieldQ30Map;
                    std::map<common::ReadNumber, stats::JsonNumber*> readQscoreSumMap;

                    stats::JsonArray* jsonReadMetrics = NULL;
                    stats::JsonArray* jsonIndexMetrics = NULL;
                    std::string sampleName;
                    if (sampleStats.first != "all" && projectStats.first != "all")
                    {
                        common::SampleNumber sampleNumber;
                        getSampleNameAndNumber(laneInfo,
                                               projectStats.first,
                                               sampleStats.first,
                                               sampleName,
                                               sampleNumber);

                        if (sampleNumber == 0)
                        {
                            jsonConversionResultsForLane->add("Undetermined", std::shared_ptr<stats::JsonType>(jsonDemuxSample = new stats::JsonObject));
                        }
                        else
                        {
                            jsonDemuxResults->add(std::shared_ptr<stats::JsonType>(jsonDemuxSample = new stats::JsonObject()));

                            jsonDemuxSample->add("SampleId", std::shared_ptr<stats::JsonType>(new stats::JsonString(sampleStats.first)));
                            jsonDemuxSample->add("SampleName", std::shared_ptr<stats::JsonType>(new stats::JsonString(sampleName)));
                        }

                        if (sampleStats.second.size() > 1 &&
                            !(sampleStats.second.begin()->first == layout::Barcode::DEFAULT_BARCODE || (++sampleStats.second.begin())->first == layout::Barcode::DEFAULT_BARCODE))

                        {
                            jsonDemuxSample->add("IndexMetrics", std::shared_ptr<stats::JsonType>(jsonIndexMetrics = new stats::JsonArray));
                        }

                        jsonDemuxSample->add("NumberReads", std::shared_ptr<stats::JsonType>(numReads = new stats::JsonNumber(0)));
                        jsonDemuxSample->add("Yield", std::shared_ptr<stats::JsonType>(yield = new stats::JsonNumber(0)));

                        jsonDemuxSample->add("ReadMetrics", std::shared_ptr<stats::JsonType>(jsonReadMetrics = new stats::JsonArray));
                    }

                    BOOST_FOREACH( const BarcodeTileReadBarcodeStats::value_type &barcodeStats, sampleStats.second)
                    {
                        if (sampleStats.first != "all" &&
                            projectStats.first != "all" && barcodeStats.first != "all" &&
                            jsonIndexMetrics != NULL )
                        {
                            recordBarcodeStats(barcodeLaneStatsMap_[laneInfo.getNumber()],
                                               barcodeStats.first,
                                               *jsonIndexMetrics);
                        }


                        BOOST_FOREACH( const TileReadBarcodeStats::value_type &tileStats, barcodeStats.second)
                        {
                            if (barcodeStats.first == "all")
                            {
                                if (sampleStats.first != "all" && projectStats.first != "all")
                                {
                                    totalClustersRaw->add(tileStats.second.first[stats::TileBarcodeStats::RAW].clustersCount);
                                    totalClustersPF->add(tileStats.second.first[stats::TileBarcodeStats::PF].clustersCount);

                                    numReads->add(tileStats.second.first[stats::TileBarcodeStats::PF].clustersCount);
                                }

                            }

                            conversionStatsXml.addLaneTile( flowcellStats.first,
                                                            projectStats.first,
                                                            sampleStats.first == "all" ? "all" : uniqueSampleNameMap[sampleStats.first],
                                                            barcodeStats.first,
                                                            tileStats.first,
                                                            laneNumber,
                                                            tileStats.second.first );

                            recordReadStats(conversionStatsXml,
                                            tileStats.second.second,
                                            laneNumber,
                                            tileStats.first,
                                            barcodeStats.first,
                                            projectStats.first,
                                            sampleStats.first,
                                            sampleName,
                                            uniqueSampleNameMap,
                                            readYieldMap,
                                            readYieldQ30Map,
                                            readQscoreSumMap,
                                            trimmedCountMap,
                                            yield,
                                            jsonReadMetrics);

                            if (sampleStats.first != "all" && projectStats.first != "all")
                            {
                                totalYieldValue = yield->getValue();
                            }
                        }
                    }
                    totalYield->add(totalYieldValue);
                }
            }
        }
    }

    recordUnknownBarcodes(jsonStats,
                          conversionStatsXml);

    std::ofstream os(conversionStatsXmlFile_.string().c_str());
    if (!os) {
        BOOST_THROW_EXCEPTION(common::IoError(errno, "Unable to open file for writing: " + conversionStatsXmlFile_.string()));
    }
    if (!(os << conversionStatsXml)) {
        BOOST_THROW_EXCEPTION(common::IoError(errno, "Failed to store conversion statistics in : " + conversionStatsXmlFile_.string()));
    }

    adapterTrimmingOstr << trimLengthOstr.str();

    std::ofstream jsonStatsStr(statsJsonFile_.string().c_str());
    jsonStats.stream(jsonStatsStr);
}

void ConverterStats::recordBarcodeStats(const BarcodeStatsMap& barcodeStatsMap,
                                        const std::string& barcodeName,
                                        stats::JsonArray& jsonIndexMetrics)
{
    std::vector<stats::JsonNumber*> mismatchCountVec;
    stats::JsonObject* mismatchCounts = NULL;

    stats::JsonObject* indexMetric = NULL;
    jsonIndexMetrics.add(std::shared_ptr<stats::JsonType>(indexMetric = new stats::JsonObject));
    indexMetric->add("IndexSequence", std::shared_ptr<stats::JsonType>(new stats::JsonString(barcodeName)));
    indexMetric->add("MismatchCounts", std::shared_ptr<stats::JsonType>(mismatchCounts = new stats::JsonObject));

    size_t mismatches = 0;
    for (const common::TotalClustersCount count : barcodeStatsMap.at(barcodeName).getBarcodeMismatchCounts())
    {
        if (count == 0 && mismatches > 1)
        {
            break;
        }

        if (mismatchCountVec.size() <= mismatches)
        {
            stats::JsonNumber* mismatchCount = NULL;
            mismatchCounts->add(boost::lexical_cast<std::string>(mismatches), std::shared_ptr<stats::JsonType>(mismatchCount = new stats::JsonNumber(count)));
            mismatchCountVec.push_back(mismatchCount);
        }
        else
        {
            mismatchCountVec[mismatches]->add(count);
        }

        ++mismatches;
    }
}

void ConverterStats::recordReadStats(stats::ConversionStatsXml& conversionStatsXml,
                                     const AnnotatedTileReadBarcodeStats& tileStats,
                                     unsigned int laneNumber,
                                     common::TileNumber tileNumber,
                                     const std::string& barcodeName,
                                     const std::string& project,
                                     const std::string& sampleId,
                                     const std::string& sampleName,
                                     std::map<std::string, std::string>& uniqueSampleNameMap,
                                     std::map<common::ReadNumber, stats::JsonNumber*>& readYieldMap,
                                     std::map<common::ReadNumber, stats::JsonNumber*>& readYieldQ30Map,
                                     std::map<common::ReadNumber, stats::JsonNumber*>& readQscoreSumMap,
                                     std::map<std::pair<std::string, std::string>, std::vector<size_t>>& trimmedCountMap,
                                     stats::JsonNumber* yield,
                                     stats::JsonArray* jsonReadMetrics)
{
    for (const auto& readStats : tileStats)
    {
        if (barcodeName == "all" &&
            sampleId != "all" &&
            project != "all")
        {
            yield->add(readStats.second[stats::TileBarcodeStats::PF].yield);

            stats::JsonNumber*& readYield = readYieldMap[readStats.first];
            stats::JsonNumber*& readYieldQ30 = readYieldQ30Map[readStats.first];
            stats::JsonNumber*& readQscoreSum = readQscoreSumMap[readStats.first];
            if (!readYield)
            {
                stats::JsonObject* readMetric = NULL;
                jsonReadMetrics->add(std::shared_ptr<stats::JsonType>(readMetric = new stats::JsonObject));
                readMetric->add("ReadNumber", std::shared_ptr<stats::JsonType>(new stats::JsonNumber(readStats.first+1)));
                readMetric->add("Yield", std::shared_ptr<stats::JsonType>(readYield = new stats::JsonNumber(0)));
                readMetric->add("YieldQ30", std::shared_ptr<stats::JsonType>(readYieldQ30 = new stats::JsonNumber(0)));
                readMetric->add("QualityScoreSum", std::shared_ptr<stats::JsonType>(readQscoreSum = new stats::JsonNumber(0)));
                readMetric->add("TrimmedBases", std::shared_ptr<stats::JsonType>(new stats::JsonNumber(trimmedCountMap[std::make_pair(sampleId, sampleName)][readStats.first])));
            }

            readYield->add(readStats.second[stats::TileBarcodeStats::PF].yield);
            readYieldQ30->add(readStats.second[stats::TileBarcodeStats::PF].yieldQ30);
            readQscoreSum->add(readStats.second[stats::TileBarcodeStats::PF].qualityScoreSum);

        }

        conversionStatsXml.addLaneTileRead( layout_.getFlowcellInfo().getFlowcellId(),
                                            project,
                                            sampleId == "all" ? "all" : uniqueSampleNameMap[sampleId],
                                            barcodeName,
                                            tileNumber,
                                            readStats.first + 1,
                                            laneNumber,
                                            readStats.second );
    }
}

void ConverterStats::recordUnknownBarcodes(stats::JsonObject& jsonStats,
                                           stats::ConversionStatsXml& conversionStatsXml) const
{
    stats::JsonArray* jsonUnknownBarcodes = NULL;
    jsonStats.add("UnknownBarcodes", std::shared_ptr<stats::JsonType>(jsonUnknownBarcodes = new stats::JsonArray));

    BOOST_FOREACH (const UnknownPair &stats, std::make_pair(topUnknownBarcodes_.begin(), topUnknownBarcodes_.end()))
    {
        stats::JsonObject* jsonUnknownBarcodesForLane = NULL;
        jsonUnknownBarcodes->add(std::shared_ptr<stats::JsonType>(jsonUnknownBarcodesForLane = new stats::JsonObject));

        layout::LaneInfo laneInfo = stats.first;
        unsigned int laneNumber = laneInfo.getNumber();
        const stats::BarcodeHits::Popular& popularBarcodes = stats.second.getPopularBarcodes();

        jsonUnknownBarcodesForLane->add("Lane", std::shared_ptr<stats::JsonType>(new stats::JsonNumber(laneNumber)));

        stats::JsonObject* jsonBarcodes = NULL;
        jsonUnknownBarcodesForLane->add("Barcodes", std::shared_ptr<stats::JsonType>(jsonBarcodes = new stats::JsonObject));
        for (const auto& barcodeHit : popularBarcodes)
        {
            jsonBarcodes->add(barcodeHit.first.toString(), std::shared_ptr<stats::JsonType>(new stats::JsonNumber(barcodeHit.second)));
        }

        conversionStatsXml.addFlowcellLane( layout_.getFlowcellInfo().getFlowcellId(),
                                            laneNumber,
                                            popularBarcodes.size() <= 10 ?
                                                popularBarcodes :
                                                stats::BarcodeHits::Popular(popularBarcodes.begin(), popularBarcodes.begin()+10) );
    }
}

namespace detail
{

FastqSummaryStats::FastqSummaryStats(common::SampleNumber  sampleNumber,
                                     const std::string&    sampleName,
                                     common::TileNumber    tileNumber,
                                     common::TotalClustersCount numberReadsRaw,
                                     common::TotalClustersCount numberReadsPF)
: sampleNumber_(sampleNumber)
, sampleName_(sampleName == "Undetermined" ? "None" : sampleName)
, tileNumber_(tileNumber)
, numberReadsRaw_(numberReadsRaw)
, numberReadsPF_(numberReadsPF)
{
}

bool FastqSummaryStats::operator<(const FastqSummaryStats& other) const
{
    if (tileNumber_ < other.tileNumber_) return true;
    if (tileNumber_ > other.tileNumber_) return false;
    return sampleNumber_ < other.sampleNumber_;
}

std::ostream& operator<<(std::ostream& os, const FastqSummaryStats& stats)
{
    os << stats.sampleNumber_ << "\t"
       << stats.tileNumber_ << "\t"
       << stats.numberReadsRaw_ << "\t"
       << stats.numberReadsPF_;

    return os;
}

}

void ConverterStats::writeAdapterStats(const layout::LaneInfo&              laneInfo,
                                       const ProjectSampleReadBarcodeStats& stats,
                                       std::ostream&                        adapterOstr,
                                       std::ostream&                        trimLengthOstr,
                                       std::map<std::pair<std::string, std::string>, std::vector<size_t>>& trimmedBaseMap) const
{
    BOOST_FOREACH(const ProjectSampleReadBarcodeStats::value_type& sampleBarcodeTileStats, stats)
    {
        BOOST_FOREACH(const SampleReadBarcodeStats::value_type& barcodeTileStats, sampleBarcodeTileStats.second)
        {
            const layout::SampleInfo& sampleInfo = laneInfo.getSampleInfos().at(barcodeTileStats.second.first);

            writeAdapterStats(laneInfo.getNumber(),
                              sampleBarcodeTileStats.first,
                              barcodeTileStats.first,
                              sampleInfo.getName(),
                              sampleInfo.getNumber(),
                              barcodeTileStats.second.second[stats::TileBarcodeStats::PF].trimLengthCountForReads,
                              adapterOstr,
                              trimLengthOstr,
                              trimmedBaseMap);
        }
    }
}

void ConverterStats::writeSummaryStats(const ProjectSampleBarcodeTileReadBarcodeStats& stats,
                                       const layout::LaneInfo&                         laneInfo) const
{
    detail::FastqSummaryStatsContainer fastqSummaryStats;

    for (const ProjectSampleBarcodeTileReadBarcodeStats::value_type& sampleBarcodeTileStats : stats)
    {
        if (sampleBarcodeTileStats.first == "all") continue;

        for (const SampleBarcodeTileReadBarcodeStats::value_type& barcodeTileStats : sampleBarcodeTileStats.second)
        {
            if (barcodeTileStats.first == "all") continue;

            BarcodeTileReadBarcodeStats::const_iterator pos = barcodeTileStats.second.find("all");
            BOOST_ASSERT(pos != barcodeTileStats.second.end());

            for (const TileReadBarcodeStats::value_type& tileStats : pos->second)
            {
                std::string sampleName;
                common::SampleNumber sampleNumber;
                getSampleNameAndNumber(laneInfo,
                                       sampleBarcodeTileStats.first,
                                       barcodeTileStats.first,
                                       sampleName,
                                       sampleNumber);

                fastqSummaryStats.push_back(detail::FastqSummaryStats(
                    sampleNumber,
                    sampleName,
                    tileStats.first, // tileNumber
                    tileStats.second.first[stats::TileBarcodeStats::RAW].clustersCount,
                    tileStats.second.first[stats::TileBarcodeStats::PF].clustersCount));
            }
        }
    }

    std::sort(fastqSummaryStats.begin(),
              fastqSummaryStats.end());

    writeFastqSummaryStats(laneInfo.getNumber(), fastqSummaryStats);
    writeDemuxSummaryStats(laneInfo.getNumber(), fastqSummaryStats);
}

void ConverterStats::writeFastqSummaryStats(common::LaneNumber                        laneNumber,
                                            const detail::FastqSummaryStatsContainer& fastqSummaryStats) const
{
    boost::filesystem::path fastqSummaryFileName(fastqSummaryPartialFileName_.string() + boost::lexical_cast<std::string>(laneNumber) + ".txt");
    boost::filesystem::ofstream fastqSummaryOstr(fastqSummaryFileName);
    fastqSummaryOstr << "SampleNumber\tTile\tNumberOfReadsRaw\tNumberOfReadsPF" << std::endl;

    BOOST_FOREACH(const detail::FastqSummaryStats& summaryStats, fastqSummaryStats)
    {
        fastqSummaryOstr << summaryStats << std::endl;
    }
}

void ConverterStats::writeDemuxSummaryStats(common::LaneNumber                        laneNumber,
                                            const detail::FastqSummaryStatsContainer& fastqSummaryStats) const
{
    if (fastqSummaryStats.empty())
    {
        return;
    }

    boost::filesystem::path demuxSummaryFileName(demuxSummaryPartialFileName_.string() + boost::lexical_cast<std::string>(laneNumber) + ".txt");
    boost::filesystem::ofstream demuxSummaryOstr(demuxSummaryFileName);

    std::stringstream sampleNamesStr;
    sampleNamesStr << "SampleName";
    demuxSummaryOstr << "SampleNumber";

    size_t totalReads = 0;
    bool firstLoop = true;
    std::vector<common::ClustersCount> rawReadsPerTile;
    BOOST_FOREACH(const detail::FastqSummaryStats& stats, fastqSummaryStats)
    {
        if (stats.sampleNumber_ == fastqSummaryStats.begin()->sampleNumber_)
        {
            if (!firstLoop)
            {
                rawReadsPerTile.push_back(totalReads);
            }
            else
            {
                firstLoop = false;
            }
            totalReads = 0;
        }

        if (rawReadsPerTile.empty())
        {
            demuxSummaryOstr << "\t" << stats.sampleNumber_;
            sampleNamesStr << "\t" << stats.sampleName_;
        }
        totalReads += stats.numberReadsPF_;
    }
    rawReadsPerTile.push_back(totalReads);
    demuxSummaryOstr << std::endl << sampleNamesStr.str();

    int tileIndex = -1;
    BOOST_FOREACH(const detail::FastqSummaryStats& stats, fastqSummaryStats)
    {
        if (stats.sampleNumber_ == fastqSummaryStats.begin()->sampleNumber_)
        {
            ++tileIndex;
            demuxSummaryOstr << std::endl << "L" << laneNumber << "T" << stats.tileNumber_;
        }

        demuxSummaryOstr
            << "\t" << std::setprecision(7)
            << (rawReadsPerTile.at(tileIndex) == 0 ? 0 : 100.0 * stats.numberReadsPF_ / rawReadsPerTile.at(tileIndex));
    }

    demuxSummaryOstr
        << std::endl << std::endl
        << "### Most Popular Unknown Index Sequences" << std::endl
        << "### Columns: Index_Sequence Hit_Count" << std::endl;

    BOOST_FOREACH (const UnknownPair &stats, std::make_pair(topUnknownBarcodes_.begin(), topUnknownBarcodes_.end()))
    {
        layout::LaneInfo laneInfo = stats.first;
        if (laneNumber == laneInfo.getNumber())
        {
            BOOST_FOREACH(const stats::BarcodeHits::Popular::value_type& barcodeHits, stats.second.getPopularBarcodes())
            {
                demuxSummaryOstr << barcodeHits.first << "\t" << barcodeHits.second << std::endl;
            }
            break;
        }
    }
}

// This function is inefficient, but it's not important for the overall runtime.
void ConverterStats::getSampleNameAndNumber(const layout::LaneInfo& laneInfo,
                                            const std::string&      project,
                                            const std::string&      sampleId,
                                            std::string&            sampleName,
                                            common::SampleNumber&   sampleNumber) const
{
    BOOST_FOREACH(const layout::SampleInfo& sampleInfo, std::make_pair(laneInfo.sampleInfosBegin(), laneInfo.sampleInfosEnd()))
    {
        if (project == sampleInfo.getProject() &&
            sampleId == sampleInfo.getId())
        {
            sampleName = sampleInfo.getName();
            sampleNumber = sampleInfo.getNumber();
            return;
        }
    }

    BCL2FASTQ_LOG(common::LogLevel::WARNING)
        << (boost::format("Sample number was not found for project: '%s' and sample id: '%s'") % project % sampleId).str()
        << std::endl;
}

void ConverterStats::writeAdapterStats(common::LaneNumber                                   laneNumber,
                                       const std::string&                                   project,
                                       const std::string&                                   sampleId,
                                       const std::string&                                   sampleName,
                                       common::SampleNumber                                 sampleNumber,
                                       const stats::AllReadsStats::TrimLengthCountForReads& trimLengthCountForReads,
                                       std::ostream&                                        trimBasesOstr,
                                       std::ostream&                                        trimLengthOstr,
                                       std::map<std::pair<std::string, std::string>, std::vector<size_t>>& trimmedBaseMap) const
{
    common::ReadNumber readNumber = 0;
    BOOST_FOREACH(const stats::AllReadsStats::TrimLengthCount& trimLengthCount, trimLengthCountForReads)
    {
        size_t totalTrimBasesForRead = 0;
        size_t trimLength = 0;
        size_t basesPerRead = trimLengthCount.size() - 1;
        size_t numReads = 0;

        BOOST_FOREACH(size_t count, trimLengthCount)
        {
            numReads += count;
            totalTrimBasesForRead += trimLength++ * count;
        }

        trimmedBaseMap[std::make_pair(sampleId, sampleName)].push_back(totalTrimBasesForRead);

        trimLengthOstr
            << std::endl
            << "Lane: " << laneNumber
            << "\tRead: " << readNumber+1
            << "\tProject: " << project
            << "\tSample Id: " << sampleId
            << "\tSample Name: " << sampleName
            << "\tSample Number: " << sampleNumber
            << std::endl;

        trimLength = 0;
        for (size_t count : boost::adaptors::reverse(trimLengthCount))
        {
            trimLengthOstr << trimLength++ << "\t" << static_cast<double>(count) / (numReads == 0 ? 1 : numReads) << std::endl;
        }

        trimBasesOstr
           << laneNumber << "\t"
           << ++readNumber << "\t"
           << project << "\t"
           << sampleId << "\t"
           << sampleName << "\t"
           << sampleNumber << "\t"
           << totalTrimBasesForRead << "\t"
           << (100.0*totalTrimBasesForRead) / (numReads == 0 ? 1 : numReads*basesPerRead)
           << std::endl;
    }
}

} // namespace conversion
} // namespace bcl2fastq


