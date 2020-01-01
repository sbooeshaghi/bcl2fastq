/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file FastqCreator.hpp
 *
 * \brief Declaration of FastqCreator.
 *
 * \author Aaron Day
 */


#ifndef BCL2FASTQ_CONVERSION_FASTQ_CREATOR_HPP
#define BCL2FASTQ_CONVERSION_FASTQ_CREATOR_HPP


#include "conversion/FastqCreator.hh"

#include "common/FastIo.hh"
#include "io/GzipCompressor.hh"
#include "conversion/BclBaseConversion.hh"

#include <boost/foreach.hpp>

namespace bcl2fastq {
namespace conversion {

template<common::NumBasesPerByte numBasesPerByte>
FastqCreateTask<numBasesPerByte>::FastqCreateTask(
    std::shared_ptr<FastqCreateTaskManager>& taskManager,
    const data::BclBufferVec &bclBuffers,
    size_t cycleIndex,
    size_t cycleIndexEnd,
    const layout::FlowcellInfo &flowcellInfo,
    const layout::LaneInfo &laneInfo,
    const layout::ReadInfo& currentReadInfo,
    SampleIndex::FastqOffsetsContainer::const_iterator offsetsBegin,
    SampleIndex::FastqOffsetsContainer::const_iterator offsetsEnd,
    std::size_t maskShortAdapterReads,
    float adapterStringency,
    const boost::ptr_vector<AdapterLocator<numBasesPerByte>>& maskAdapters,
    const boost::ptr_vector<AdapterLocator<numBasesPerByte>>& trimAdapters,
    bool generateReverseComplementFastqs,
    bool includeNonPfClusters,
    bool useBgzf,
    int compressionLevel,
    bool findAdaptersWithSlidingWindow,
    FastqCreateStats& stats,
    FastqBuffer::FastqsContainer::value_type::value_type::value_type &outputBuffer
)
: Task(taskManager)
, bclBuffers_(bclBuffers)
, cycleIndex_(cycleIndex + currentReadInfo.getBclBufferOffset())
, cycleIndexEnd_(cycleIndex_ + currentReadInfo.cycleInfos().size())
, flowcellInfo_(flowcellInfo)
, laneInfo_(laneInfo)
, readInfo_(currentReadInfo)
, sampleMetadata_(bclBuffers_.at(offsetsBegin->first).samples_.at(offsetsBegin->second).first)
, offsetsBegin_(offsetsBegin)
, offsetsEnd_(offsetsEnd)
, outputBuffer_(outputBuffer)
, buffer_()
, useBgzf_(useBgzf)
, compressionLevel_(compressionLevel)
, findAdaptersWithSlidingWindow_(findAdaptersWithSlidingWindow)
, minimumTrimmedReadLength_(laneInfo_.getMinimumTrimmedReadLength())
, maskShortAdapterReads_(maskShortAdapterReads)
, adapterStringency_(adapterStringency)
, maskAdapters_(maskAdapters)
, trimAdapters_(trimAdapters)
, generateReverseComplementFastqs_(generateReverseComplementFastqs)
, includeNonPfClusters_(includeNonPfClusters)
, stats_(stats)
, umiCycles_()
{
    common::CycleNumber umiCycleOffset = 0;
    for (const auto& readInfo : laneInfo_.readInfos())
    {
        // We do not load BCL files for cycles we don't use.
        // The index into the BCL buffer needs to account for this.

        const common::CycleRange& umiCycles = readInfo.getUmiCycles();
        if (umiCycles.second != 0)
        {
            // These are the indexes into the bcl buffer
            umiCycles_.push_back(common::CycleRange(umiCycles.first + umiCycleOffset, umiCycles.second + umiCycleOffset));
        }

        umiCycleOffset += readInfo.cyclesToLoad().size();
    }
}

template<common::NumBasesPerByte numBasesPerByte>
bool FastqCreateTask<numBasesPerByte>::execute(int32_t threadNum)
{
    io::GzipCompressor compressor(outputBuffer_,
                                  useBgzf_,
                                  boost::iostreams::gzip_params(compressionLevel_, // config param with default=4
                                                                boost::iostreams::zlib::deflated, // default
                                                                15, // default
                                                                9));

    for (auto offset = offsetsBegin_; offset != offsetsEnd_; ++offset)
    {
        fastqCreate(*offset,
                    threadNum,
                    compressor);
    }

    return true;
}


namespace detail {


static const char fastqHeaderStart = '@';
static const char fastqHeaderSeparator = '\n';
static const char fastqHeaderValuesSeparator = ':';
static const char fastqHeaderPositionsSeparator = ' ';
static const char fastqHeaderUmiSeparator = '+';
static const std::string fastqBasesSeparator("\n+\n");
static const char fastqQualitiesSeparator('\n');
static const char fastqMaskedBase('N');
static const char fastqMaskedQuality('#');


} // namespace detail

template<common::NumBasesPerByte numBasesPerByte>
bool FastqCreateTask<numBasesPerByte>::computeStatistics( const FastqConstIterator<numBasesPerByte>& basesBegin,
                                                          const FastqConstIterator<numBasesPerByte>& basesEnd,
                                                          SampleIndex::FastqOffsetsContainer::value_type offset,
                                                          int32_t threadNum,
                                                          bool filterFlag,
                                                          size_t trimmedCount,
                                                          const std::vector<uint32_t>& remappedQscores )
{
    layout::LaneInfo::TileInfosContainer::size_type tileIndex =  std::distance( laneInfo_.getTileInfos().begin(),
                                                                                bclBuffers_[offset.first].tileInfo_ );
    stats::ReadStats readStats;
    readStats.yield = std::distance(basesBegin, basesEnd);
    readStats.yieldQ30 = remappedQscores.empty() ?
        std::count_if(basesBegin,
                      basesEnd,
                      [](unsigned char bcl)
                          { const unsigned char quality = bcl >> 2;
                            return quality >= 30; }) :
        std::count_if(basesBegin,
                      basesEnd,
                      [&remappedQscores](unsigned char bcl)
                          { const unsigned char quality = remappedQscores[bcl>>2];
                            return quality >= 30; });

    readStats.qualityScoreSum = remappedQscores.empty() ?
        std::accumulate(basesBegin,
                        basesEnd,
                        common::TotalClustersCount(0),
                        [](common::TotalClustersCount total, unsigned char bcl)
                            { return total + (bcl>>2); }) :
        std::accumulate(basesBegin,
                        basesEnd,
                        common::TotalClustersCount(0),
                        [&remappedQscores](common::TotalClustersCount total, unsigned char bcl)
                            { return total + remappedQscores[bcl>>2]; });

    auto& statsForBarcode = stats_.at(threadNum)->conversionStats_.at(sampleMetadata_.sampleIndex_).at(tileIndex).at(sampleMetadata_.barcodeIndex_);
    statsForBarcode.first[stats::TileBarcodeStats::RAW].RecordAdapterStats(trimmedCount, readInfo_.getNumber());
    statsForBarcode.second.at(readInfo_.getNumber() - 1)[stats::TileBarcodeStats::RAW] += readStats;
    if( filterFlag )
    {
        if( 1 == readInfo_.getNumber() )  { statsForBarcode.first[stats::TileBarcodeStats::PF].clustersCount += 1; }
        statsForBarcode.first[stats::TileBarcodeStats::PF].RecordAdapterStats(trimmedCount, readInfo_.getNumber());
        statsForBarcode.second.at(readInfo_.getNumber() - 1)[stats::TileBarcodeStats::PF] += readStats;
    }

    return (0 != sampleMetadata_.sampleIndex_);
}

template<common::NumBasesPerByte numBasesPerByte>
void FastqCreateTask<numBasesPerByte>::fastqCreate(
    SampleIndex::FastqOffsetsContainer::value_type offset,
    int32_t threadNum,
    io::GzipCompressor& compressor
)
{
    const data::BclBuffer &bclBuffer = bclBuffers_.at(offset.first);
    layout::LaneInfo::TileInfosContainer::size_type tileIndex =  std::distance( laneInfo_.getTileInfos().begin(),
                                                                                bclBuffer.tileInfo_ );
    auto& statsForBarcode = stats_.at(threadNum)->conversionStats_.at(sampleMetadata_.sampleIndex_).at(tileIndex).at(sampleMetadata_.barcodeIndex_);
    if( 1 == readInfo_.getNumber() )  { statsForBarcode.first[stats::TileBarcodeStats::RAW].clustersCount += 1; }

    const bool filterFlag = bclBuffer.filters_[offset.second].data_;
    const FastqConstIterator<numBasesPerByte> begin = FastqConstIterator<numBasesPerByte>(bclBuffer, offset.second) + cycleIndex_;
    const FastqConstIterator<numBasesPerByte> end = FastqConstIterator<numBasesPerByte>(bclBuffer, offset.second) + cycleIndexEnd_;

    size_t trimmedCount = 0;
    if (includeNonPfClusters_ || filterFlag)
    {
        buffer_.clear();
        createHeader(offset,
                     filterFlag);

        createBasesAndQualities(offset,
                                compressor,
                                trimmedCount);
    }

    if (!filterFlag && !bclBuffer.bcls_.back().includeNonPf_)
    {
        return;
    }

    if (readInfo_.isDataRead())
    {
        if (!computeStatistics(begin, end, offset, threadNum, filterFlag, trimmedCount, bclBuffer.bcls_[0].remappedQscores_))
        {
            if( 1 == readInfo_.getNumber() && (includeNonPfClusters_ || filterFlag) )  { stats_.at(threadNum)->unknownBarcodes_.record( bclBuffer.samples_[offset.second].second ); }
        }
    }

}

template<common::NumBasesPerByte numBasesPerByte>
void FastqCreateTask<numBasesPerByte>::writeHeaderElement(const std::string& element)
{
    size_t bufferSize = buffer_.size();
    buffer_.resize(bufferSize + element.size() + 1);
    std::copy(element.begin(), element.end(), buffer_.begin() + bufferSize);
    buffer_.back() = detail::fastqHeaderValuesSeparator;
}

template<common::NumBasesPerByte numBasesPerByte>
void FastqCreateTask<numBasesPerByte>::createHeader(SampleIndex::FastqOffsetsContainer::value_type offset,
                                                    bool                                           filterFlag)
{
    const auto& bclBuffer = bclBuffers_.at(offset.first);
    const common::ControlFlag controlFlag = bclBuffer.controls_.at(offset.second).data_;

    // HEADER
    buffer_.push_back(detail::fastqHeaderStart);
    // instrument
    writeHeaderElement(flowcellInfo_.getInstrument());
    // run number
    writeHeaderElement(flowcellInfo_.getRunNumber());
    // flow cell ID
    writeHeaderElement(flowcellInfo_.getFlowcellId());

    // lane number
    BOOST_ASSERT(laneInfo_.getNumber() < 1000);

    if (laneInfo_.getNumber() >= 100)
    {
        buffer_.push_back(char((laneInfo_.getNumber() / 100) + '0'));
        int32_t remainder = laneInfo_.getNumber() % 100;
        buffer_.push_back(char((remainder / 10) + '0'));
        buffer_.push_back(char((remainder % 10) + '0'));
    }
    else if (laneInfo_.getNumber() >= 10)
    {
        buffer_.push_back(char((laneInfo_.getNumber() / 10) + '0'));
        buffer_.push_back(char((laneInfo_.getNumber() % 10) + '0'));
    }
    else
    {
        buffer_.push_back(char(laneInfo_.getNumber() + '0'));
    }
    buffer_.push_back(detail::fastqHeaderValuesSeparator);
    // tile number
    common::putUnsignedInteger(bclBuffer.tileInfo_->getNumber(), buffer_);
    buffer_.push_back(detail::fastqHeaderValuesSeparator);
    // x position
    common::putUnsignedInteger( bclBuffer.positions_->size() > offset.second
                              ? (*bclBuffer.positions_)[offset.second].x_
                              : 0, buffer_ );
    buffer_.push_back(detail::fastqHeaderValuesSeparator);
    // y position
    common::putUnsignedInteger( bclBuffer.positions_->size() > offset.second
                              ? (*bclBuffer.positions_)[offset.second].y_
                              : 0, buffer_ );
    // UMI
    if (!umiCycles_.empty())
    {
        buffer_.push_back(detail::fastqHeaderValuesSeparator);
        size_t i = 0;
        for (const auto& cycle : umiCycles_)
        {
            if (i > 0)
            {
                buffer_.push_back(detail::fastqHeaderUmiSeparator);
            }
            std::transform(
                FastqConstIterator<numBasesPerByte>(bclBuffer, offset.second) + cycle.first,
                FastqConstIterator<numBasesPerByte>(bclBuffer, offset.second) + cycle.second,
                std::back_inserter(buffer_),
                bclBuffer.bcls_[0].numBitsPerQscore_ == 0 ? convertBcl2FastqBaseNoQscore : &convertBcl2FastqBase
            );

            ++i;
        }
    }
    buffer_.push_back(detail::fastqHeaderPositionsSeparator);

    // read number
    BOOST_ASSERT(readInfo_.getNumber() <= 9);
    buffer_.push_back(char(readInfo_.getNumber() + '0'));
    buffer_.push_back(detail::fastqHeaderValuesSeparator);
    // is-filtered
    buffer_.push_back(filterFlag ? 'N' : 'Y');
    buffer_.push_back(detail::fastqHeaderValuesSeparator);
    // control number
    common::putUnsignedInteger(controlFlag, buffer_);
    buffer_.push_back(detail::fastqHeaderValuesSeparator);

    // barcode
    if (bclBuffer.samples_[offset.second].second.componentsBegin() == bclBuffer.samples_[offset.second].second.componentsEnd())
    {
        common::putUnsignedInteger((laneInfo_.sampleInfosBegin()+sampleMetadata_.sampleIndex_)->getNumber(), buffer_);
    }
    else
    {
        bclBuffer.samples_[offset.second].second.output(buffer_);
    }
    buffer_.push_back(detail::fastqHeaderSeparator);
}

template<common::NumBasesPerByte numBasesPerByte>
void FastqCreateTask<numBasesPerByte>::createBasesAndQualities(SampleIndex::FastqOffsetsContainer::value_type offset,
                                                               io::GzipCompressor& compressor,
                                                               size_t& trimmedCount)
{
    const auto& bclBuffer = bclBuffers_.at(offset.first);
    const FastqConstIterator<numBasesPerByte> begin = FastqConstIterator<numBasesPerByte>(bclBuffer, offset.second) + cycleIndex_;

    // The last cycle we want to put in the FASTQ file
    const FastqConstIterator<numBasesPerByte> end = FastqConstIterator<numBasesPerByte>(bclBuffer, offset.second) + cycleIndexEnd_;

    // BASES
    FastqConstIterator<numBasesPerByte> maskBegin(end);
    for (const AdapterLocator<numBasesPerByte>& adapterLocator : maskAdapters_)
    {
        const FastqConstIterator<numBasesPerByte> adapterPos = adapterLocator.identifyAdapter(begin, end);
        maskBegin = std::min(maskBegin, adapterPos);
    }

    FastqConstIterator<numBasesPerByte> trimBegin(end);
    for (const AdapterLocator<numBasesPerByte>& adapterLocator : trimAdapters_)
    {
        const FastqConstIterator<numBasesPerByte> adapterPos = adapterLocator.identifyAdapter(begin, end);
        trimBegin = std::min(trimBegin, adapterPos);
    }

    trimmedCount = end - std::min(trimBegin, maskBegin);
    if (readInfo_.isDataRead())
    {
        if ((trimBegin - begin) < typename FastqConstIterator<numBasesPerByte>::difference_type(minimumTrimmedReadLength_))
        {
            maskBegin = std::min(trimBegin, maskBegin);
            trimBegin = begin + minimumTrimmedReadLength_;
        }

        if ((maskBegin - begin) < typename FastqConstIterator<numBasesPerByte>::difference_type(maskShortAdapterReads_))
        {
            maskBegin = begin;
        }
    }

    const typename FastqConstIterator<numBasesPerByte>::difference_type maskedCount = trimBegin - maskBegin;
    FastqConstIterator<numBasesPerByte> basesEnd(std::min(maskBegin, trimBegin));
    size_t bufferSize = buffer_.size();

    size_t numBases = std::distance(begin, basesEnd);
    size_t numBasesPlusMask = numBases;
    if (maskedCount > 0)
    {
        numBasesPlusMask += maskedCount;
    }

    // Resize for bases + BasesSeparator + quality + QualitySeparator
    buffer_.resize(bufferSize + 2*numBasesPlusMask + 1 + detail::fastqBasesSeparator.size());
    std::vector<char>::iterator iter(buffer_.begin());
    std::advance(iter, bufferSize);
    if (generateReverseComplementFastqs_)
    {
        if (maskedCount > 0)
        {
            std::fill_n(iter, maskedCount, detail::fastqMaskedBase);
            std::advance(iter, maskedCount);
        }

        std::transform(
            FastqConstReverseIterator<numBasesPerByte>(basesEnd),
            FastqConstReverseIterator<numBasesPerByte>(begin),
            iter,
            bclBuffer.bcls_[0].numBitsPerQscore_ == 0 ? &convertBcl2FastqBaseComplementNoQscore : &convertBcl2FastqBaseComplement
        );

        std::advance(iter, numBases);
    }
    else
    {
        std::transform(
            begin,
            basesEnd,
            iter,
            bclBuffer.bcls_[0].numBitsPerQscore_ == 0 ? &convertBcl2FastqBaseNoQscore : &convertBcl2FastqBase
        );

        std::advance(iter, numBases);
        if (maskedCount > 0)
        {
            std::fill_n(iter, maskedCount, detail::fastqMaskedBase);
            std::advance(iter, maskedCount);
        }
    }

    std::copy(detail::fastqBasesSeparator.begin(), detail::fastqBasesSeparator.end(), iter);
    std::advance(iter, detail::fastqBasesSeparator.size());

    // QUALITIES
    if (generateReverseComplementFastqs_)
    {
        if (maskedCount > 0)
        {
            std::fill_n(iter, maskedCount, detail::fastqMaskedQuality);
            std::advance(iter, maskedCount);
        }

        if ( bclBuffers_[offset.first].bcls_[0].remappedQscores_.empty())
        {
            std::transform(FastqConstReverseIterator<numBasesPerByte>(basesEnd),
                           FastqConstReverseIterator<numBasesPerByte>(begin),
                           iter,
                           &convertBcl2FastqQuality);
        }
        else
        {
            std::transform(FastqConstReverseIterator<numBasesPerByte>(basesEnd),
                           FastqConstReverseIterator<numBasesPerByte>(begin),
                           iter,
                           [this, &offset] (unsigned char bcl) { return convertBcl2FastqQuality2(bcl, bclBuffers_[offset.first].bcls_[0].remappedQscores_); });
        }

        std::advance(iter, numBases);
    }
    else
    {
        if (bclBuffers_[offset.first].bcls_[0].remappedQscores_.empty())
        {
            std::transform(begin, basesEnd, iter, &convertBcl2FastqQuality);
        }
        else
        {
            std::transform(begin,
                           basesEnd,
                           iter,
                           [this, &offset] (unsigned char bcl) { return convertBcl2FastqQuality2(bcl, bclBuffers_[offset.first].bcls_[0].remappedQscores_); });
        }

        std::advance(iter, numBases);
        if (maskedCount > 0)
        {
            std::fill_n(iter, maskedCount, detail::fastqMaskedQuality);
            std::advance(iter, maskedCount);
        }
    }
    *iter = detail::fastqQualitiesSeparator;

    std::streamsize bytesWritten = compressor.write(&*buffer_.begin(), buffer_.size());
    BCL2FASTQ_ASSERT_MSG( static_cast<FastqBuffer::FastqsContainer::value_type::value_type::value_type::size_type>(bytesWritten) == buffer_.size(),
                          "Only " << bytesWritten << " of " << buffer_.size() << " bytes have been written" );
}

template<common::NumBasesPerByte numBasesPerByte>
FastqCreator<numBasesPerByte>::FastqCreator(
    TaskQueue& taskQueue,
    ThreadSafeQueue<std::shared_ptr<DemuxBuffer>>& inputBuffersToUse,
    ThreadSafeQueue<std::shared_ptr<DemuxBuffer>>& inputBuffersToRecycle,
    ThreadSafeMap<uint32_t, std::shared_ptr<FastqBuffer>>& outputBuffersToSubmit,
    ThreadSafeQueue<std::shared_ptr<FastqBuffer>>& outputBuffersToUse,
    const layout::Layout &layout,
    const layout::LaneInfo &laneInfo,
    std::size_t maskShortAdapterReads,
    float adapterStringency,
    bool generateReverseComplementFastqs,
    bool includeNonPfClusters,
    bool createFastqsForIndexReads,
    bool useBgzf,
    int compressionLevel,
    bool findAdaptersWithSlidingWindow,
    uint32_t numProcessingThreads,
    ConversionStatsForSampleTileBarcode &summaryTileStats,
    stats::BarcodeHits &unknownBarcodeHits
)
: Stage("Fastq creating", taskQueue)
, inputBuffersToUse_(inputBuffersToUse)
, inputBuffersToRecycle_(inputBuffersToRecycle)
, outputBuffersToSubmit_(outputBuffersToSubmit)
, outputBuffersToUse_(outputBuffersToUse)
, layout_(layout)
, laneInfo_(laneInfo)
, maskShortAdapterReads_(maskShortAdapterReads)
, adapterStringency_(adapterStringency)
, maskAdapters_()
, trimAdapters_()
, generateReverseComplementFastqs_(generateReverseComplementFastqs)
, includeNonPfClusters_(includeNonPfClusters)
, createFastqsForIndexReads_(createFastqsForIndexReads)
, useBgzf_(useBgzf)
, compressionLevel_(compressionLevel)
, findAdaptersWithSlidingWindow_(findAdaptersWithSlidingWindow)
, summaryTileStats_(summaryTileStats)
, unknownBarcodesHits_(unknownBarcodeHits)
, statsForEachThread_(numProcessingThreads)
, prevGroupNumber_(0)
, numBuffersCreated_(0)
{
    createTileStats();

    for (const auto& readInfo : laneInfo_.readInfos())
    {
        if (readInfo.shouldCreateFastqForRead(createFastqsForIndexReads_))
        {
            maskAdapters_.push_back(new boost::ptr_vector<AdapterLocator<numBasesPerByte>>());
            trimAdapters_.push_back(new boost::ptr_vector<AdapterLocator<numBasesPerByte>>());

            createAdapters(maskAdapters_.back(),
                           readInfo.maskAdaptersBegin(),
                           readInfo.maskAdaptersEnd());

            createAdapters(trimAdapters_.back(),
                           readInfo.trimAdaptersBegin(),
                           readInfo.trimAdaptersEnd());
        }
    }
}

template<common::NumBasesPerByte numBasesPerByte>
FastqCreator<numBasesPerByte>::~FastqCreator()
{
}

template<common::NumBasesPerByte numBasesPerByte>
void FastqCreator<numBasesPerByte>::createTileStats()
{
    size_t nonIndexReadLength = 0;
    size_t numNonIndexReads = 0;
    for (const auto& readInfo : laneInfo_.readInfos())
    {
        if (readInfo.isDataRead())
        {
            nonIndexReadLength = std::max(nonIndexReadLength, readInfo.cycleInfos().size());
            numNonIndexReads = std::max(numNonIndexReads, readInfo.getNumber());
        }
    }

    layout::LaneInfo::TileInfosContainer::size_type tilesCount = laneInfo_.getTileInfos().size();

    summaryTileStats_.resize(laneInfo_.getSampleInfos().size(),
                             ConversionStatsForTileBarcode(tilesCount));

    const auto& sampleInfos = laneInfo_.getSampleInfos();
    size_t sampleIndex = 0;
    for (auto& summaryStatsForSample : summaryTileStats_)
    {
        size_t tileIndex = 0;
        for (auto& summaryStatsForTile : summaryStatsForSample)
        {
            size_t barcodeIndex = 0;
            summaryStatsForTile.resize(std::max(static_cast<size_t>(1UL), sampleInfos.at(sampleIndex).getBarcodes().size()));
            for (auto& summaryStatsForBarcode : summaryStatsForTile)
            {
                summaryStatsForBarcode.first[stats::TileBarcodeStats::RAW].init(nonIndexReadLength, numNonIndexReads);
                summaryStatsForBarcode.first[stats::TileBarcodeStats::PF].init(nonIndexReadLength, numNonIndexReads);

                layout::BarcodeTranslationTable::SampleMetadata sampleMetadata(sampleIndex,barcodeIndex);
                for (const auto& readInfo : laneInfo_.readInfos())
                {
                    if (readInfo.isDataRead())
                    {
                        summaryStatsForBarcode.second.push_back(stats::ReadBarcodeStats(readInfo.getNumber()));
                    }
                }

                ++barcodeIndex;
            }

            ++tileIndex;
        }

        ++sampleIndex;
    }

    for (auto& statsForThread : statsForEachThread_)
    {
        statsForThread = std::make_shared<AllConversionStats>(
            summaryTileStats_,
            unknownBarcodesHits_);
    }
}

template<common::NumBasesPerByte numBasesPerByte>
bool FastqCreator<numBasesPerByte>::startNewTasks()
{
    static int sequentialId = 0; ++sequentialId;
    BCL2FASTQ_LOG(common::LogLevel::DEBUG) << "FastqCreator::startNewTask " << sequentialId  << std::endl;

    TaskQueue &taskQueue = this->getTaskQueue();
    if (taskQueue.isTerminated() || inputBuffersToUse_.isTerminated() )
    {
        // Only occurs if there was an error
        outputBuffersToSubmit_.terminate();
        inputBuffersToUse_.terminate();

        return false;
    }

    auto taskManager = getNewTaskManager();
    if (taskManager.get() == NULL)
    {
        return false;
    }

    FastqBuffer &outputBuffer = taskManager->getOutputBuffer();
    const auto& inputBuffer = taskManager->getInputBuffer().bclBuffers_;
    const auto& sampleIndex = *taskManager->getInputBuffer().sampleIndex_;

    layout::LaneInfo::SampleInfosContainer::size_type sampleInfosCounts = laneInfo_.getSampleInfos().size();

    // Check if tile file is missing. 
    DemuxBuffer & demuxBuffer  = taskManager->getInputBuffer();
    if (!demuxBuffer.hasCycleData()) 
    {
        // 'clear" the outputBuffer;
        for( layout::LaneInfo::SampleInfosContainer::size_type sampleInfoIndex = 0;
             sampleInfoIndex < sampleInfosCounts;
             ++sampleInfoIndex )
        {
// TODO: What's this supposed to do? It's a loop over reads, resizing the 0th element every time.
            const layout::ReadInfosContainer::size_type readInfoIndex = 0;
            for (const auto& readInfo : laneInfo_.readInfos())
            {
                outputBuffer.fastqs_.at(sampleInfoIndex).at(readInfoIndex).resize(0);
            }
        }
        return true;
    }

    BCL2FASTQ_ASSERT_MSG(inputBuffer.front().bcls_.front().bcls_.size() > 0, "No BCL cycles to process");
    BCL2FASTQ_LOG(common::LogLevel::DEBUG) << "Clusters/BCL: " << inputBuffer.front().bcls_.front().getNumClustersPerByte() * inputBuffer.front().bcls_.front().bcls_.size() << std::endl;

    uint32_t taskNum = 0;
    for( layout::LaneInfo::SampleInfosContainer::size_type sampleInfoIndex = 0;
         sampleInfoIndex < sampleInfosCounts;
         ++sampleInfoIndex )
    {
        size_t cycleIndex = 0;

        layout::ReadInfosContainer::size_type readInfoIndex = 0;
        for (const auto& readInfo : laneInfo_.readInfos())
        {
            // For index reads, use unmasked cycles.
            size_t nextCycleIndex = cycleIndex + readInfo.getNumCyclesToLoad();

           // cycleIndex += readInfo.getBclBufferOffset();
            if (readInfo.shouldCreateFastqForRead(createFastqsForIndexReads_))
            {
                SampleIndex::FastqOffsetsContainer::size_type offsetsIndex = 0;
                const SampleIndex::FastqOffsetsContainer::const_iterator offsetsEnd = sampleIndex.offsetsEnd(sampleInfoIndex);
                SampleIndex::FastqOffsetsContainer::const_iterator offset = sampleIndex.offsetsBegin(sampleInfoIndex);

                outputBuffer.fastqs_.at(sampleInfoIndex).at(readInfoIndex).resize((offsetsEnd-offset + data::ClustersPerTask - 1) / data::ClustersPerTask);

                while(offset != offsetsEnd)
                {
                    BCL2FASTQ_ASSERT_MSG(offsetsEnd - offset > 0, "Invalid offset iterator");
                    const SampleIndex::FastqOffsetsContainer::const_iterator nextOffset =
                        offset
                        +
                        (static_cast<std::size_t>(offsetsEnd - offset) < data::ClustersPerTask ? offsetsEnd - offset : data::ClustersPerTask);

                    taskQueue.addData(std::make_shared<FastqCreateTask<numBasesPerByte>>(
                        taskManager,
                        inputBuffer,
                        cycleIndex,
                        nextCycleIndex,
                        layout_.getFlowcellInfo(),
                        laneInfo_,
                        readInfo,
                        offset,
                        nextOffset,
                        maskShortAdapterReads_,
                        adapterStringency_,
                        maskAdapters_[readInfoIndex],
                        trimAdapters_[readInfoIndex],
                        generateReverseComplementFastqs_,
                        includeNonPfClusters_,
                        useBgzf_,
                        compressionLevel_,
                        findAdaptersWithSlidingWindow_,
                        statsForEachThread_,
                        outputBuffer.fastqs_.at(sampleInfoIndex).at(readInfoIndex).at(offsetsIndex)
                    ));

                    ++taskNum;
                    ++offsetsIndex;
                    offset = nextOffset;
                }
                ++readInfoIndex;
            }

            cycleIndex = nextCycleIndex;
        }
    }

    return true;
}

template<common::NumBasesPerByte numBasesPerByte>
void FastqCreator<numBasesPerByte>::terminate()
{
    this->getTaskQueue().terminate();
    inputBuffersToUse_.terminate();
}

template<common::NumBasesPerByte numBasesPerByte>
std::shared_ptr<FastqCreateTaskManager> FastqCreator<numBasesPerByte>::getNewTaskManager()
{
    std::shared_ptr<DemuxBuffer> inputBuffer;
    if (!inputBuffersToUse_.tryGetData(inputBuffer))
    {
        FastqCreateTaskManager::waitForAllTasksToFinish();

        // We're done. Notify the next stage (fastq writing).
        outputBuffersToSubmit_.setFinished(prevGroupNumber_);
        return std::shared_ptr<FastqCreateTaskManager>();
    }

    prevGroupNumber_ = inputBuffer->getGroupNumber();
    std::shared_ptr<FastqBuffer> outputBuffer;

    //const bool allBuffersWaiting = (numBuffersCreated_ >= 3 && outputBuffersToUse_.isEmpty());
    if (numBuffersCreated_ < 3) // || allBuffersWaiting)
    {
        if (!outputBuffersToUse_.tryGetDataCheckEmpty(outputBuffer))
        {
            // There aren't any buffers available. Make a new one.
            // The memory allocation isn't too bad here, since we're recycling the buffers.

            outputBuffer = std::make_shared<FastqBuffer>();
            ++numBuffersCreated_;

            layout::LaneInfo::SampleInfosContainer::size_type sampleInfosCounts = laneInfo_.getSampleInfos().size();

            outputBuffer->fastqs_.resize(sampleInfosCounts);

            layout::ReadInfosContainer::size_type readsCount = std::count_if(laneInfo_.readInfos().begin(),
                                                                             laneInfo_.readInfos().end(),
                                                                             [this] (const layout::ReadInfo& readInfo)
                                                                                 { return !readInfo.cycleInfos().empty() &&
                                                                                          (createFastqsForIndexReads_ || readInfo.isDataRead()); });

            for( layout::LaneInfo::SampleInfosContainer::size_type sampleInfoIndex = 0;
                 sampleInfoIndex < sampleInfosCounts;
                 ++sampleInfoIndex )
            {
                outputBuffer->fastqs_[sampleInfoIndex].resize(readsCount);
            }
        }
    }
    else
    {
        if (!outputBuffersToUse_.tryGetData(outputBuffer))
        {
            return std::shared_ptr<FastqCreateTaskManager>();
        }
    }

    outputBuffer->setGroupNumber(inputBuffer->getGroupNumber());

    return std::make_shared<FastqCreateTaskManager>(inputBuffer,
                                                    inputBuffersToRecycle_,
                                                    outputBuffer,
                                                    outputBuffersToSubmit_);
}

template<common::NumBasesPerByte numBasesPerByte>
void FastqCreator<numBasesPerByte>::createAdapters(boost::ptr_vector<AdapterLocator<numBasesPerByte>>& adapters,
                                                   layout::ReadInfo::AdaptersContainer::const_iterator begin,
                                                   layout::ReadInfo::AdaptersContainer::const_iterator end) const
{
    size_t minimumTrimmedReadLength = laneInfo_.getMinimumTrimmedReadLength();
    BOOST_FOREACH (const layout::ReadInfo::AdaptersContainer::value_type& adapter, std::make_pair(begin, end))
    {
        findAdaptersWithSlidingWindow_ ?
            adapters.push_back(new AdapterLocatorSlidingWindow<numBasesPerByte>(adapter, adapterStringency_)) :
            adapters.push_back(new AdapterLocatorWithIndels<numBasesPerByte>(adapter, adapterStringency_, minimumTrimmedReadLength));
    }
}


} // namespace conversion
} // namespace bcl2fastq

#endif // BCL2FASTQ_CONVERSION_FASTQ_CREATOR_HPP

