/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file FastqCreator.hh
 *
 * \brief Declaration of FASTQ creator.
 *
 * \author Marek Balint
 * \author Mauricio Varea
 */


#ifndef BCL2FASTQ_CONVERSION_FASTQCREATOR_HH
#define BCL2FASTQ_CONVERSION_FASTQCREATOR_HH

#include "layout/Layout.hh"
#include "stats/TileStats.hpp"
#include "stats/BarcodeHits.hh"
#include "data/BclBuffer.hh"
#include "conversion/FastqBuffer.hh"
#include "conversion/Stage.hh"
#include "conversion/Task.hh"
#include "conversion/FastqIterator.hh"
#include "conversion/AdapterLocator.hh"
#include "conversion/SampleIndex.hh"

#include <boost/ptr_container/ptr_vector.hpp>

namespace bcl2fastq
{

namespace io
{
    class GzipCompressor;
}

namespace conversion
{

class FastqCreateTaskManager;

/// \brief BaseCalling statistics (one per tile, sample, barcode).
typedef std::pair< stats::TileBarcodeStats,
                   std::vector<stats::ReadBarcodeStats> > ConversionStats;

// First index is tile, second is sample, third is barcode
typedef std::vector<std::vector<std::vector<ConversionStats>>> ConversionStatsForSampleTileBarcode;

typedef ConversionStatsForSampleTileBarcode::value_type ConversionStatsForTileBarcode;

// There will be 1 of these objects for every processing thread.
struct AllConversionStats : private boost::noncopyable
{
    AllConversionStats(ConversionStatsForSampleTileBarcode& summaryConversionStats,
                       stats::BarcodeHits& summaryUnknownBarcodes);

    // Record to the summary stats in the destructor
    ~AllConversionStats();

private:
    ConversionStatsForSampleTileBarcode& summaryConversionStats_;
    stats::BarcodeHits& summaryUnknownBarcodes_;

public:
    ConversionStatsForSampleTileBarcode conversionStats_;
    stats::BarcodeHits unknownBarcodes_;
};

// Stats for all tasks (index is task number)
    typedef std::vector<std::shared_ptr<AllConversionStats>> FastqCreateStats;

    typedef std::shared_ptr<FastqCreateStats> FastqCreateStatsPtr;


/// \brief Task: Convert BCLs to FASTQs.
template<common::NumBasesPerByte numBasesPerByte>
class FastqCreateTask : public Task
{
public:

    /// \brief Constructor.
    /// \param bclBuffer BCL buffer.
    /// \param cyclesBegin Beginning of BCL cycles to process.
    /// \param cyclesEnd End of BCL cycles to process.
    /// \param flowcellInfo Flowcell meta data.
    /// \param laneInfo Lane meata data.
    /// \param readInfo Read meta data.
    /// \param offsetNumber Offset number.
    /// \param offsetsBegin Beginning of FASTQ offsets range.
    /// \param offsetsEnd End of FASTQ offsets range.
    /// \param maskShortAdapterReads Maximum number of useful bases in read after adapter trimming for which whole read is masked.
    /// \param adapterStringency Adapter stringency.
    /// \param maskAdapters Mask adapters.
    /// \param trimAdapters Trim adapters.
    /// \param generateReverseComplementFastqs Generate reverse complement FASTQs flag.
    /// \param includeNonPfClusters Include non-PF clusters in created FASTQ files flag.
    /// \param layout Flowcell layout.
    /// \param useBgzf If true, use BGZF compression
    /// \param compressionLevel Compression level used by zlib
    /// \param findAdaptersWithSlidingWindow If true, find adapters with the sliding window algorithm.
    /// \param tileStats Statistics collected for a tile
    /// \param outputBuffer FASTQ output buffer.
    FastqCreateTask(
        std::shared_ptr<FastqCreateTaskManager>& taskManager,
        const data::BclBufferVec& bclBuffers,
        size_t cyclesIndex,
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
        FastqCreateStats &stats,
        FastqBuffer::FastqsContainer::value_type::value_type::value_type &outputBuffer
    );

    virtual PriorityLevel getPriority() const { return PriorityLevel::Three; }

    virtual bool execute(int32_t threadNum);

private:

    /// \brief Update tileStats_ accordingly.
    /// \param basesBegin Beginning of data.
    /// \param basesEnd End of data.
    /// \param barcode Index for this read.
    /// \param filterFlag Whether data passes chastity filter or not.
    /// \param trimmedCount Number of trimmed bases.
    /// \return Whether or not the barcode has been found in the table.
    bool computeStatistics(
        const FastqConstIterator<numBasesPerByte>& basesBegin,
        const FastqConstIterator<numBasesPerByte>& basesEnd,
        SampleIndex::FastqOffsetsContainer::value_type offset,
        int32_t threadNum,
        bool filterFlag,
        size_t trimmedCount,
        const std::vector<uint32_t>& remappedQscores
    );

    /// \brief Crate FASTQ.
    /// \param offset Offset of the FASTQ in the input buffer.
    /// \param compressor Gzip compressor.
    void fastqCreate(
        SampleIndex::FastqOffsetsContainer::value_type offset,
        int32_t threadNum,
        io::GzipCompressor& compressor
    );

    /// \brief Create the barcode strings for the current read.
    /// \param barcodeStrings Vector of barcode strings.
    /// \param offset Offset of the FASTQ in the input buffer.
    void createBarcodeStrings(std::vector<std::string>&                      barcodeStrings,
                              SampleIndex::FastqOffsetsContainer::value_type offset);

    /// \brief Write a header element.
    /// \param element Element to write.
    void writeHeaderElement(const std::string& element);

    /// \brief Create the header.
    /// \param offset Offset of the FASTQ in the input buffer.
    /// \param filterFlag True if filter passed.
    void createHeader(SampleIndex::FastqOffsetsContainer::value_type offset,
                      bool                                           filterFlag);

    /// \brief Create the bases and quality scores.
    /// \param offset Offset of the FASTQ in the input buffer.
    /// \param compressor Gzip compressor.
    /// \param trimmedCount Number of trimmed bases.
    void createBasesAndQualities(SampleIndex::FastqOffsetsContainer::value_type offset,
                                 io::GzipCompressor& compressor,
                                 size_t& trimmedCount);

private:

    /// \brief BCL buffer.
    const data::BclBufferVec& bclBuffers_;

    /// \brief Beginning of BCL cycles to process.
    const size_t cycleIndex_;

    /// \brief End of BCL cycles to process.
    const size_t cycleIndexEnd_;

    /// \brief Flowcell meta data.
    const layout::FlowcellInfo &flowcellInfo_;

    /// \brief Lane meta data.
    const layout::LaneInfo &laneInfo_;

    /// \brief Read meta data.
    const layout::ReadInfo &readInfo_;

    /// \brief Sample number.
    const layout::BarcodeTranslationTable::SampleMetadata sampleMetadata_;

    /// \brief Beginning of FASTQ offsets range.
    const SampleIndex::FastqOffsetsContainer::const_iterator offsetsBegin_;

    /// \brief End of FASTQ offsets range.
    const SampleIndex::FastqOffsetsContainer::const_iterator offsetsEnd_;

    /// \brief Buffer to write FASTQs to.
    FastqBuffer::FastqsContainer::value_type::value_type::value_type &outputBuffer_;

    /// \brief Intermediate buffer.
    FastqBuffer::FastqsContainer::value_type::value_type::value_type buffer_;

    /// \brief Use BGZF compression if true.
    bool useBgzf_;

    /// \brief Compression level.
    int compressionLevel_;

    /// \brief Find adapters with the sliding window algorithm if true.
    bool findAdaptersWithSlidingWindow_;

    /// \brief Minimum read length after adapter trimming.
    std::size_t minimumTrimmedReadLength_;

    /// \brief Maximum number of useful bases in read after adapter trimming for which whole read is masked.
    std::size_t maskShortAdapterReads_;

    /// \brief Adapter stringency.
    float adapterStringency_;

    /// \brief Mask adapters.
    const boost::ptr_vector<AdapterLocator<numBasesPerByte>>& maskAdapters_;

    /// \brief Trim adapters.
    const boost::ptr_vector<AdapterLocator<numBasesPerByte>>& trimAdapters_;

    /// \brief Generate reverse complement FASTQs flag.
    bool generateReverseComplementFastqs_;

    /// \brief Include non-PF clusters in created FASTQ files.
    bool includeNonPfClusters_;

    /// \brief Conversion statistics.
    FastqCreateStats &stats_;

    /// \brief UMI cycles.
    std::vector< common::CycleRange > umiCycles_;
};

class FastqCreateTaskManager : public TaskManager
{
public:
    FastqCreateTaskManager(std::shared_ptr<DemuxBuffer>& inputBuffer,
                           ThreadSafeQueue<std::shared_ptr<DemuxBuffer>>& inputBuffersToRecycle,
                           std::shared_ptr<FastqBuffer>& outputBuffer,
                           ThreadSafeMap<uint32_t, std::shared_ptr<FastqBuffer>>& outputBufferMap);

    virtual ~FastqCreateTaskManager();

    static void waitForAllTasksToFinish();

    DemuxBuffer& getInputBuffer() { return *inputBuffer_; }
    FastqBuffer& getOutputBuffer() { return *outputBuffer_; }

private:
    std::shared_ptr<DemuxBuffer> inputBuffer_;
    ThreadSafeQueue<std::shared_ptr<DemuxBuffer>>& inputBuffersToRecycle_;
    std::shared_ptr<FastqBuffer> outputBuffer_;
    ThreadSafeMap<uint32_t, std::shared_ptr<FastqBuffer>>& outputBufferMap_;

    static std::atomic<uint32_t> numTaskManagers_;
    static std::condition_variable cvAllTaskManagersDone_;
    static std::mutex mut_;
};

/// \brief FASTQ creator.
template<common::NumBasesPerByte numBasesPerByte>
class FastqCreator : public Stage
{
public:

    /// \brief Constructor.
    /// \param threadsCount Number of threads.
    /// \param layout Flowcell layout.
    /// \param laneInfo Lane meata data.
    /// \param maskShortAdapterReads Maximum number of useful bases in read after adapter trimming for which whole read is masked.
    /// \param adapterStringency Adapter stringency.
    /// \param generateReverseComplementFastqs Generate reverse complement FASTQs.
    /// \param includeNonPfClusters Include non-PF clusters in created FASTQ files flag.
    /// \param createFastqsForIndexReads Create FASTQ files also for index reads flag.
    /// \param useBgzf If true, use BGZF compression
    /// \param compressionLevel Compression level used by zlib
    /// \param findAdaptersWithSlidingWindow If true, find adapters with the sliding window algorithm.
    /// \param tileStats Statistics collected for a tile
    /// \param unknownBarcodeHits Counts of unknown barcodes
    FastqCreator(
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
    );

    virtual ~FastqCreator();

public:

    virtual bool startNewTasks();

private:

    virtual void terminate();

    virtual std::shared_ptr<FastqCreateTaskManager> getNewTaskManager();

    void createTileStats();

    void createAdapters(boost::ptr_vector<AdapterLocator<numBasesPerByte>>& adapters,
                        layout::ReadInfo::AdaptersContainer::const_iterator begin,
                        layout::ReadInfo::AdaptersContainer::const_iterator end) const;

private:

    ThreadSafeQueue<std::shared_ptr<DemuxBuffer>>& inputBuffersToUse_;
    ThreadSafeQueue<std::shared_ptr<DemuxBuffer>>& inputBuffersToRecycle_;
    ThreadSafeMap<uint32_t, std::shared_ptr<FastqBuffer>>& outputBuffersToSubmit_;
    ThreadSafeQueue<std::shared_ptr<FastqBuffer>>& outputBuffersToUse_;

    /// \brief Layout.
    const layout::Layout &layout_;

    /// \brief Current lane.
    const layout::LaneInfo &laneInfo_;

    /// \brief Maximum number of useful bases in read after adapter trimming for which whole read is masked.
    std::size_t maskShortAdapterReads_;

    /// \brief Adapter stringency.
    float adapterStringency_;

    /// \brief Mask adapters.
    boost::ptr_vector< boost::ptr_vector<AdapterLocator<numBasesPerByte>> > maskAdapters_;

    /// \brief Trim adapters.
    boost::ptr_vector< boost::ptr_vector<AdapterLocator<numBasesPerByte>> > trimAdapters_;

    /// \brief Generate reverse complement FASTQs flag.
    bool generateReverseComplementFastqs_;

    /// \brief Include non-PF clusters in created FASTQ files.
    bool includeNonPfClusters_;

    /// \brief Create FASTQ files also for index reads.
    bool createFastqsForIndexReads_;

    /// \brief Use BGZF compression if true.
    bool useBgzf_;

    /// \brief Compression level used by zlib
    int compressionLevel_;

    /// \brief Find adapters with the sliding window algorithm if true.
    bool findAdaptersWithSlidingWindow_;

    /// \brief Conversion statistics (summary).
    ConversionStatsForSampleTileBarcode &summaryTileStats_;

    /// \brief Unknown barcodes statistics (summary).
    stats::BarcodeHits &unknownBarcodesHits_;

    FastqCreateStats statsForEachThread_;

    uint32_t prevGroupNumber_;

    uint32_t numBuffersCreated_;
};


} // namespace conversion
} // namespace bcl2fastq

#include "conversion/FastqCreator.hpp"

#endif // BCL2FASTQ_CONVERSION_FASTQCREATOR_HH


