/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file Demultiplexer.hh
 *
 * \brief Declaration of demultiplexer.
 *
 * \author Marek Balint
 * \author Mauricio Varea
 */


#ifndef BCL2FASTQ_CONVERSION_DEMULTIPLEXER_HH
#define BCL2FASTQ_CONVERSION_DEMULTIPLEXER_HH


#include "data/BclBuffer.hh"
#include "layout/Layout.hh"
#include "layout/LaneInfo.hh"
#include "layout/BarcodeTranslationTable.hh"
#include "stats/DemultiplexingStatsXml.hh"
#include "stats/BarcodeStats.hh"
#include "conversion/Stage.hh"
#include "conversion/Task.hh"
#include "conversion/SampleIndex.hh"

#include "conversion/FastqIterator.hh"

namespace bcl2fastq
{
namespace conversion
{

class DemuxTaskManager;

/// \brief Task: Demultiplex.
class DemultiplexTask : public Task
{
public:

    /// \brief Demultiplexing statistics. One stats::BarcodeStats per sample, barcode, tile.    
    typedef std::vector<std::vector<std::vector<stats::BarcodeStats>>> DemuxStats;

    class DemuxStatsSubset
    {
    public:
        DemuxStatsSubset(DemuxStats& summaryDemuxStats);
        ~DemuxStatsSubset();

        DemuxStats& getStats() { return stats_; }

    private:
        DemuxStats& summaryDemuxStats_;
        DemuxStats stats_;
    };

public:

    /// \brief Constructor.
    /// \param inputBuffer Input buffer.
    /// \param readsInfos Read infos.
    /// \param offsetBegin Starting offset.
    /// \param offsetEnd Ending offset (one past the end).
    /// \param barcodeTranslationTable Barcode translation table.
    /// \param tileInfosBegin tile info
    /// \param numIndexReads number of index reads
    /// \param includeNonPfClusters If true, include non PF clusters.
    /// \param numBasesPerByte Number of bases in a byte.
    /// \param demuxStats Demultiplexing statistics.
    /// \param outputBuffer Output buffer.
    DemultiplexTask(
        std::shared_ptr<DemuxTaskManager>& taskManager,
        const data::BclBuffer &inputBuffer,
        const layout::ReadInfosContainer& readInfos,
        data::BclBuffer::BclCycleContainer::size_type offsetBegin,
        data::BclBuffer::BclCycleContainer::size_type offsetEnd,
        const layout::BarcodeTranslationTable &barcodeTranslationTable,
        layout::LaneInfo::TileInfosContainer::const_iterator tileInfosBegin,
        size_t numIndexReads,
        bool includeNonPfClusters,
        common::NumBasesPerByte numBasesPerByte,
        std::vector<std::shared_ptr<DemultiplexTask::DemuxStatsSubset>>& demuxStats,
        data::BclBuffer::SamplesContainer &outputBuffer
    );

    virtual ~DemultiplexTask() { }

    // Demux has a serial step. Give it the highest priority to avoid waiting on this.
    virtual PriorityLevel getPriority() const { return PriorityLevel::Zero; }

    virtual bool execute(int32_t threadNum);

    template<common::NumBasesPerByte numBasesPerByte>
    bool execute(int32_t threadNum);

private:

    /// \brief Input buffer.
    const data::BclBuffer &inputBuffer_;

    /// \brief Read infos.
    const layout::ReadInfosContainer& readInfos_;

    /// \brief Starting offset.
    const data::BclBuffer::BclCycleContainer::size_type offsetBegin_;

    /// \brief Ending offset (one past the end).
    const data::BclBuffer::BclCycleContainer::size_type offsetEnd_;

    /// \brief Barcode to sample translation table.
    const layout::BarcodeTranslationTable &barcodeTranslationTable_;

    /// \brief Reference to the first tile, used to calculate index.
    layout::LaneInfo::TileInfosContainer::const_iterator tileInfosBegin_;

    /// \brief Number of index reads.
    size_t numIndexReads_;

    /// \brief Include non-PF clusters
    bool includeNonPfClusters_;

    /// \brief Number of bases per byte.
    common::NumBasesPerByte numBasesPerByte_;

    /// \brief Demultiplexing statistics.
    std::vector<std::shared_ptr<DemultiplexTask::DemuxStatsSubset>>& demuxStats_;

    /// \brief Output buffer.
    data::BclBuffer::SamplesContainer &outputBuffer_;
};

class DemuxTaskManager : public TaskManager
{
public:
    DemuxTaskManager(std::shared_ptr<DemuxBuffer>& buffer,
                     ThreadSafeQueue<std::shared_ptr<DemuxBuffer>>& outputBufferQueue,
                     TaskManager::GroupSerializer &groupSerializer);

    virtual ~DemuxTaskManager();

    DemuxBuffer& getBuffer() { return *buffer_; }

    static void waitForAllTasksToFinish();

private:
    std::shared_ptr<DemuxBuffer> buffer_;
    ThreadSafeQueue<std::shared_ptr<DemuxBuffer>>& outputBufferQueue_;

    TaskManager::GroupSerializer &groupSerializer_;

    static std::atomic<uint32_t> numTaskManagers_;
    static std::condition_variable cvAllTaskManagersDone_;
    static std::mutex mut_;
};

/// \brief Demultiplexer
class Demultiplexer : public Stage
{
public:

    /// \brief Maximum number of allowed mismatches in a barcode component type definition.
    typedef std::size_t MismatchesCount;

    /// \brief Maximum allowed mismatches mismatches in a barcode component container type definition.
    typedef std::vector<MismatchesCount> BarcodeMismatchCountsContainer;

public:

    /// \brief Constructor.
    /// \param threadsCount Number of threads.
    /// \param layout Flowcell layout.
    /// \param laneInfo Lane meata data.
    /// \param barcodeTranslationTable Table to translate barcodes to samples.
    /// \includeNonPfClusters If true, include non PF clusters.
    /// \param interopDir Interop directory.
    Demultiplexer(
        TaskQueue& taskQueue,
        ThreadSafeQueue<std::shared_ptr<DemuxBuffer>>& inputBuffersToUse,
        ThreadSafeQueue<std::shared_ptr<DemuxBuffer>>& outputBuffersToSubmit,
        const layout::Layout &layout,
        const layout::LaneInfo &laneInfo,
        const layout::BarcodeTranslationTable &barcodeTranslationTable,
        bool includeNonPfClusters,
        int32_t numThreads,
        DemultiplexTask::DemuxStats &summaryDemuxStats
    );

    virtual ~Demultiplexer();

    virtual bool startNewTasks();

    void calcFilterOffsets(data::BclBufferVec& buffers);

private:

    virtual void terminate();

    virtual std::shared_ptr<DemuxTaskManager> getNewTaskManager();

    void createTileStats(const layout::BarcodeTranslationTable::SampleMetadata &sampleMetadata);

    ThreadSafeQueue<std::shared_ptr<DemuxBuffer>>& inputBuffersToUse_;
    ThreadSafeQueue<std::shared_ptr<DemuxBuffer>>& outputBuffersToSubmit_;

    /// \brief Layout.
    const layout::Layout &layout_;

    /// \brief Current lane.
    const layout::LaneInfo &laneInfo_;

    /// \brief Barcode to sample translation table.
    const layout::BarcodeTranslationTable &barcodeTranslationTable_;

    /// \brief Include non-PF clusters
    bool includeNonPfClusters_;

    /// \brief Number of threads.
    common::ThreadVector::size_type threadsCount_;

    /// \brief Demultiplexing statistics (summary).
    DemultiplexTask::DemuxStats &summaryDemuxStats_;

    /// \brief Flag determining whether demultiplexing should take place.
    bool doDemultiplex_;

    /// \brief Number of index reads.
    size_t numIndexReads_;

    std::vector<std::shared_ptr<DemultiplexTask::DemuxStatsSubset>> stats_;

    TaskManager::GroupSerializer groupSerializer_;
};


} // namespace conversion
} // namespace bcl2fastq

#include "conversion/Demultiplexer.hpp"

#endif // BCL2FASTQ_CONVERSION_DEMULTIPLEXER_HH


