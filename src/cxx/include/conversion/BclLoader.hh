/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file BclLoader.hh
 *
 * \brief Declaration of BCL loader.
 *
 * \author Marek Balint
 * \author Mauricio Varea
 */


#ifndef BCL2FASTQ_CONVERSION_BCLLOADER_HH
#define BCL2FASTQ_CONVERSION_BCLLOADER_HH


#include "conversion/BclReader.hh"

#include "layout/Layout.hh"
#include "data/BclFile.hh"
#include "data/CycleBCIFile.hh"
#include "data/PositionsFile.hh"
#include "data/BclBuffer.hh"
#include "conversion/Stage.hh"
#include "conversion/Task.hh"
#include "conversion/SampleIndex.hh"

#include <boost/filesystem/path.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <vector>
#include <boost/shared_ptr.hpp>
#include <mutex>
#include <atomic>

namespace bcl2fastq {


namespace conversion {

class BclLoaderTaskManager;

/// \brief Task: Load BCL data into buffer.
class BclLoadTask : public Task
{
public:

    /// \brief Constructor.
    /// \param bclFile BCL file to be loaded.
    /// \param cycleBciFile Cycle BCI file.
    /// \param aggregateTilesFlag All tiles aggregated into single file flag.
    /// \param laneInfo Lane meata data.
    /// \param cycleInfo Cycle meta data.
    /// \param clustersCount Number of clusters to be loaded.
    /// \param ignoreMissingBcls Ignore missing BCLs flag.
    /// \param bclIdx bcl index
    /// \param outputBuffer Buffer for data to be loaded to.
    BclLoadTask(
        std::shared_ptr<BclLoaderTaskManager>& taskManager,
        data::RawBclBufferGroup& bclData,
        std::mutex* aggregatedBclFileMutex,
        common::TileAggregationMode aggregateTilesFlag,
        const layout::LaneInfo &laneInfo,
        const layout::CycleInfo &cycleInfo,
        layout::LaneInfo::TileInfosContainer::const_iterator tileInfoIter,
        bool ignoreMissingBcls,
        size_t bclIdx,
        DemuxBuffer& outputBuffer
    );

public:

    virtual PriorityLevel getPriority() const { return PriorityLevel::Two; }

    virtual bool execute(int32_t threadNum);

private:

    void checkBytesRead(int32_t bytesRead,
                        int32_t expectedBytes,
                        int32_t tileIdx,
                        const boost::filesystem::path& bclPath);

    /// \brief Read all tiles (each tile from a separate file).
    bool readSeparateTiles();

    /// \brief Read all tiles (each tile from the same file).
    bool readAggregatedTiles();

    /// \brief Read all tiles (each tile from a CBCL file).
    bool readCbclTiles();

    std::stringstream &getTileDescriptor(std::stringstream &msg, data::PerCycleData perCycleData);

    /// \brief BCL file to be loaded.
    data::RawBclBufferGroup& bclData_;

    std::mutex* aggregatedBclFileMutex_;

    /// \brief ALl tiles aggregated into single file flag.
    common::TileAggregationMode aggregateTilesMode_;

    /// \brief Lane meta data.
    const layout::LaneInfo &laneInfo_;

    /// \brief Cycle meta data.
    const layout::CycleInfo &cycleInfo_;

    /// \brief Tile meta data.
    layout::LaneInfo::TileInfosContainer::const_iterator tileInfoIter_;

    /// \brief Ignore missing BCLs flag.
    bool ignoreMissingBcls_;

    /// \brief Index into bcl buffer
    size_t bclIdx_;

    /// \brief Buffer to load data to.
    DemuxBuffer &outputBuffer_;
};


/// \brief Task: Load positions data into buffer.
class PositionsLoadTask : public Task
{
public:

    /// \brief Constructor.
    /// \param positionsFile Positions file.
    /// \param aggregateTilesFlag All tiles aggregated into single file flag.
    /// \param laneInfo Lane meata data.
    /// \param tileInfo Tile meta data.
    /// \param clustersCount Number of clusters to be loaded.
    /// \param ignoreMissingPositions Ignore missing positions flag.
    /// \param outputBuffer Buffer for data to be loaded to.
    PositionsLoadTask(
        std::shared_ptr<BclLoaderTaskManager>& taskManager,
        data::RawBclBufferGroup& bclData,
        bool aggregateTilesFlag,
        const layout::LaneInfo &laneInfo,
        layout::LaneInfo::TileInfosContainer::const_iterator &tileInfo,
        const std::vector<common::ClustersCount>& clustersCounts,
        bool ignoreMissingPositions,
        data::BclBufferVec& outputBuffer,
        std::shared_ptr<data::BclBuffer::PatternedPositionsContainer>& patternedFlowcellPositions
    );

public:

    // Since the positions could be shared if using a patterned flowcell, the priority is a bit higher.
    virtual PriorityLevel getPriority() const { return PriorityLevel::One; }

    virtual bool execute(int32_t threadNum);

    bool execute(const common::RawDataBuffer& inputBuffer,
                 data::BclBuffer::PositionsContainer& outputBuffer);

private:

    /// \brief BCL file to be loaded.
    data::RawBclBufferGroup& bclData_;

    /// \brief All tiles aggregated into single file flag.
    bool aggregateTilesFlag_;

    /// \brief Lane meta data.
    const layout::LaneInfo &laneInfo_;

    /// \brief Tile meta data.
    layout::LaneInfo::TileInfosContainer::const_iterator tileInfoIter_;

    /// \brief Number of clusters to load.
    const std::vector<common::ClustersCount>& clustersCounts_;

    /// \brief Ignore missing positions flag.
    bool ignoreMissingPositions_;

    /// \brief Buffer to load data to.
    data::BclBufferVec &outputBuffer_;

    /// \brief Buffer for patterned flowcell positions
    std::shared_ptr<data::BclBuffer::PatternedPositionsContainer>& patternedFlowcellPositions_;
};


/// \brief Task: Load filter data into buffer.
class FilterLoadTask : public Task
{
public:

    /// \brief Constructor.
    /// \param filterFile Filter file.
    /// \param aggregateTilesFlag All tiles aggregated into single file flag.
    /// \param laneInfo Lane meata data.
    /// \param tileInfo Tile meta data.
    /// \param clustersCount Number of clusters to be loaded.
    /// \param ignoreMissingFilters Ignore missing filters flag.
    /// \param outputBuffer Buffer for data to be loaded to.
    /// \param controlsOutputBuffer Buffer for controls data to be loaded to.
    FilterLoadTask(
        std::shared_ptr<BclLoaderTaskManager>& taskManager,
        data::RawBclBufferGroup& bclData,
        common::TileAggregationMode tileAggregationMode,
        const layout::LaneInfo &laneInfo,
        layout::LaneInfo::TileInfosContainer::const_iterator& tileInfo,
        const std::vector<common::ClustersCount>& clustersCounts,
        bool ignoreMissingFilters,
        data::BclBufferVec& outputBuffer
    );

public:

    virtual PriorityLevel getPriority() const { return PriorityLevel::Two; }

    virtual bool execute(int32_t threadNum);

private:

    /// \brief Read the aggregate tiles filter file.
    /// \param filePath Path to file.
    void readAggregateTilesFilterFile(const layout::TileInfo& tileInfo,
                                      common::ClustersCount clustersCount,
                                      data::BclBuffer& outputBuffer,
                                      const boost::filesystem::path& filePath = boost::filesystem::path());

    /// \brief Read the filter file for a single tile.
    /// \param filePath Path to file.
    void readTileFilterFile(const boost::filesystem::path& filePath,
                            const layout::TileInfo& tileInfo,
                            data::BclBuffer& outputBuffer);

    /// \brief Read the filter file.
    /// \param filterFile File to read from.
    void readFilterFile(common::RawDataBuffer& filterData,
                        const layout::TileInfo& tileInfo,
                        data::BclBuffer& outputBuffer,
                        bool isAggregatedTiles,
                        bool skipHeader);

    /// \brief Validate the filter records
    /// \param recordsRead Number of records read from file.
    /// \param filePath Path to file.
    void validateFilterRecords(std::size_t                    recordsRead,
                               common::ClustersCount          clustersCount,
                               const boost::filesystem::path& filePath,
                               data::BclBuffer&               outputBuffer);

    data::RawBclBufferGroup& bclData_;

    /// \brief File aggregation type.
    common::TileAggregationMode tileAggregationMode_;
    /// \brief All tiles aggregated into single file flag.
    bool aggregateTilesFlag_;

    /// \brief Lane meta data.
    const layout::LaneInfo &laneInfo_;

    /// \brief Tile meta data.
    layout::LaneInfo::TileInfosContainer::const_iterator tileInfoIter_;

    /// \brief Number of clusters to load.
    const std::vector<common::ClustersCount>& clustersCounts_;

    /// \brief Ignore missing filters flag.
    bool ignoreMissingFilters_;

    /// \brief Buffer to load control data to.
    data::BclBufferVec& outputBuffer_;
};


class ControlLoadTask : public Task
{
public:
    ControlLoadTask(
        std::shared_ptr<BclLoaderTaskManager>& taskManager,
        data::RawBclBufferGroup& bclData,
        const layout::LaneInfo &laneInfo,
        bool ignoreMissingControls,
        data::BclBufferVec& outputBuffer);

    virtual bool execute(int32_t threadNum);

    virtual PriorityLevel getPriority() const { return PriorityLevel::Two; }

private:
    data::RawBclBufferGroup& bclData_;

    /// \brief Lane meta data.
    const layout::LaneInfo &laneInfo_;

    /// \brief Tile meta data.
    layout::LaneInfo::TileInfosContainer::const_iterator tileInfoIter_;

    /// \brief Ignore missing controls flag.
    bool ignoreMissingControls_;

    /// \brief Buffer to load data to.
    data::BclBufferVec& outputBuffer_;
};

class BclLoaderTaskManager : public TaskManager
{
public:
    BclLoaderTaskManager(bool ignoreMissingBcls,
                         bool ignoreMissingFilters,
                         bool ignoreMissingPositions,
                         bool ignoreMissingControls,
                         common::TileAggregationMode aggregateTilesFlag,
                         std::shared_ptr<data::RawBclBufferGroup>& inputBuffer,
                         ThreadSafeQueue<std::shared_ptr<data::RawBclBufferGroup>>& inputBuffersToRecycle,
                         std::shared_ptr<DemuxBuffer>& outputBuffer,
                         ThreadSafeQueue<std::shared_ptr<DemuxBuffer>>& outputBufferQueue,
                         std::shared_ptr<data::BclBuffer::PatternedPositionsContainer>& patternedFlowcellPositions,
                         TaskManager::GroupSerializer &groupSerializer);

    virtual ~BclLoaderTaskManager();

    std::shared_ptr<data::RawBclBufferGroup> getInputBuffer() { return inputBuffer_;}
    std::shared_ptr<DemuxBuffer> getOutputBuffer() { return outputBuffer_; }

    static void waitForAllTasksToFinish();

private:
    void postExecute(DemuxBuffer& postExBuffer);

    void bclMismatchCount(const std::string& fileType,
                          common::CycleNumber cycleNumber,
                          common::TileNumber tileNumber,
                          common::RawDataBuffer::size_type realSize,
                          common::RawDataBuffer::size_type expectedSize);

    void createUniqueFakePositions(data::BclBuffer& outputBuffer,
                                   data::PerCycleData::BclsContainer::size_type bufferSize);

    void validateCycleMetadata(data::PerCycleData& cycleData,
                               data::PerCycleData& prevCycleData);

    bool ignoreMissingBcls_;
    bool ignoreMissingFilters_;
    bool ignoreMissingPositions_;
    bool ignoreMissingControls_;

    common::TileAggregationMode aggregateTilesFlag_;

    /// \brief Unique index per read if reading the Posisions file failed.
    data::PositionsFile::Record::ClusterCoordinate uniqueFailedReadIndex_;

    std::shared_ptr<data::RawBclBufferGroup> inputBuffer_;
    ThreadSafeQueue<std::shared_ptr<data::RawBclBufferGroup>>& inputBuffersToRecycle_;
    std::shared_ptr<DemuxBuffer> outputBuffer_;
    ThreadSafeQueue<std::shared_ptr<DemuxBuffer>>& outputBufferQueue_;
    std::shared_ptr<data::BclBuffer::PatternedPositionsContainer> patternedFlowcellPositions_;

    GroupSerializer &groupSerializer_;

    static std::atomic<uint32_t> numTaskManagers_;
    static std::condition_variable cvAllTaskManagersDone_;
    static std::mutex mut_;
};


/// \brief BCL loader.
class BclLoader : public Stage
{
public:

    /// \brief Constructor.
    /// \param threadsCount Number of threads.
    /// \param layout Flowcell layout.
    /// \param laneInfo Lane meata data.
    /// \param ignoreMissingBcls Ignore missing BCLs flag.
    /// \param ignoreMissingFilters Ignore missing filters flag.
    /// \param ignoreMissingPositions Ignore missing positions flag.
    BclLoader(
        TaskQueue& taskQueue,
        ThreadSafeQueue<std::shared_ptr<data::RawBclBufferGroup>>& inputBuffersToUse,
        ThreadSafeQueue<std::shared_ptr<data::RawBclBufferGroup>>& inputBuffersToRecycle,
        ThreadSafeQueue<std::shared_ptr<DemuxBuffer>>& outputBuffersToSubmit,
        ThreadSafeQueue<std::shared_ptr<DemuxBuffer>>& outputBuffersToUse,
        const layout::Layout &layout,
        const layout::LaneInfo &laneInfo,
        bool ignoreMissingBcls,
        bool ignoreMissingFilters,
        bool ignoreMissingPositions,
        bool ignoreMissingControls,
        uint32_t numTilesPerBuffer
    );

    virtual std::shared_ptr<BclLoaderTaskManager> getNewTaskManager();

    virtual bool startNewTasks();

private:

    virtual void terminate();

    /// \brief BCL files container type definition.
    typedef std::vector<std::shared_ptr<data::BclFile>> BclFilesContainer;

    /// \brief Cycle BCI files container type definition.
    typedef boost::ptr_vector<data::CycleBCIFile> CycleBciFilesContainer;

    typedef boost::ptr_vector<std::mutex> MutexContainer;

private:

    ThreadSafeQueue<std::shared_ptr<data::RawBclBufferGroup>>& inputBuffersToUse_;
    ThreadSafeQueue<std::shared_ptr<data::RawBclBufferGroup>>& inputBuffersToRecycle_;
    ThreadSafeQueue<std::shared_ptr<DemuxBuffer>>& outputBuffersToSubmit_;
    ThreadSafeQueue<std::shared_ptr<DemuxBuffer>>& outputBuffersToUse_;

    /// \brief Layout.
    const layout::Layout &layout_;

    /// \brief Current lane.
    const layout::LaneInfo &laneInfo_;

    /// \brief Ignore missing BCLs flag.
    bool ignoreMissingBcls_;

    /// \brief Ignore missing filters flag.
    bool ignoreMissingFilters_;

    /// \brief Ignore missing positions flag.
    bool ignoreMissingPositions_;

    /// \brief Ignore missing controls.
    bool ignoreMissingControls_;

    /// \brief Tile being currently processed.
    layout::LaneInfo::TileInfosContainer::const_iterator currentTileInfo_;

    MutexContainer aggregatedBclMutexes_;

    /// \brief Cycle BCI files (one per cycle).
    CycleBciFilesContainer cycleBciFiles_;

    /// \brief Positions for a patterned flowcell.
    std::shared_ptr<data::BclBuffer::PatternedPositionsContainer> patternedFlowcellPositions_;

    /// \brief Cluster counts for the current set of tiles.
    std::vector<common::ClustersCount> clustersCounts_;

    uint32_t numBuffersProcessed_;

    uint32_t numBuffersCreated_;

    const uint32_t numTilesPerBuffer_;

    std::shared_ptr<data::RawBclBufferGroup> prevInputBuffer_;

    bool terminated_;

    TaskManager::GroupSerializer groupSerializer_;
};


} // namespace conversion
} // namespace bcl2fastq


#endif // BCL2FASTQ_CONVERSION_BCLLOADER_HH


