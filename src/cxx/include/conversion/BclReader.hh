/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file BclReader.hh
 *
 * \brief Declaration of BCL reader.
 *
 * \author Aaron Day
 */


#ifndef BCL2FASTQ_CONVERSION_BCLREADER_HH
#define BCL2FASTQ_CONVERSION_BCLREADER_HH


#include "data/BclFileReader.hh"
#include "data/RawBclBuffer.hh"
#include "io/SyncFile.hh"
#include "io/GzipDecompressor.hh"
#include "io/FileBufWithReopen.hh"
#include "layout/Layout.hh"
#include "data/CycleBCIFile.hh"
#include "data/BclFile.hh"
#include "conversion/Stage.hh"
#include "conversion/Task.hh"
#include "conversion/SampleIndex.hh"
#include "common/Types.hh"

#include <vector>
#include <boost/filesystem/path.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <atomic>

namespace bcl2fastq {
namespace conversion {

class BclReaderTaskManager;


class PatternedPositionsContainer : public common::RawDataBuffer
{
public:
    PatternedPositionsContainer() : common::RawDataBuffer(), mut_(), cv_(), isReady_(false) { }

    virtual ~PatternedPositionsContainer() { }

    void setReady() { isReady_ = true; cv_.notify_all(); }

    void waitForReady() { if (!isReady_) { std::unique_lock<std::mutex> lock(mut_); cv_.wait(lock, [this] { return isReady_.load(); }); } }

private:
    std::mutex mut_;
    std::condition_variable cv_;
    std::atomic<bool> isReady_;
};

/// \brief Task: Read BCL data into buffer.
class BclReadTask : public Task
{
public:

    /// \brief Constructor.
    /// \param taskManager Manages all tasks.
    /// \param bclFile BCL file to be loaded.
    /// \param outputBuffer Buffer for data to be loaded to.
    BclReadTask(
        std::shared_ptr<BclReaderTaskManager>& taskManager,
        data::BclFileReader& bclFileReader,
        data::RawBclBufferGroup& outputBuffer
    );

public:

    virtual bool execute(int32_t threadNum);

private:

    /// \brief BCL file reader.
    data::BclFileReader &bclFileReader_;

    /// \brief Tile meta data.
    layout::LaneInfo::TileInfosContainer::const_iterator tileInfoIter_;

    /// \brief Buffer to load data to.
    data::RawBclBufferGroup &outputBuffer_;
};


/// \brief Task: Read positions data into buffer.
class PositionsReadTask : public Task
{
public:

    /// \brief Constructor.
    /// \param positionsFile Positions file.
    /// \param intensitiesDir Path to intensities directory.
    /// \param aggregateTilesFlag All tiles aggregated into single file flag.
    /// \param laneInfo Lane meata data.
    /// \param tileInfo Tile meta data.
    /// \param ignoreMissingPositions Ignore missing positions flag.
    /// \param outputBuffer Buffer for data to be loaded to.
    PositionsReadTask(
        std::shared_ptr<BclReaderTaskManager>& taskManager,
        std::shared_ptr<io::SyncFile> &positionsFile,
        const boost::filesystem::path &intensitiesDir,
        bool aggregateTilesFlag,
        bool isPatternedFlowcell,
        const layout::LaneInfo &laneInfo,
        layout::LaneInfo::TileInfosContainer::const_iterator &tileInfo,
        bool ignoreMissingPositions,
        data::RawBclBufferGroup& outputBuffer,
        std::shared_ptr<PatternedPositionsContainer>& patternedFlowcellPositions
    );

public:

    virtual bool execute(int32_t threadNum);

private:

    bool readEntirePosFile(common::RawDataBuffer& outputBuffer);

    bool readPartOfAggregatedPosFile(common::RawDataBuffer& outputBuffer);

    /// \brief Positions file.
    std::shared_ptr<io::SyncFile>& positionsFile_;

    /// \brief Intensities directory path.
    const boost::filesystem::path &intensitiesDir_;

    /// \brief All tiles aggregated into single file flag.
    bool aggregateTilesFlag_;

    /// \brief Single location file flag.
    bool isPatternedFlowcell_;

    /// \brief Lane meta data.
    const layout::LaneInfo &laneInfo_;

    /// \brief Tile meta data.
    layout::LaneInfo::TileInfosContainer::const_iterator tileInfoIter_;

    /// \brief Ignore missing positions flag.
    bool ignoreMissingPositions_;

    /// \brief Buffer to load data to.
    data::RawBclBufferGroup &outputBuffer_;

    /// \brief Buffer for patterned flowcell positions
    std::shared_ptr<PatternedPositionsContainer> patternedFlowcellPositions_;
};


/// \brief Task: Read filter data into buffer.
class FilterReadTask : public Task
{
public:

    /// \brief Constructor.
    /// \param filterFile Filter file.
    /// \param inputDir Path to input directory.
    /// \param aggregateTilesFlag All tiles aggregated into single file flag.
    /// \param laneInfo Lane meata data.
    /// \param tileInfo Tile meta data.
    /// \param ignoreMissingFilters Ignore missing filters flag.
    /// \param outputBuffer Buffer for data to be loaded to.
    /// \param controlsOutputBuffer Buffer for controls data to be loaded to.
    FilterReadTask(
        std::shared_ptr<BclReaderTaskManager>& taskManager,
        std::shared_ptr<io::SyncFile> &filterFile,
        const boost::filesystem::path &inputDir,
        common::TileAggregationMode tileAggregationMode,
        const layout::LaneInfo &laneInfo,
        layout::LaneInfo::TileInfosContainer::const_iterator& tileInfo,
        bool ignoreMissingFilters,
        data::RawBclBufferGroup& outputBuffer
    );

public:

    virtual bool execute(int32_t threadNum);

private:

    bool readEntireFilterFile(data::RawBclBuffer& rawBclBuffer);

    bool readPartOfAggregatedFilterFile(common::RawDataBuffer& outputBuffer);

    /// \brief Validate the filter records
    /// \param recordsRead Number of records read from file.
    /// \param filePath Path to file.
    void validateFilterRecords(std::size_t                    recordsRead,
                               const boost::filesystem::path& filePath,
                               data::RawBclBuffer&            outputBuffer);

    /// \brief Filter file.
    std::shared_ptr<io::SyncFile>& filterFile_;

    /// \brief Input directory path.
    const boost::filesystem::path &inputDir_;

    /// \brief File aggregation type.
    common::TileAggregationMode tileAggregationMode_;
    /// \brief All tiles aggregated into single file flag.
    bool aggregateTilesFlag_;


    /// \brief Lane meta data.
    const layout::LaneInfo &laneInfo_;

    /// \brief Tile meta data.
    layout::LaneInfo::TileInfosContainer::const_iterator tileInfoIter_;

    /// \brief Ignore missing filters flag.
    bool ignoreMissingFilters_;

    /// \brief Buffer to load control data to.
    data::RawBclBufferGroup& outputBuffer_;
};


class ControlReadTask : public Task
{
public:
    ControlReadTask(std::shared_ptr<BclReaderTaskManager>& taskManager,
        const boost::filesystem::path &inputDir,
        const layout::LaneInfo &laneInfo,
        layout::LaneInfo::TileInfosContainer::const_iterator &tileInfo,
        bool ignoreMissingControls,
        data::RawBclBufferGroup& outputBuffer);

    virtual bool execute(int32_t threadNum);

private:
    bool readEntireControlFile(common::RawDataBuffer& outputBuffer);

    /// \brief Input directory path.
    const boost::filesystem::path &inputDir_;

    /// \brief Lane meta data.
    const layout::LaneInfo &laneInfo_;

    /// \brief Tile meta data.
    layout::LaneInfo::TileInfosContainer::const_iterator tileInfoIter_;

    /// \brief Ignore missing controls flag.
    bool ignoreMissingControls_;

    /// \brief Buffer to load data to.
    data::RawBclBufferGroup& outputBuffer_;
};

class BclReaderTaskManager : public TaskManager
{
public:
    BclReaderTaskManager(bool ignoreMissingBcls,
                         bool ignoreMissingFilters,
                         bool ignoreMissingPositions,
                         bool ignoreMissingControls,
                         std::shared_ptr<data::RawBclBufferGroup>& outputBuffer,
                         ThreadSafeQueue<std::shared_ptr<data::RawBclBufferGroup>>& outputBufferQueue,
                         std::shared_ptr<PatternedPositionsContainer>& patternedFlowcellPositions,
                         TaskManager::GroupSerializer &groupSerializer);

    virtual ~BclReaderTaskManager();

    data::RawBclBufferGroup& getOutputBuffer() { return *outputBuffer_; }

    static void waitForAllTasksToFinish();

private:
    void postExecute();

    bool ignoreMissingBcls_;
    bool ignoreMissingFilters_;
    bool ignoreMissingPositions_;
    bool ignoreMissingControls_;

    std::shared_ptr<data::RawBclBufferGroup> outputBuffer_;
    ThreadSafeQueue<std::shared_ptr<data::RawBclBufferGroup>>& outputBufferQueue_;
    std::shared_ptr<PatternedPositionsContainer> patternedFlowcellPositions_;

    static std::atomic<uint32_t> numTaskManagers_;
    static std::condition_variable cvAllTaskManagersDone_;
    static std::mutex mut_;

    TaskManager::GroupSerializer &groupSerializer_;
};


/// \brief BCL loader.
class BclReader : public Stage
{
public:

    /// \brief Constructor.
    /// \param threadsCount Number of threads.
    /// \param layout Flowcell layout.
    /// \param laneInfo Lane meata data.
    /// \param ignoreMissingBcls Ignore missing BCLs flag.
    /// \param ignoreMissingFilters Ignore missing filters flag.
    /// \param ignoreMissingPositions Ignore missing positions flag.
    /// \param inputDir Input directory.
    /// \param intensitiesDir Intensities directory.
    BclReader(
        TaskQueue& taskQueue,
        ThreadSafeQueue<std::shared_ptr<data::RawBclBufferGroup>>& outputBuffersToSubmit,
        ThreadSafeQueue<std::shared_ptr<data::RawBclBufferGroup>>& outputBuffersToUse,
        const layout::Layout &layout,
        const layout::LaneInfo &laneInfo,
        bool ignoreMissingBcls,
        bool ignoreMissingFilters,
        bool ignoreMissingPositions,
        bool ignoreMissingControls,
        const boost::filesystem::path &inputDir,
        const boost::filesystem::path &intensitiesDir
    );

uint32_t getNumTilesPerBuffer() const { return numTilesPerBuffer_; }

private:

    virtual void terminate();

    void calculateNumTilesPerBuffer();

    virtual std::shared_ptr<BclReaderTaskManager> getNewTaskManager();

    virtual bool startNewTasks();

    /// \brief BCL files container type definition.
    typedef boost::ptr_vector<io::SyncFile> BclFilesContainer;

    typedef boost::ptr_vector<data::BclFileReader> BclFileReadersContainer;

    /// \brief Cycle BCI files container type definition.
    typedef boost::ptr_vector<data::CycleBCIFile> CycleBciFilesContainer;

private:

    ThreadSafeQueue<std::shared_ptr<data::RawBclBufferGroup>>& outputBuffersToSubmit_;
    ThreadSafeQueue<std::shared_ptr<data::RawBclBufferGroup>>& outputBuffersToUse_;

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

    /// \brief Input directory.
    boost::filesystem::path inputDir_;

    /// \brief Intensities directory.
    boost::filesystem::path intensitiesDir_;

    BclFileReadersContainer bclFileReaders_;

    /// \brief Positions file.
    std::shared_ptr<io::SyncFile> positionsFile_;

    /// \brief Filter file.
    std::shared_ptr<io::SyncFile> filterFile_;

    /// \brief Positions for a patterned flowcell.
    std::shared_ptr<PatternedPositionsContainer> patternedFlowcellPositions_;

    uint32_t numTilesPerBuffer_;

    uint32_t numBuffersCreated_;

    uint32_t groupNumber_;

    TaskManager::GroupSerializer groupSerializer_;
}; // end class BclReader

} // namespace conversion
} // namespace bcl2fastq

#endif // BCL2FASTQ_CONVERSION_BCLREADER_HH

