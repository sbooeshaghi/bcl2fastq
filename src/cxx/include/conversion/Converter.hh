/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file Converter.hh
 *
 * \brief Declaration of BCL to FASTQ converter.
 *
 * \author Marek Balint
 * \author Mauricio Varea
 */


#ifndef BCL2FASTQ_CONVERSION_CONVERTER_HH
#define BCL2FASTQ_CONVERSION_CONVERTER_HH


#include <vector>
#include <map>

//#include <boost/unordered_map.hpp>
#include <boost/noncopyable.hpp>

#include "data/InteropFile.hh"
#include "common/Threads.hh"
#include "config/Bcl2FastqOptions.hh"
#include "layout/Layout.hh"
#include "stats/TileStats.hpp"
#include "stats/BarcodeStats.hh"
#include "conversion/FastqBuffer.hh"
#include "conversion/BclReader.hh"
#include "conversion/BclLoader.hh"
#include "conversion/Demultiplexer.hh"
#include "conversion/FastqCreator.hh"
#include "conversion/FastqWriter.hh"
#include "conversion/ThreadPool.hh"


namespace bcl2fastq
{

// Forward declarations
namespace stats
{
class ConversionStatsXml;
class JsonObject;
class JsonArray;
class JsonNumber;
}

namespace conversion
{

const uint32_t MAX_BCL_BUFFER_QUEUE_DEPTH = 5;

namespace detail
{

class FastqSummaryStats
{
public:
    FastqSummaryStats(common::SampleNumber  sampleNumber,
                      const std::string&    sampleName,
                      common::TileNumber    tileNumber,
                      common::TotalClustersCount numberReadsRaw,
                      common::TotalClustersCount numberReadsPF);

    bool operator<(const FastqSummaryStats& other) const;

    common::SampleNumber  sampleNumber_;
    std::string           sampleName_;
    common::TileNumber    tileNumber_;
    common::TotalClustersCount numberReadsRaw_;
    common::TotalClustersCount numberReadsPF_;
};

typedef std::vector<FastqSummaryStats> FastqSummaryStatsContainer;

}

/// \todo: consider managing interOp stats from here as well
class ConverterStats
{
    typedef std::pair<layout::LaneInfo,stats::BarcodeHits> UnknownPair;
    typedef std::map<layout::LaneInfo,stats::BarcodeHits> UnknownBarcodeHits;

    typedef std::pair<layout::LaneInfo,DemultiplexTask::DemuxStats> DemuxPair;
    typedef std::map<layout::LaneInfo,DemultiplexTask::DemuxStats> LaneBarcodeStats;

    typedef std::map<layout::LaneInfo,ConversionStatsForSampleTileBarcode> LaneTileStats;

    typedef std::map<std::string, stats::BarcodeStats> SampleLaneBarcodeStats;
    typedef std::map<std::string, SampleLaneBarcodeStats> ProjectSampleLaneBarcodeStats;
    typedef std::map<std::string, ProjectSampleLaneBarcodeStats> FlowcellProjectSampleLaneBarcodeStats;

    typedef std::map<common::ReadNumber, stats::ReadBarcodeStats > AnnotatedTileReadBarcodeStats;
    typedef std::pair<stats::TileBarcodeStats, AnnotatedTileReadBarcodeStats > TileReadStats;
    typedef std::map<common::TileNumber, TileReadStats > TileReadBarcodeStats;
    typedef std::map<std::string, TileReadBarcodeStats> BarcodeTileReadBarcodeStats;
    typedef std::map<std::string, BarcodeTileReadBarcodeStats> SampleBarcodeTileReadBarcodeStats;
    typedef std::map<std::string, SampleBarcodeTileReadBarcodeStats> ProjectSampleBarcodeTileReadBarcodeStats;
    typedef std::map<std::string, ProjectSampleBarcodeTileReadBarcodeStats> FlowcellProjectSampleBarcodeTileReadBarcodeStats;

    typedef std::map<std::string, stats::BarcodeStats> BarcodeStatsMap;
    typedef std::map<int, BarcodeStatsMap> BarcodeLaneStatsMap;

    // map the sampleId to sample index and the stats object
    typedef std::map<std::string, std::pair<size_t, stats::TileBarcodeStats>> SampleReadBarcodeStats;

    typedef std::map<std::string, SampleReadBarcodeStats> ProjectSampleReadBarcodeStats;

public:
    ConverterStats( const layout::Layout &layout,
                    const boost::filesystem::path &fastqSummaryPartialFileName,
                    const boost::filesystem::path &demuxSummaryPartialFileName,
                    const boost::filesystem::path &adapterTrimmingTxtFile,
                    const boost::filesystem::path &interopFile,
                    const boost::filesystem::path &demuxStatsXmlFile,
                    const boost::filesystem::path &conversionStatsXmlFile,
                    const boost::filesystem::path &statsJsonFile);

    boost::filesystem::path getDemuxStatsXmlPath() const {return boost::filesystem::canonical( demuxStatsXmlFile_ );}
    boost::filesystem::path getConversionStatsXmlPath() const {return boost::filesystem::canonical( conversionStatsXmlFile_ );}

    DemultiplexTask::DemuxStats &getDemuxStats(const layout::LaneInfo &laneInfo)            {return demuxStats_[laneInfo];}
    ConversionStatsForSampleTileBarcode &getConversionStats(const layout::LaneInfo &laneInfo)  {return conversionStats_[laneInfo];}
    stats::BarcodeHits &getUnknownBarcodeStat(const layout::LaneInfo &laneInfo)             {return topUnknownBarcodes_[laneInfo];}

    void dumpInteropStats() const;
    void dumpDemuxStats();
    void dumpConversionStats();

private:
    void recordBarcodeStats(const BarcodeStatsMap& barcodeStatsMap,
                            const std::string& barcodeName,
                            stats::JsonArray& jsonIndexMetrics);

    void recordReadStats(stats::ConversionStatsXml& conversionStatsXml,
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
                         stats::JsonArray* jsonReadMetrics);

    void recordUnknownBarcodes(stats::JsonObject& jsonStats,
                               stats::ConversionStatsXml& conversionStatsXml) const;

    void writeAdapterStats(const layout::LaneInfo&              laneInfo,
                           const ProjectSampleReadBarcodeStats& stats,
                           std::ostream&                        adapterOstr,
                           std::ostream&                        trimLengthOstr,
                           std::map<std::pair<std::string, std::string>, std::vector<size_t>>& trimmedBaseMap) const;

    void getSampleNameAndNumber(const layout::LaneInfo& laneInfo,
                                const std::string&      project,
                                const std::string&      sampleId,
                                std::string&            sampleName,
                                common::SampleNumber&   sampleNumber) const;

    void writeSummaryStats(const ProjectSampleBarcodeTileReadBarcodeStats& stats,
                           const layout::LaneInfo&                         laneInfo) const;

    void writeFastqSummaryStats(common::LaneNumber                          laneNumber,
                                const detail::FastqSummaryStatsContainer& fastqSummaryStats) const;

    void writeDemuxSummaryStats(common::LaneNumber                        laneNumber,
                                const detail::FastqSummaryStatsContainer& fastqSummaryStats) const;

    void writeAdapterStats(common::LaneNumber                                   laneNumber,
                           const std::string&                                   project,
                           const std::string&                                   sampleId,
                           const std::string&                                   sampleName,
                           common::SampleNumber                                 sampleNumber,
                           const stats::AllReadsStats::TrimLengthCountForReads& trimLengthCountForReads,
                           std::ostream&                                        adapterOstr,
                           std::ostream&                                        trimLengthOstr,
                           std::map<std::pair<std::string, std::string>, std::vector<size_t>>& trimmedBaseMap) const;

    void deserializeMetadata(const layout::LaneInfo &laneInfo,
                             size_t sampleIndex,
                             size_t barcodeIndex,
                             std::string &indexName,
                             std::string &sampleId,
                             std::string &projectName,
                             std::string &sampleName,
                             common::SampleNumber &sampleNumber) const;

    const layout::Layout &layout_;

    /// \brief Path to FastqSummaryF1* file.
    const boost::filesystem::path fastqSummaryPartialFileName_;

    /// \brief Path to DemuxSummaryF1* file.
    const boost::filesystem::path demuxSummaryPartialFileName_;

    /// \brief Path to AdapterTrimming.txt file
    const boost::filesystem::path adapterTrimmingTxtFile_;

    /// \brief Path to IndexMetricsOut.bin file
    const boost::filesystem::path interopFile_;

    /// \brief Path to DemultiplexingStats.xml file
    const boost::filesystem::path demuxStatsXmlFile_;

    /// \brief Path to ConversionStats.xml file
    const boost::filesystem::path conversionStatsXmlFile_;

    /// \brief Path to Stats.json file.
    const boost::filesystem::path statsJsonFile_;

    /// \brief Demultiplexing statistics (one per lane).
    LaneBarcodeStats demuxStats_;

    /// \brief Basecalling statistics (one per lane).
    LaneTileStats conversionStats_;

    /// \brief Table of unknown barcodes (one per lane).
    UnknownBarcodeHits topUnknownBarcodes_;

    std::vector<stats::AllReadsStats::TrimLengthCountForReads> trimLengthCountForReadsPerLane_;

    /// \brief Map the index name to the BarcodeStats object
    BarcodeLaneStatsMap barcodeLaneStatsMap_;
};


/// \brief BCL to FASTQ converter.
///
/// \dot
/// digraph LAYOUT {
///     fontsize=12;
///
///     source [ label="BCLs" shape=folder ];
///     stage1 [ label="load BCLs" URL="\ref BclLoader" shape=component ];
///     buffer1 [ label="BCL Buffer" URL="\ref BclBuffer" shape=box3d ];
///     stage2 [ label="demultiplex" URL="\ref Demultiplexer" shape=component ];
///     buffer2 [ label="BCL Buffer" URL="\ref BclBuffer" shape=box3d ];
///     stage3 [ label="create FASTQs" URL="\ref FastqCreator" shape=component ];
///     buffer3 [ label="FASTQ Buffer" URL="\ref FastqBuffer" shape=box3d ];
///     stage4 [ label="write FASTQs" URL="\ref FastqWriter" shape=component ];
///     sink [ label="FASTQs" shape=folder ];
///
///     {
///         rank=same;
///         source -> stage1;
///     }
///     {
///         rank=same;
///         stage1 -> buffer1;
///     }
///     {
///         rank=same;
///         buffer1 -> stage2;
///     }
///     {
///         rank=same;
///         stage2 -> buffer2;
///     }
///     {
///         rank=same;
///         buffer2 -> stage3;
///     }
///     {
///         rank=same;
///         stage3 -> buffer3;
///     }
///     {
///         rank=same;
///         buffer3 -> stage4;
///     }
///     {
///         rank=same;
///         stage4 -> sink;
///     }
///
/// \enddot
class Converter : private boost::noncopyable
{
public:
    /// \brief Constructor.
    /// \param options Program options.
    /// \param layout Flowcell layout.
    /// \param laneInfo Lane metadata.
    Converter(
        const bcl2fastq::config::Bcl2FastqOptions &options,
        const layout::Layout &layout,
        const layout::LaneInfo &laneInfo,
        bcl2fastq::conversion::ConverterStats &allStats
    );

public:

    static std::unique_ptr<Stage> createNewFastqCreator(TaskQueue& processTaskQueue,
                                                        ThreadSafeQueue<std::shared_ptr<DemuxBuffer>>& inUseDemuxBufferQueue,
                                                        ThreadSafeQueue<std::shared_ptr<DemuxBuffer>>& recycledBclBufferQueue,
                                                        ThreadSafeMap<uint32_t, std::shared_ptr<FastqBuffer>>& inUseFastqBufferMap,
                                                        ThreadSafeQueue<std::shared_ptr<FastqBuffer>>& recycledFastqBufferQueue,
                                                        const layout::Layout& layout,
                                                        const layout::LaneInfo& laneInfo,
                                                        const bcl2fastq::config::Bcl2FastqOptions& options,
                                                        bcl2fastq::conversion::ConverterStats& allStats);

    /// \brief Run conversion.
    void run();

private:

    /// \brief Stage references container type definition.
    typedef std::vector<boost::reference_wrapper<Stage> > StageRefsContainer;

    /// \brief Thread function executing given stage.
    /// \param stages Individual stages to be executed.
    /// \param threadNum Thread number.
    void stageThreadFunction(StageRefsContainer &stages, common::ThreadVector::size_type threadNum);

private:

    class StatsFinalizer : private boost::noncopyable
    {
    public:
        StatsFinalizer(bool isOnly1Sample,
                       const ConversionStatsForSampleTileBarcode& conversionStats,
                       DemultiplexTask::DemuxStats& demuxStats);

        ~StatsFinalizer();

    private:
        DemultiplexTask::DemuxStats& demuxStats_;
        const ConversionStatsForSampleTileBarcode& conversionStats_;
        bool isOnly1Sample_;
    };

    /// \brief Flowcell layout.
    const layout::Layout &layout_;

    /// \brief lane info
    const layout::LaneInfo& laneInfo_;

    /// \brief Updates stats upon destruction.
    StatsFinalizer statsFinalizer_;

    /// \brief Stats
    ConverterStats allStats_;

    ThreadSafeQueue<std::shared_ptr<data::RawBclBufferGroup>> inUseRawBclBufferQueue_;
    ThreadSafeQueue<std::shared_ptr<data::RawBclBufferGroup>> recycledRawBclBufferQueue_;

    ThreadSafeQueue<std::shared_ptr<DemuxBuffer>> inUseBclBufferQueue_;
    ThreadSafeQueue<std::shared_ptr<DemuxBuffer>> recycledBclBufferQueue_;

    ThreadSafeQueue<std::shared_ptr<DemuxBuffer>> inUseDemuxBufferQueue_;

    ThreadSafeMap<uint32_t, std::shared_ptr<FastqBuffer>> inUseFastqBufferMap_;
    ThreadSafeQueue<std::shared_ptr<FastqBuffer>> recycledFastqBufferQueue_;

    TaskQueue processTaskQueue_;
    TaskQueue ioReadTaskQueue_;
    TaskQueue ioWriteTaskQueue_;

    Terminator terminator_;

    /// \brief Barcode to sample translation table.
    const layout::BarcodeTranslationTable barcodeTranslationTable_;

    /// \brief BCL reader.
    BclReader bclReader_;

    /// \brief BCL loader.
    BclLoader bclLoader_;

    /// \brief Demultiplexer.
    Demultiplexer demultiplexer_;

    /// \brief FASTQ creator.
    std::unique_ptr<Stage> fastqCreator_;

    /// \brief FASTQ writer.
    FastqWriter fastqWriter_;

    ThreadPool processThreadPool_;
    ThreadPool ioReadThreadPool_;
    ThreadPool ioWriteThreadPool_;
};


} // namespace conversion
} // namespace bcl2fastq


#endif // BCL2FASTQ_CONVERSION_CONVERTER_HH


