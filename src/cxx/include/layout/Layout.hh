/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file Layout.hh
 *
 * \brief Declaration of data layout.
 *
 * \author Marek Balint
 * \author Mauricio Varea
 */


#ifndef BCL2FASTQ_LAYOUT_LAYOUT_HH
#define BCL2FASTQ_LAYOUT_LAYOUT_HH


#include <vector>

#include <boost/filesystem/path.hpp>
#include <boost/optional.hpp>
#include <boost/noncopyable.hpp>
#include <boost/regex.hpp>

#include "layout/FlowcellInfo.hh"
#include "layout/LaneInfo.hh"
#include "layout/ReadInfo.hh"
#include "common/Types.hh"
#include "config/Bcl2FastqOptions.hh"

namespace bcl2fastq {

namespace config {
class SampleSheetCsv;
}

namespace layout {



/// \brief Data %layout.
///
/// \dot
/// digraph LAYOUT {
///     nodesep=0.5;
///     rankdir=LR;
///     fontsize=12;
///
///     layout [ label="layout" URL="\ref Layout" shape=box ];
///
///         flowcell [ label="flowcell" URL="\ref FlowcellInfo" shape=box ];
///         layout -> flowcell [ label="1:1" URL="\ref Layout::flowcellInfo_" ];
///
///             flowcellId [ label="ID" URL="\ref FlowcellInfo::flowcellId_" shape=ellipse ];
///             flowcell -> flowcellId
///
///             runNumber [ label="run number" URL="\ref FlowcellInfo::runNumber_" shape=ellipse ];
///             flowcell -> runNumber
///
///             instrument [ label="instrument" URL="\ref FlowcellInfo::instrument_" shape=ellipse ];
///             flowcell -> instrument
///
///             aggregateTilesFlag [ label="aggregated tiles flag" URL="\ref FlowcellInfo::aggregateTilesFlag_" shape=ellipse ];
///             flowcell -> aggregateTilesFlag
///
///         lane [ label="lane" URL="\ref LaneInfo" shape=box ];
///         layout -> lane [ label="1:n" URL="\ref Layout::laneInfos_" ];
///
///             laneNumber [ label="number" URL="\ref LaneInfo::number_" shape=ellipse ];
///             lane -> laneNumber
///
///             sample [ label="sample" URL="\ref SampleInfo" shape=box ];
///             lane -> sample [ label="1:n" URL="\ref LaneInfo::sampleInfos_" ]
///
///                 sampleNumber [ label="number" URL="\ref SampleInfo::number_" shape=ellipse ];
///                 sample -> sampleNumber
///
///                 sampleId [ label="ID" URL="\ref SampleInfo::sampleId_" shape=ellipse ];
///                 sample -> sampleId
///
///                 sampleName [ label="name" URL="\ref SampleInfo::sampleName_" shape=ellipse ];
///                 sample -> sampleName
///
///                 barcode [ label="barcode" URL="\ref Barcode" shape=box ];
///                 sample -> barcode [ label="1:n" URL="\ref SampleInfo::barcodes_" ];
///
///                     component [ label="component" URL="\ref Barcode::Component" shape=box ];
///                     barcode -> component [ label="1:n" URL="\ref Barcode::components_" ];
///
///                         componentLength [ label="length" URL="\ref Barcode::Component::length_" shape=ellipse ];
///                         component -> componentLength;
///
///                         componentValue [ label="value" URL="\ref Barcode::Component::value_" shape=ellipse ];
///                         component -> componentValue;
///
///             tile [ label="tile" URL="\ref TileInfo" shape=box ];
///             lane -> tile [ label="1:n" URL="\ref LaneInfo::tileInfos_" ]
///
///             tileNumber [ label="number" URL="\ref TileInfo::number_" shape=ellipse ];
///             tile -> tileNumber
///
///             tileIndex [ label="index" URL="\ref TileInfo::index_" shape=ellipse ];
///             tile -> tileIndex
///
///             tileClustersCount [ label="clusters count" URL="\ref TileInfo::clustersCount_" shape=ellipse ];
///             tile -> tileClustersCount
///
///             tileHaveClustersCount [ label="have clusters count flag" URL="\ref TileInfo::haveClustersCount_" shape=ellipse ];
///             tile -> tileHaveClustersCount
///
///             tileSkippedTilesCount [ label="skipped tiles count" URL="\ref TileInfo::skippedTiles" shape=ellipse ];
///             tile -> tileSkippedTilesCount
///
///             tileSkippedClustersCount [ label="skipped clusters count" URL="\ref TileInfo::skippedClusters_" shape=ellipse ];
///             tile -> tileSkippedClustersCount
///
///         read [ label="read" URL="\ref ReadInfo" shape=box ];
///         layout -> read [ label="1:n" URL="\ref Layout::readInfos_" ];
///
///             readNumber [ label="number" URL="\ref ReadInfo::number_" shape=ellipse ];
///             read -> readNumber
///
///             indexReadFlag [ label="index read flag" URL="\ref ReadInfo::indexReadFlag_" shape=ellipse ];
///             read -> indexReadFlag
///
///             maskAdapters [ label="mask adapters" URL="\ref ReadInfo::maskAdapters_" shape=ellipse ];
///             read -> maskAdapters
///
///             trimAdapters [ label="trim adapters" URL="\ref ReadInfo::trimAdapters_" shape=ellipse ];
///             read -> trimAdapters
/// }
/// \enddot
///
/// - %layout
///     - flowcell
///         - ID <- RunInfo.xml
///         - run number <- RunInfo.xml
///         - instrument <- RunInfo.xml
///         - aggregated tiles flag <- aggregated-tiles cmdline argument
///     - lanes
///         - number <- RunInfo.xml (1..n)
///         - samples
///             - number <- SampleSheet.csv (0,1..n)
///             - ID <- SampleSheet.csv
///             - name <- SampleSheet.csv
///             - barcodes <- SampleSheet.csv
///                 - barcode components
///         - tiles
///             - number <- BCI file (Nova), config.xml (HiSeq/MiSeq)
///             - index <- index in tiles container (0..n)
///             - clusters count <- BCI file (Nova), BCL file (HiSeq/MiSeq)
///             - have clusters count flag <- true only in BCI file is present
///             - number of skipped tiles immediately before this tile <- tiles cmdline argument
///             - number of skipped clusters immediately before this tile <- tiles cmdline argument
///     - reads
///         - number <- RunInfo.xml (1..n)
///         - index read flag <- RunInfo.xml
///         - mask adapters <- SampleSheet.csv
///         - trim adapters <- SampleSheet.csv
///         - cycles
///             - number <- RunInfo.xml (FirstCycle..LastCycle or NumCycles)
class Layout : private boost::noncopyable
{
public:
/// \brief Detect layout.
/// \param inputDir Path to input directory.
/// \param reportsDir Path to reports directory.
/// \param statsDir Path to stats directory.
/// \param tilesAggregationFlag All tiles in same BCL file vs. separate BCL file for each tile.
/// \param tilesFilterList List of regexps for tile filtering.
/// \return Initialized layout object.
    Layout(const boost::filesystem::path& intensitiesDir,
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
           config::BarcodeMismatchCountsContainer& maxMismatches,
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
           bool includeNonPfClusters);

    /// \brief Get flowcell metadata.
    /// \return Flowcell metadata.
    const FlowcellInfo & getFlowcellInfo() const;

    /// \brief Get beginning of lanes.
    /// \return Iterator to beginning of lanes.
    LaneInfosContainer::const_iterator laneInfosBegin() const;

    /// \brief Get end of lanes.
    /// \return Iterator to end of lanes.
    LaneInfosContainer::const_iterator laneInfosEnd() const;

    /// \brief Get the map of tile number to file name.
    /// \return Map of tile number to file name.
    const common::TileFileMap& getTileFileMap() const { return tileFileMap_; }

    common::NumBasesPerByte getNumBasesPerByte() const { return numBasesPerByte_; }

private:

    void detectFlowcellInfo(const boost::filesystem::path& intensitiesDir,
                            const RunInfoXml&              runInfoXml,
                            bool                           hasBci,
                            bool                           hasCbcl);

    /// \brief Flowcell metadata.
    FlowcellInfo flowcellInfo_;

    /// \brief Lanes metadata container.
    LaneInfosContainer laneInfos_;

    /// \brief Map of tile number to file name.
    common::TileFileMap tileFileMap_;

    /// \brief Number of bases per byte.
    common::NumBasesPerByte numBasesPerByte_;
};

class TileLayoutDetector
{
public:

    static void detectTileLayout(
        const std::vector<std::string>& tilesFilterList,
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
        bool                            ignoreMissingBcls);

protected:

	static void readTilesFromCbclHeaders(
        const boost::filesystem::path&   inputDir,
        const config::SampleSheetCsv&    sampleSheet,
        const std::vector<boost::regex>& tileRegexps,
        const std::string&               flowcellId,
        std::vector<LaneInfo>&           lanes,
        common::TileFileMap&             tileFileMap,
        common::NumBasesPerByte&         numBasesPerByte,
        common::CycleNumber              read1StartCycle,
        common::CycleNumber              read1EndCycle,
        bool                             ignoreMissingBcls);

}; // end class TileLayoutDetector

// TODO: Refactor this into a class.
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
                      const std::vector<std::vector<size_t>>& barcodeLengthsForLane);

} // namespace layout


} // namespace bcl2fastq


#endif // BCL2FASTQ_LAYOUT_LAYOUT_HH


// http://www.graphviz.org/content/generate-directory-tree-dot
// digraph LAYOUT {
//     nodesep=0.5;
//     rankdir=LR;
//     fixedsize=true;
//     node [ concentrate=true ];
//
//     layout [ label="layout" URL="\ref Layout" shape=box ];
//
//         flowcell_point [ shape = point ];
//         flowcell [ label="flowcell" URL="\ref FlowcellInfo" shape=box ];
//         flowcell_point -> flowcell [ label="1:1" URL="\ref Layout::flowcellInfo_" ];
//
//             flowcellId_point [ shape = point ];
//             flowcellId [ label="ID" URL="\ref FlowcellInfo::flowcellId_" shape=ellipse ];
//             flowcellId_point -> flowcellId;
//
//             runNumber_point [ shape = point ];
//             runNumber [ label="run number" URL="\ref FlowcellInfo::runNumber_" shape=ellipse ];
//             runNumber_point -> runNumber;
//
//             instrument_point [ shape = point ];
//             instrument [ label="instrument" URL="\ref FlowcellInfo::instrument_" shape=ellipse ];
//             instrument_point -> instrument;
//
//             aggregateTilesFlag_point [ shape = point ];
//             aggregateTilesFlag [ label="aggregated tiles flag" URL="\ref FlowcellInfo::aggregateTilesFlag_" shape=ellipse ];
//             aggregateTilesFlag_point -> aggregateTilesFlag;
//
//         lane_point [ shape = point ];
//         lane [ label="lane" URL="\ref LaneInfo" shape=box ];
//         lane_point -> lane [ label="1:n" URL="\ref Layout::laneInfos_" ];
//
//             laneNumber_point [ shape = point ];
//             laneNumber [ label="number" URL="\ref LaneInfo::number_" shape=ellipse ];
//             laneNumber_point -> laneNumber;
//
//             sample_point [ shape = point ];
//             sample [ label="sample" URL="\ref SampleInfo" shape=box ];
//             sample_point -> sample [ label="1:n" URL="\ref LaneInfo::sampleInfos_" ];
//
//                 sampleNumber_point [ shape = point ];
//                 sampleNumber [ label="number" URL="\ref SampleInfo::number_" shape=ellipse ];
//                 sampleNumber_point -> sampleNumber;
//
//                 sampleId_point [ shape = point ];
//                 sampleId [ label="ID" URL="\ref SampleInfo::sampleId_" shape=ellipse ];
//                 sampleId_point -> sampleId;
//
//                 sampleName_point [ shape = point ];
//                 sampleName [ label="name" URL="\ref SampleInfo::sampleName_" shape=ellipse ];
//                 sampleName_point -> sampleName;
//
//                 barcode_point [ shape = point ];
//                 barcode [ label="barcode" URL="\ref Barcode" shape=box ];
//                 barcode_point -> barcode [ label="1:n" URL="\ref SampleInfo::barcodes_" ];
//
//                     component_point [ shape = point ];
//                     component [ label="component" URL="\ref Barcode::Component" shape=box ];
//                     component_point -> component [ label="1:n" URL="\ref Barcode::components_" ];
//
//                         componentLength_point [ shape = point ];
//                         componentLength [ label="length" URL="\ref Barcode::Component::length_" shape=ellipse ];
//                         componentLength_point -> componentLength;
//
//                         componentValue_point [ shape = point ];
//                         componentValue [ label="value" URL="\ref Barcode::Component::value_" shape=ellipse ];
//                         componentValue_point -> componentValue;
//
//             tile_point [ shape = point ];
//             tile [ label="tile" URL="\ref TileInfo" shape=box ];
//             tile_point -> tile [ label="1:n" URL="\ref LaneInfo::tileInfos_" ];
//
//             tileNumber_point [ shape = point ];
//             tileNumber [ label="number" URL="\ref TileInfo::number_" shape=ellipse ];
//             tileNumber_point -> tileNumber;
//
//             tileIndex_point [ shape = point ];
//             tileIndex [ label="index" URL="\ref TileInfo::index_" shape=ellipse ];
//             tileIndex_point -> tileIndex;
//
//             tileClustersCount_point [ shape = point ];
//             tileClustersCount [ label="clusters count" URL="\ref TileInfo::clustersCount_" shape=ellipse ];
//             tileClustersCount_point -> tileClustersCount;
//
//             tileHaveClustersCount_point [ shape = point ];
//             tileHaveClustersCount [ label="have clusters count flag" URL="\ref TileInfo::haveClustersCount_" shape=ellipse ];
//             tileHaveClustersCount_point -> tileHaveClustersCount;
//
//             tileSkippedTilesCount_point [ shape = point ];
//             tileSkippedTilesCount [ label="skipped tiles count" URL="\ref TileInfo::skippedTiles" shape=ellipse ];
//             tileSkippedTilesCount_point -> tileSkippedTilesCount;
//
//             tileSkippedClustersCount_point [ shape = point ];
//             tileSkippedClustersCount [ label="skipped clusters count" URL="\ref TileInfo::skippedClusters_" shape=ellipse ];
//             tileSkippedClustersCount_point -> tileSkippedClustersCount;
//
//         read_point [ shape = point ];
//         read [ label="read" URL="\ref ReadInfo" shape=box ];
//         read_point -> read [ label="1:n" URL="\ref Layout::readInfos_" ];
//
//             readNumber_point [ shape = point ];
//             readNumber [ label="number" URL="\ref ReadInfo::number_" shape=ellipse ];
//             readNumber_point -> readNumber;
//
//             indexReadFlag_point [ shape = point ];
//             indexReadFlag [ label="index read flag" URL="\ref ReadInfo::indexReadFlag_" shape=ellipse ];
//             indexReadFlag_point -> indexReadFlag;
//
//             maskAdapters_point [ shape = point ];
//             maskAdapters [ label="mask adapters" URL="\ref ReadInfo::maskAdapters_" shape=ellipse ];
//             maskAdapters_point -> maskAdapters;
//
//             trimAdapters_point [ shape = point ];
//             trimAdapters [ label="trim adapters" URL="\ref ReadInfo::trimAdapters_" shape=ellipse ];
//             trimAdapters_point -> trimAdapters;
//
//     {
//         rank=same;
//         layout -> flowcell_point -> lane_point -> read_point [ arrowhead=none ];
//     }
//     {
//         rank=same;
//         flowcell -> flowcellId_point -> runNumber_point -> instrument_point -> aggregateTilesFlag_point [ arrowhead=none ];
//     }
//     {
//         rank=same;
//         lane -> laneNumber_point -> sample_point -> tile_point [ arrowhead=none ];
//     }
//     {
//         rank=same;
//         sample -> sampleNumber_point -> sampleId_point -> sampleName_point -> barcode_point [ arrowhead=none ];
//     }
//     {
//         rank=same;
//         barcode -> component_point [ arrowhead=none ];
//     }
//     {
//         rank=same;
//         component -> componentLength_point -> componentValue_point [ arrowhead=none ];
//     }
//     {
//         rank=same;
//         tile -> tileNumber_point -> tileIndex_point -> tileClustersCount_point -> tileHaveClustersCount_point -> tileSkippedTilesCount_point -> tileSkippedClustersCount_point [ arrowhead=none ];
//     }
//     {
//         rank=same;
//         read -> readNumber_point -> indexReadFlag_point -> maskAdapters_point -> trimAdapters_point [ arrowhead=none ];
//     }


