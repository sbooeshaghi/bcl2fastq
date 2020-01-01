/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file SampleSheetCsv.hh
 *
 * \brief Declaration of sample sheet helper.
 *
 * \author Marek Balint
 */


#ifndef BCL2FASTQ_CONFIG_SAMPLESHEETCSV_HH
#define BCL2FASTQ_CONFIG_SAMPLESHEETCSV_HH


#include <vector>
#include <string>

#include <boost/filesystem/path.hpp>
#include <boost/noncopyable.hpp>
#include "common/CsvGrammar.hh"
#include "common/Types.hh"
#include "common/SampleMetadata.hh"

#include <boost/shared_ptr.hpp>

namespace bcl2fastq {


namespace config {

/// \brief Sample sheet helper.
class SampleSheetCsv : private boost::noncopyable
{
public:

    class ExcludedTileRange
    {
    public:
        ExcludedTileRange(const std::string& flowcell,
                          common::LaneNumber laneNumber,
                          common::TileNumber firstTile,
                          common::TileNumber lastTile);

        bool isExcluded(const std::string& flowcell,
                        common::LaneNumber laneNumber,
                        common::TileNumber tileNumber) const;

    private:
        std::string        flowcell_;
        common::LaneNumber laneNumber_;
        common::TileNumber firstTile_;
        common::TileNumber lastTile_;
    };

    /// \brief Adapters container type definition.
    typedef std::vector<std::pair<common::ReadNumber,std::string> > AdaptersContainer;

    /// \brief Excluded tiles container type definition
    typedef std::vector<ExcludedTileRange> ExcludedTilesContainer;

public:

    /// \brief Constructor.
    /// \param sampleSheetData Sample sheet CSV data.
    /// \param directoryValidator Use to validate that there are no clashes
    SampleSheetCsv(const common::CsvGrammarAttribute &sampleSheetData,
                   const boost::filesystem::path&    outputDir);

public:

    // I would like to use boost::logic::tribool, but this expression bizarrely returns false:
    // boost::logic::indeterminate(boost::tribool::indeterminate_value)
    enum TriState { TRUE=0, FALSE, INDETERMINATE };

    /// \brief Get beginning of adapters to be masked.
    /// \return Iterator to beginning of adapters to be masked.
    SampleSheetCsv::AdaptersContainer::const_iterator maskAdaptersBegin() const;

    /// \brief Get end of adapters to be masked.
    /// \return Iterator to end of adapters to be trimmed.
    SampleSheetCsv::AdaptersContainer::const_iterator maskAdaptersEnd() const;

    /// \brief Get beginning of adapters to be trimmed.
    /// \return Iterator to beginning of adapters to be trimmed.
    SampleSheetCsv::AdaptersContainer::const_iterator trimAdaptersBegin() const;

    /// \brief Get end of adapters to be trimmed.
    /// \return Iterator to end of adapters to be trimmed.
    SampleSheetCsv::AdaptersContainer::const_iterator trimAdaptersEnd() const;

    /// \brief Get sample metadata container.
    /// \return Sample metadata container.
    const common::SampleMetadataContainer& getSampleMetadata() const { return samples_; }

    /// \brief Get the setting for createFastqForIndexReads
    /// \return INDETERMINATE if the setting is not in the sample sheet.
    TriState createFastqForIndexReads() const { return createFastqForIndexReads_; }

    /// \brief Get the setting for findAdaptersWithIndels
    /// \return INDETERMINATE if the setting is not in the sample sheet.
    TriState findAdaptersWithIndels() const { return findAdaptersWithIndels_; }

    /// \brief Generate reverse complement fastqs
    /// \return INDTERMINATE if the setting is not in the sample sheet.
    TriState generateReverseComplementFastqs() const { return generateReverseComplementFastqs_; }

    /// \brief Get the read 1 end cycle
    /// \return Read 1 end cycle number.
    common::CycleNumber getRead1EndCycle() const { return read1EndWithCycle_; }

    /// \brief Get the read 2 end cycle
    /// \return Read 2 end cycle number.
    common::CycleNumber getRead2EndCycle() const { return read2EndWithCycle_; }

    /// \brief Get the read 1 start cycle
    /// \return Read 1 start cycle number.
    common::CycleNumber getRead1StartCycle() const { return read1StartFromCycle_; }

    /// \brief Get the read 1 UMI length
    /// \return Read 1 UMI length.
    common::CycleNumber getRead1UmiLength() const { return read1UmiLength_; }

    /// \brief Get the read 2 UMI length
    /// \return Read 2 UMI length.
    common::CycleNumber getRead2UmiLength() const { return read2UmiLength_; }

    /// \brief Get the read 1 UMI start cycle
    /// \return Read 1 UMI start cycle
    common::CycleNumber getRead1UmiStartFromCycle() const { return read1UmiStartFromCycle_; }

    /// \brief Get the read 2 UMI start cycle
    /// \return Read 2 UMI start cycle
    common::CycleNumber getRead2UmiStartFromCycle() const { return read2UmiStartFromCycle_; }

    /// \brief Trim UMI if true
    /// \return INDTERMINATE if the setting is not in the sample sheet.
    TriState trimUmi() const { return trimUmi_; }

    /// \brief Get the read 2 start cycle
    /// \return Read 2 start cycle number.
    common::CycleNumber getRead2StartCycle() const { return read2StartFromCycle_; }

    /// \brief Return true if the tile number is excluded.
    /// \param flowcell Flowcell id.
    /// \param laneNumber Lane number.
    /// \param tileNumber Tile number.
    /// \return True if tile is excluded.
    bool isTileExcluded(const std::string& flowcell,
                        common::LaneNumber laneNumber,
                        common::TileNumber tileNumber) const;

private:

    /// \brief Adapters to be masked.
    AdaptersContainer maskAdapters_;

    /// \brief Adapters to be trimmed.
    AdaptersContainer trimAdapters_;

    /// \brief Sample metadata.
    common::SampleMetadataContainer samples_;

    /// \brief Create fastq for index reads.
    TriState createFastqForIndexReads_;

    /// \brief Find adapters with indels.
    TriState findAdaptersWithIndels_;

    /// \brief Generate reverse complement fastqs
    TriState generateReverseComplementFastqs_;

    /// \brief End cycle number for read 1
    common::CycleNumber read1EndWithCycle_;

    /// \brief End cycle number for read 2
    common::CycleNumber read2EndWithCycle_;

    /// \brief Start cycle number for read 1
    common::CycleNumber read1StartFromCycle_;

    /// \brief Start cycle number for read 2
    common::CycleNumber read2StartFromCycle_; 

    /// \brief UMI length for read 1
    common::CycleNumber read1UmiLength_;

    /// \brief UMI length for read 2
    common::CycleNumber read2UmiLength_;

    /// \brief Read 1 UMI start cycle
    common::CycleNumber read1UmiStartFromCycle_;

    /// \brief Read 2 UMI start cycle
    common::CycleNumber read2UmiStartFromCycle_;

    /// \brief Trim UMI if true
    TriState trimUmi_;

    /// \brief Excluded tiles
    ExcludedTilesContainer excludedTiles_;
};

typedef std::shared_ptr<SampleSheetCsv> SampleSheetCsvPtr;

/// \brief Create sample sheet object.
/// \param runfolderDir Path to runfolder.
/// \param directoryValidator Used to validate there are no clashes with existing directories
/// \return Sample sheet object constructed with data from given runfolder.
SampleSheetCsvPtr createSampleSheetCsv(const boost::filesystem::path& sampleSheetPath,
                                       const boost::filesystem::path& outputDir);


} // namespace config


} // namespace bcl2fastq


#endif //BCL2FASTQ_CONFIG_SAMPLESHEETCSV_HH


