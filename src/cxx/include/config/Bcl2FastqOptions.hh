/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file Bcl2FastqOptions.hh
 *
 * \brief Declaration of the Bcl2Fastq options processing.
 *
 * \author Mauricio Varea
 * \author Marek Balint
 */

#ifndef BCL2FASTQ_CONFIG_BCL2FASTQOPTIONS_HH
#define BCL2FASTQ_CONFIG_BCL2FASTQOPTIONS_HH


#include <iosfwd>
#include <string>

#include <boost/filesystem/path.hpp>
#include <boost/optional.hpp>

#include "config/SampleSheetCsv.hh"

#include "common/Options.hh"
#include "common/Threads.hh"


namespace bcl2fastq {
namespace config {

/// \brief Maximum number of allowed mismatches in a barcode component type definition.
typedef std::size_t MismatchesCount;

/// \brief Maximum allowed mismatches mismatches in a barcode component container type definition.
typedef std::vector<MismatchesCount> BarcodeMismatchCountsContainer;

/// \brief This class circumvents a bug in boost that produces an error when
/// spaces occur in boost::filesystem::path used with boost::program_options::value
class BoostPath
{
public:
    BoostPath() : path_() { }
    BoostPath(const boost::filesystem::path& path) : path_(path) { }

    friend std::istream& operator>>(std::istream& is, BoostPath& boostPath);

    boost::filesystem::path path_;
};


/// \brief Bcl2Fastq custom options. 
class Bcl2FastqOptions : public common::Options
{
public:

    /// \brief UseBasesMasks container type definition.
    typedef std::vector<std::string> UseBasesMaskContainer;

public:

    /// \brief Default constructor.
    Bcl2FastqOptions();

public:

    /// \brief Get input directory.
    /// \return Path to the directory containing input BCL files.
    const boost::filesystem::path& getInputDir() const { return inputDir_.path_; }

    /// \brief Get runfolder directory.
    /// \return Path to the project's root directory.
    const boost::filesystem::path& getRunfolderDir() const { return runfolderDir_.path_; }

    /// \brief Get intensities directory.
    /// \return Path to the intensities directory.
    const boost::filesystem::path& getIntensitiesDir() const { return intensitiesDir_.path_; }

    /// \brief Get output directory.
    /// \return Path to demultiplexed output directory.
    const boost::filesystem::path& getOutputDir() const { return outputDir_.path_; }

    /// \brief Get interop directory.
    /// \return Path to demultiplexing statistics directory.
    const boost::filesystem::path& getInteropDir() const { return interopDir_.path_; }

    /// \brief Get stats directory.
    /// \return Path to human-readable demultiplexing statistics directory.
    const boost::filesystem::path& getStatsDir() const { return statsDir_.path_; }

    /// \brief Get reports directory.
    /// \return Path to reporting directory.
    const boost::filesystem::path& getReportsDir() const { return reportsDir_.path_; }

    /// \brief Get sample sheet path .
    /// \return Path to human-readable demultiplexing statistics directory.
    const boost::filesystem::path& getSampleSheetPath() const { return sampleSheetPath_.path_; }

    /// \brief Get number of BCL loading threads.
    /// \return Number of threds to be used for loading BCL data.
    common::ThreadVector::size_type getBclLoaderThreadsCount() const;

    /// \brief Get number of FASTQ creating threads.
    /// \return Number of threads to be used for creation of FASTQ data.
    common::ThreadVector::size_type getFastqCreatorThreadsCount() const;

    /// \brief Get number of FASTQ writing threads.
    /// \return Number of threads to be used for writing FASTQ data.
    common::ThreadVector::size_type getFastqWriterThreadsCount() const;

    /// \brief Get list of tiles to be processed.
    /// \return List of tiles to be processed.
    std::vector<std::string> getTilesFilterList() const;

    /// \brief Get minimum read length after adapter trimming.
    /// \return Minimum read length after adapter trimming.
    std::size_t getMinimumTrimmedReadLength() const;

    /// \brief Get the masks used to select which cycled are used.
    /// \return The masks
    const UseBasesMaskContainer& getUseBasesMasks() const;

    /// \brief Get maximum number of useful bases in read after adapter trimming for which whole read is masked.
    /// \return Maximum number of useful bases in read after adapter trimming for which whole read is masked.
    std::size_t getMaskShortAdapterReads() const;

    /// \brief Get adapter stringency.
    /// \return Adapter stringency.
    float getAdapterStringency() const;

    /// \brief Get flag for ignore missing BCLs.
    /// \retval true Ignore missing BCLs.
    /// \retval false Fail on missing BCLs.
    bool getIgnoreMissingBcls() const;

    /// \brief Get flag for ignore missing filters.
    /// \retval true Ignore missing filters.
    /// \retval false Fail on missing filters.
    bool getIgnoreMissingFilters() const;

    /// \brief Get flag for ignore missing positions.
    /// \retval true Ignore missing positions.
    /// \retval false Fail on missing positions.
    bool getIgnoreMissingPositions() const;

    /// \brief Get flag for ignore missing controls.
    /// \retval true Ignore missing controls.
    /// \retval false Fail on missing controls.
    bool getIgnoreMissingControls() const;

    /// \brief Get flag for generate reverse complement FASTQs.
    /// \retval true Generate reverse complement FASTQs.
    /// \retval false Do not generate reverse complement FASTQs.
    bool getGenerateReverseComplementFastqs() const;

    /// \brief Get flag for include non-PF clusters.
    /// \retval true Include non-PF clusters.
    /// \retval false Do not include non-PF clusters.
    bool getIncludeNonPfClusters() const;

    /// \brief Get flag for create FASTQ files also for index reads.
    /// \retval true Create FASTQ files also for index reads.
    /// \retval false Do not create FASTQ files for index reads.
    bool getCreateFastqsForIndexReads() const;

    /// \brief Find adapters with indels.
    /// \retval true if should find adapters with sliding window.
    bool findAdaptersWithSlidingWindow() const { return findAdaptersWithSlidingWindow_; }

    /// \brief Get flag for bgzf compression use.
    /// retval True if bgzf compression should be used for fastq files
    bool useBgzfCompression() const;

    /// \brief Get compression level used by zlib for fastq files
    /// retval compression level to be used by zlib
    int getFastqCompressionLevel() const;

    /// \brief Get flag for determining if it should automatically set allowed barcode mismatches to 0 on collision.
    /// retval True if it should automatically set allowed barcode mismatches to 0 on collision.
    bool getAutoSetToZeroBarcodeMismatches() const { return autoSetToZeroBarcodeMismatches_; }

    /// \brief Get allowed barcode mismatch counts.
    /// \return Barcode mismatch counts.
    BarcodeMismatchCountsContainer& getBarcodeMismatches() { return barcodeMismatchCounts_; }

    /// \brief Get allowed barcode mismatch counts.
    /// \return Barcode mismatch counts.
    const BarcodeMismatchCountsContainer& getBarcodeMismatches() const { return barcodeMismatchCounts_; }

    /// \brief Return true if fastq files should not be split by lane
    /// \return True if fastq files should not be split by lane.
    bool noLaneSplitting() const { return noLaneSplitting_; }

    const config::SampleSheetCsv& getSampleSheet() const { return *sampleSheet_; }

public:

    common::ThreadVector::size_type getAppliedFastqWriterThreadsCount() const;

private:

    virtual void initNamedOptions(boost::program_options::options_description &options);

    virtual void initUnnamedOptions(boost::program_options::options_description &options);

    virtual void initPositionalOptions(boost::program_options::positional_options_description &options);

private:

    virtual std::string usagePrefix() const;

    virtual std::string usageSuffix() const;

private:

    /// \brief Store the options from the sample sheet.
    /// \param vm Options map where the settings are stored.
    virtual void storeSampleSheetSettings(boost::program_options::variables_map &vm);

    /// \brief Post-process parsed command line arguments.
    /// \param vm Command line arguments to be post-process.
    virtual void postProcess(boost::program_options::variables_map &vm);

private:

    /// \brief Path to BaseCalls directory.
    BoostPath inputDir_;

    /// \brief Path to BaseCalls directory.
    BoostPath runfolderDir_;

    /// \brief Path to intensities directory.
    BoostPath intensitiesDir_;

    /// \brief Path to demultiplexed output.
    BoostPath outputDir_;

    /// \brief Path to demultiplexing stats directory.
    BoostPath interopDir_;

    /// \brief Path to human-readable demultiplexing stats directory.
    BoostPath statsDir_;

    /// \brief Path to reporting directory.
    BoostPath reportsDir_;

    /// \brief Path to sample sheet.
    BoostPath sampleSheetPath_;

    /// \brief Number of threads to be used for loading BCL data.
    common::ThreadVector::size_type bclLoaderThreadsCount_;

    /// \brief Number of threads to be used for demultiplexing.
    common::ThreadVector::size_type demultiplexerThreadsCount_;

    /// \brief Number of threads to be used for creation of FASTQ data.
    common::ThreadVector::size_type fastqCreatorThreadsCount_;

    /// \brief Number of threads to be used for writing FASTQ data.
    common::ThreadVector::size_type fastqWriterThreadsCount_;

    /// \brief List of tiles to be processed.
    std::string tilesFilterString_;

    /// \brief Minimum read length after adapter trimming.
    std::size_t minimumTrimmedReadLength_;

    /// \brief Identifies which cycles should be used
    UseBasesMaskContainer useBasesMasks_;

    /// \brief Maximum number of useful bases in read after adapter trimming for which whole read is masked.
    std::size_t maskShortAdapterReads_;

    /// \brief Adapter stringency.
    float adapterStringency_;

    /// \brief Generate reverse complement FASTQs flag.
    bool generateReverseComplementFastqs_;

    /// \brief Include non-PF clusters in created FASTQ files.
    bool includeNonPfClusters_;

    /// \brief Ignore missing BCLs (assume 'N'/'#' for missing calls).
    bool ignoreMissingBcls_;

    /// \brief Ignore missing filters (assume 'true' for missing filters).
    bool ignoreMissingFilters_;

    /// \brief Ignore missing positions (assume [0,0] for missing positions).
    bool ignoreMissingPositions_;

    /// \brief Ignore missing controls (assume 0 for missing controls).
    bool ignoreMissingControls_;

    /// \brief Create FASTQ files also for index reads.
    bool createFastqsForIndexReads_;

    /// \brief Find adapters with sliding window.
    bool findAdaptersWithSlidingWindow_;

    /// \brief Turn off bgzf compression.
    bool noBgzfCompression_;

    /// \brief Compression level used by zlib for fastq files
    int fastqCompressionLevel_;

    /// \brief Barcode mismatch counts delimited string.
    std::string barcodeMismatchCountsDelimitedString_;

    /// \brief Barcode mismatch counts.
    BarcodeMismatchCountsContainer barcodeMismatchCounts_;

    /// \brief If true, automatically set allowed barcode mismatches to 0 on collision
    bool autoSetToZeroBarcodeMismatches_;

    /// \brief If true, don't split fastq files by lane
    bool noLaneSplitting_;

    /// \brief Sample sheet.
    config::SampleSheetCsvPtr sampleSheet_;
};


} // namespace config
} // namespace bcl2fastq


#endif // BCL2FASTQ_CONFIG_BCL2FASTQOPTIONS_HH


