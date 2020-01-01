/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file Bcl2FastqOptions.cpp
 *
 * \brief Implementation of the Bcl2Fastq options processing.
 *
 * \author Mauricio Varea
 * \author Marek Balint
 */


#include <string>
#include <sstream>
#include <iostream>

#include <cmath>

#include <boost/filesystem/operations.hpp>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

#include "config.h"
#include "config/Bcl2FastqOptions.hh"
#include <boost/algorithm/string.hpp>


namespace bcl2fastq {


namespace config {

const std::string createFastqForIndexReadsSetting = "create-fastq-for-index-reads";
const std::string findAdaptersWithSlidingWindowSetting = "find-adapters-with-sliding-window";
const std::string reverseComplementSetting = "write-fastq-reverse-complement";
const std::string useBasesMaskSetting = "use-bases-mask";

// This gets around a bug in boost::program_options::value<BoostPath>
// which doesn't allow spaces in the path name
std::istream& operator>>(std::istream& is, BoostPath& boostPath)
{
    std::string pathName;

    std::getline(is, pathName);
    boostPath.path_ = pathName;

    return is;
}

Bcl2FastqOptions::Bcl2FastqOptions()
: inputDir_()
, runfolderDir_()
, intensitiesDir_()
, outputDir_()
, interopDir_()
, statsDir_()
, reportsDir_()
, sampleSheetPath_()
, bclLoaderThreadsCount_(0)
, fastqCreatorThreadsCount_(0)
, fastqWriterThreadsCount_(0)
, tilesFilterString_()
, minimumTrimmedReadLength_()
, useBasesMasks_()
, maskShortAdapterReads_()
, adapterStringency_()
, generateReverseComplementFastqs_()
, includeNonPfClusters_()
, ignoreMissingBcls_()
, ignoreMissingFilters_()
, ignoreMissingPositions_()
, ignoreMissingControls_()
, createFastqsForIndexReads_()
, findAdaptersWithSlidingWindow_(false)
, noBgzfCompression_(false)
, fastqCompressionLevel_(1)
, barcodeMismatchCountsDelimitedString_()
, barcodeMismatchCounts_()
, autoSetToZeroBarcodeMismatches_(false)
, noLaneSplitting_(false)
, sampleSheet_()
{
}

common::ThreadVector::size_type Bcl2FastqOptions::getBclLoaderThreadsCount() const
{
    return bclLoaderThreadsCount_;
}

common::ThreadVector::size_type Bcl2FastqOptions::getFastqCreatorThreadsCount() const
{
    return fastqCreatorThreadsCount_;
}

common::ThreadVector::size_type Bcl2FastqOptions::getAppliedFastqWriterThreadsCount() const 
{
    return getFastqWriterThreadsCount() == 0 
        ? getBclLoaderThreadsCount()
        : getFastqWriterThreadsCount();
}

common::ThreadVector::size_type Bcl2FastqOptions::getFastqWriterThreadsCount() const
{
    return fastqWriterThreadsCount_;
}

std::vector<std::string> Bcl2FastqOptions::getTilesFilterList() const
{
    std::vector<std::string> tilesFilterList;

    if (!tilesFilterString_.empty())
    {
        boost::split(tilesFilterList, tilesFilterString_, boost::is_any_of(","));

        std::for_each(tilesFilterList.begin(),
                      tilesFilterList.end(),
                      boost::bind(boost::algorithm::trim<std::string>, _1, std::locale()));
    }

    return tilesFilterList;
}

std::size_t Bcl2FastqOptions::getMinimumTrimmedReadLength() const
{
    return minimumTrimmedReadLength_;
}

const Bcl2FastqOptions::UseBasesMaskContainer& Bcl2FastqOptions::getUseBasesMasks() const
{
    return useBasesMasks_;
}

std::size_t Bcl2FastqOptions::getMaskShortAdapterReads() const
{
    return maskShortAdapterReads_;
}

float Bcl2FastqOptions::getAdapterStringency() const
{
    return adapterStringency_;
}

bool Bcl2FastqOptions::getIgnoreMissingBcls() const
{
    return ignoreMissingBcls_;
}

bool Bcl2FastqOptions::getIgnoreMissingFilters() const
{
    return ignoreMissingFilters_;
}

bool Bcl2FastqOptions::getIgnoreMissingPositions() const
{
    return ignoreMissingPositions_;
}

bool Bcl2FastqOptions::getIgnoreMissingControls() const
{
    return ignoreMissingControls_;
}

bool Bcl2FastqOptions::getGenerateReverseComplementFastqs() const
{
    return generateReverseComplementFastqs_;
}

bool Bcl2FastqOptions::getIncludeNonPfClusters() const
{
    return includeNonPfClusters_;
}

bool Bcl2FastqOptions::getCreateFastqsForIndexReads() const
{
    return createFastqsForIndexReads_;
}

bool Bcl2FastqOptions::useBgzfCompression() const
{
    return !noBgzfCompression_;
}

int Bcl2FastqOptions::getFastqCompressionLevel() const
{
    return fastqCompressionLevel_;
}

void Bcl2FastqOptions::initNamedOptions(boost::program_options::options_description &options)
{
    options.add_options()(
        "input-dir,i",
        boost::program_options::value<BoostPath>(&inputDir_)
             ->default_value(BoostPath(),"<runfolder-dir>/Data/Intensities/BaseCalls/"),
        "path to input directory"
    );

    options.add_options()(
        "runfolder-dir,R",
        boost::program_options::value<BoostPath>(&runfolderDir_)
             ->default_value(BoostPath(boost::filesystem::current_path()),"./"),
        "path to runfolder directory"
    );

    options.add_options()(
        "intensities-dir",
        boost::program_options::value<BoostPath>(&intensitiesDir_)
             ->default_value(BoostPath(),"<input-dir>/../"),
        "path to intensities directory\n"
        "If intensities directory is specified, --input-dir must also be specified."
    );

    options.add_options()(
        "output-dir,o",
        boost::program_options::value<BoostPath>(&outputDir_)
             ->default_value(BoostPath(),"<input-dir>"),
        "path to demultiplexed output"
    );

    options.add_options()(
        "interop-dir",
        boost::program_options::value<BoostPath>(&interopDir_)
             ->default_value(BoostPath(),"<runfolder-dir>/InterOp/"),
        "path to demultiplexing statistics directory"
    );

    options.add_options()(
        "stats-dir",
        boost::program_options::value<BoostPath>(&statsDir_)
             ->default_value(BoostPath(),"<output-dir>/Stats/"),
        "path to human-readable demultiplexing statistics directory"
    );

    options.add_options()(
        "reports-dir",
        boost::program_options::value<BoostPath>(&reportsDir_)
             ->default_value(BoostPath(),"<output-dir>/Reports/"),
        "path to reporting directory"
    );

    options.add_options()(
        "sample-sheet",
        boost::program_options::value<BoostPath>(&sampleSheetPath_)
             ->default_value(BoostPath(),"<runfolder-dir>/SampleSheet.csv"),
        "path to the sample sheet"
    );

    options.add_options()(
        "loading-threads,r",
        boost::program_options::value<common::ThreadVector::size_type>(&bclLoaderThreadsCount_)
             ->default_value(4),
        "number of threads used for loading BCL data"
    );

    options.add_options()(
        "processing-threads,p",
        boost::program_options::value<common::ThreadVector::size_type>(&fastqCreatorThreadsCount_),
        "number of threads used for processing demultiplexed data"
    );

    options.add_options()(
        "writing-threads,w",
        boost::program_options::value<common::ThreadVector::size_type>(&fastqWriterThreadsCount_)
             ->default_value(4),
        "number of threads used for writing FASTQ data.\n"
        "This should not be set higher than the number of samples. "
        "If set =0 then these threads will be placed in the "
        "same pool as the loading threads, and the number "
        "of shared threads will be determined by --loading-threads."
    );

    options.add_options()(
        "tiles",
        boost::program_options::value<std::string>(&tilesFilterString_),
        "comma-separated list of regular expressions to select only a subset of the tiles available in the flow-cell."
        " Multiple entries allowed, each applies to the corresponding base-calls.\n"
        "For example:\n"
        " * to select all the tiles ending with '5' in all lanes:\n     --tiles [0-9][0-9][0-9]5\n"
        " * to select tile 2 in lane 1 and all the tiles in the other lanes:\n     --tiles s_1_0002,s_[2-8]"
    );

    options.add_options()(
        "minimum-trimmed-read-length",
        boost::program_options::value<std::size_t>(&minimumTrimmedReadLength_)
             ->default_value(35),
        "minimum read length after adapter trimming"
    );

    options.add_options()(
        useBasesMaskSetting.c_str(),
        boost::program_options::value<UseBasesMaskContainer>(&useBasesMasks_)
             ->default_value(UseBasesMaskContainer(), ""),
        "specifies how to use each cycle."
    );

    options.add_options()(
        "mask-short-adapter-reads",
        boost::program_options::value<std::size_t>(&maskShortAdapterReads_)
             ->default_value(22),
        "smallest number of remaining bases (after masking bases below the minimum trimmed read length) below which whole read is masked"
    );

    options.add_options()(
        "adapter-stringency",
        boost::program_options::value<float>(&adapterStringency_)
             ->default_value(0.9,"0.9"),
        "adapter stringency"
    );

    options.add_options()(
        "ignore-missing-bcls",
        boost::program_options::bool_switch(&ignoreMissingBcls_)
             ->default_value(false),
        "assume 'N'/'#' for missing calls"
    );

    options.add_options()(
        "ignore-missing-filter",
        boost::program_options::bool_switch(&ignoreMissingFilters_)
             ->default_value(false),
        "assume 'true' for missing filters"
    );

    options.add_options()(
        "ignore-missing-positions",
        boost::program_options::bool_switch(&ignoreMissingPositions_)
             ->default_value(false),
        "assume [0,i] for missing positions, where i is incremented starting from 0"
    );

    options.add_options()(
        "ignore-missing-controls",
        boost::program_options::bool_switch(&ignoreMissingControls_)
             ->default_value(false),
        "(deprecated) assume 0 for missing controls"
    );

    options.add_options()(
        reverseComplementSetting.c_str(),
        boost::program_options::bool_switch(&generateReverseComplementFastqs_)
             ->default_value(false),
        "generate FASTQs containing reverse complements of actual data"
    );

    options.add_options()(
        "with-failed-reads",
        boost::program_options::bool_switch(&includeNonPfClusters_)
             ->default_value(false),
        "include non-PF clusters"
    );

    options.add_options()(
        createFastqForIndexReadsSetting.c_str(),
        boost::program_options::bool_switch(&createFastqsForIndexReads_)
             ->default_value(false),
        "create FASTQ files also for index reads"
    );

    options.add_options()(
        findAdaptersWithSlidingWindowSetting.c_str(),
        boost::program_options::bool_switch(&findAdaptersWithSlidingWindow_)
             ->default_value(false),
        "find adapters with simple sliding window algorithm"
    );

    options.add_options()(
        "no-bgzf-compression",
        boost::program_options::bool_switch(&noBgzfCompression_)
             ->default_value(false),
        "turn off BGZF compression for FASTQ files"
    );

    options.add_options()(
        "fastq-compression-level",
        boost::program_options::value<int>(&fastqCompressionLevel_)
             ->default_value(4),
        "zlib compression level (1-9) used for FASTQ files"
    );

    options.add_options()(
        "barcode-mismatches",
        boost::program_options::value<std::string>(&barcodeMismatchCountsDelimitedString_)
             ->default_value("1"),
        "number of allowed mismatches per index\n"
        "Multiple, comma delimited, entries allowed. Each entry is applied to the corresponding index; "
        "last entry applies to all remaining indices.\n"
        "Accepted values: 0, 1, 2."
    );

    options.add_options()(
        "no-lane-splitting",
        boost::program_options::bool_switch(&noLaneSplitting_)
             ->default_value(false),
        "do not split fastq files by lane."
    );
}

void Bcl2FastqOptions::initUnnamedOptions(boost::program_options::options_description &options)
{
    options.add_options()(
        "auto-set-to-zero-barcode-mismatches",
        boost::program_options::bool_switch(&autoSetToZeroBarcodeMismatches_)
             ->default_value(false),
        "automatically set allowed barcode mismatches to 0 on collision"
    );
}

void Bcl2FastqOptions::initPositionalOptions(boost::program_options::positional_options_description &options)
{
}

std::string Bcl2FastqOptions::usagePrefix() const
{
    std::ostringstream os;

    os << "Usage:" << std::endl;
    os << "      " << BCL2FASTQ_NAME_SHORT << " [options]" << std::endl;

    return os.str();
}

std::string Bcl2FastqOptions::usageSuffix() const
{
    return "";
}

void Bcl2FastqOptions::storeSampleSheetSettings(boost::program_options::variables_map &vm)
{
    /// \todo Refactoring: Is there a better way for obtaining platform-independent path separator? (Consider CMake)
    auto pathSeparator = boost::filesystem::path("/").make_preferred().native();

    if (inputDir_.path_.empty())
    {
        inputDir_.path_ =
            runfolderDir_.path_
            / boost::filesystem::path("Data")
            / boost::filesystem::path("Intensities")
            / boost::filesystem::path("BaseCalls")
            / pathSeparator
        ;
    }

    // output-dir
    if (outputDir_.path_.empty())
    {
        outputDir_.path_ = inputDir_.path_;
    }

    // sample sheet location
    if (sampleSheetPath_.path_.empty())
    {
        sampleSheetPath_.path_ =
            runfolderDir_.path_
            / boost::filesystem::path("SampleSheet.csv");
    }
    else if (!boost::filesystem::exists(sampleSheetPath_.path_))
    {
        BCL2FASTQ_LOG(common::LogLevel::WARNING) << "Sample sheet " << sampleSheetPath_.path_ << " does not exist" << std::endl;
    }

    if (!(sampleSheet_ = config::createSampleSheetCsv(sampleSheetPath_.path_,
                                                      outputDir_.path_)))
    {
        return;
    }

    // For settings not passed on the command line, look for them in the sample sheet.
    if (vm[createFastqForIndexReadsSetting].defaulted())
    {
        SampleSheetCsv::TriState createFastqForIndexReads = sampleSheet_->createFastqForIndexReads();
        if (createFastqForIndexReads != SampleSheetCsv::INDETERMINATE)
        {
            createFastqsForIndexReads_ = (createFastqForIndexReads == SampleSheetCsv::TRUE);
        }

    }
    if (vm[findAdaptersWithSlidingWindowSetting].defaulted())
    {
        SampleSheetCsv::TriState findAdaptersWithIndels = sampleSheet_->findAdaptersWithIndels();
        if (findAdaptersWithIndels != SampleSheetCsv::INDETERMINATE)
        {
            findAdaptersWithSlidingWindow_ = (findAdaptersWithIndels != SampleSheetCsv::TRUE);
        }
    }
    if (vm[reverseComplementSetting].defaulted())
    {
        SampleSheetCsv::TriState reverseComplement = sampleSheet_->generateReverseComplementFastqs();
        if (reverseComplement != SampleSheetCsv::INDETERMINATE)
        {
            generateReverseComplementFastqs_ = (reverseComplement == SampleSheetCsv::TRUE);
        }
    }
}

void Bcl2FastqOptions::postProcess(boost::program_options::variables_map &vm)
{
    /// \todo Refactoring: Is there a better way for obtaining platform-independent path separator? (Consider CMake)
    auto pathSeparator = boost::filesystem::path("/").make_preferred().native();

    BCL2FASTQ_LOG(common::LogLevel::INFO) << "Runfolder path: '" << runfolderDir_.path_.string() << "'" << std::endl;
    BCL2FASTQ_LOG(common::LogLevel::INFO) << "Input path: '" << inputDir_.path_.string() << "'" << std::endl;

    // intensities-dir
    if (intensitiesDir_.path_.empty())
    {
        intensitiesDir_.path_ = boost::filesystem::path(inputDir_.path_ / pathSeparator).parent_path().parent_path() / pathSeparator;
    }
    BCL2FASTQ_LOG(common::LogLevel::INFO) << "Intensities path: '" << intensitiesDir_.path_.string() << "'" << std::endl;

    BCL2FASTQ_LOG(common::LogLevel::INFO) << "Output path: '" << outputDir_.path_.string() << "'" << std::endl;

    // interop-dir
    if (interopDir_.path_.empty())
    {
        interopDir_.path_ = 
            runfolderDir_.path_
            / boost::filesystem::path("InterOp")
            / pathSeparator
        ;
    }
    BCL2FASTQ_LOG(common::LogLevel::INFO) << "InterOp path: '" << interopDir_.path_.string() << "'" << std::endl;

    // stats-dir
    if (statsDir_.path_.empty())
    {
        statsDir_.path_ =
            outputDir_.path_
            / boost::filesystem::path("Stats")
            / pathSeparator
        ;
    }
    BCL2FASTQ_LOG(common::LogLevel::INFO) << "Stats path: '" << statsDir_.path_.string() << "'" << std::endl;

    // reports-dir
    if (reportsDir_.path_.empty())
    {
        reportsDir_.path_ =
            outputDir_.path_
            / boost::filesystem::path("Reports")
            / pathSeparator
        ;
    }
    BCL2FASTQ_LOG(common::LogLevel::INFO) << "Reports path: '" << reportsDir_.path_.string() << "'" << std::endl;

    // thread counts
    {
        std::size_t cpusCount = boost::thread::hardware_concurrency();
        BCL2FASTQ_LOG(common::LogLevel::INFO) << "Detected CPUs: " << cpusCount << std::endl;

        BCL2FASTQ_LOG(common::LogLevel::INFO) << "Loading threads: " << bclLoaderThreadsCount_ << std::endl;
        if (fastqCreatorThreadsCount_ == 0)
        {
            fastqCreatorThreadsCount_ = cpusCount;
        }
        BCL2FASTQ_LOG(common::LogLevel::INFO) << "Processing threads: " << fastqCreatorThreadsCount_ << std::endl;
        BCL2FASTQ_LOG(common::LogLevel::INFO) << "Writing threads: " << getAppliedFastqWriterThreadsCount() << std::endl;
    }

    // barcode mismatches
    boost::tokenizer<boost::char_separator<char> > tknzr(
        barcodeMismatchCountsDelimitedString_,
        boost::char_separator<char>(",:|+")
    );
    std::transform(
        tknzr.begin(),
        tknzr.end(),
        std::back_inserter(barcodeMismatchCounts_),
        static_cast<std::size_t (*) (const std::string&)>(&boost::lexical_cast<std::size_t>)
    );

    std::stringstream barcodeMismatchCountsDelimitedString;
    std::copy(barcodeMismatchCounts_.begin(), barcodeMismatchCounts_.end(), std::ostream_iterator<MismatchesCount>(barcodeMismatchCountsDelimitedString, " "));
    BCL2FASTQ_LOG(common::LogLevel::INFO) << "Allowed barcode mismatches: " << barcodeMismatchCountsDelimitedString.str() << std::endl;

    BCL2FASTQ_LOG(common::LogLevel::INFO) << "Tiles: " << (tilesFilterString_.empty() ? "<ALL>" : tilesFilterString_) << std::endl;;
    BCL2FASTQ_LOG(common::LogLevel::INFO) << "Minimum trimmed read length: " << minimumTrimmedReadLength_ << std::endl;

    std::ostringstream useBasesMaskInfoStr;
    BOOST_FOREACH(const std::string useBasesMask, useBasesMasks_)
    {
        useBasesMaskInfoStr << std::endl << "\t" << useBasesMask;
    }
    BCL2FASTQ_LOG(common::LogLevel::INFO) << "Use bases masks: " << (useBasesMasks_.empty() ? "<NONE>" : useBasesMaskInfoStr.str()) << std::endl;

    BCL2FASTQ_LOG(common::LogLevel::INFO) << "Mask short adapter reads: " << maskShortAdapterReads_ << std::endl;
    BCL2FASTQ_LOG(common::LogLevel::INFO) << "Adapter stringency: " << adapterStringency_ << std::endl;
    const std::string trimMethodText = (findAdaptersWithSlidingWindow_) ?
        "Sliding window" : "Allow matches with indels";
    BCL2FASTQ_LOG(common::LogLevel::INFO) << "Adapter trimming method: " << trimMethodText << std::endl;
    BCL2FASTQ_LOG(common::LogLevel::INFO) << "Ignore missing BCLs: " << (ignoreMissingBcls_ ? "YES" : "NO") << std::endl;
    BCL2FASTQ_LOG(common::LogLevel::INFO) << "Ignore missing filters: " << (ignoreMissingFilters_ ? "YES" : "NO") << std::endl;
    BCL2FASTQ_LOG(common::LogLevel::INFO) << "Ignore missing positions: " << (ignoreMissingPositions_ ? "YES" : "NO") << std::endl;
    BCL2FASTQ_LOG(common::LogLevel::INFO) << "Ignore missing controls: " << (ignoreMissingControls_ ? "YES" : "NO") << std::endl;
    BCL2FASTQ_LOG(common::LogLevel::INFO) << "Include non-PF clusters: " << (includeNonPfClusters_ ? "YES" : "NO") << std::endl;
    BCL2FASTQ_LOG(common::LogLevel::INFO) << "Create FASTQs for index reads: " << (createFastqsForIndexReads_ ? "YES" : "NO") << std::endl;
    BCL2FASTQ_LOG(common::LogLevel::INFO) << "Use bgzf compression for FASTQ files: " << (noBgzfCompression_ ? "NO" : "YES") << std::endl;
    BCL2FASTQ_LOG(common::LogLevel::INFO) << "FASTQ compression level: " << fastqCompressionLevel_ << std::endl;
}


} // namespace config
} // namespace bcl2fastq


