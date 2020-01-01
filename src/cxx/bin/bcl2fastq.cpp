/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file bcl2fastq.cpp
 *
 * \brief Bcl2Fastq main.
 *
 * \author Marek Balint
 * \author Mauricio Varea
 */

#include <boost/cstdlib.hpp>

#include <algorithm>

#include <boost/bind.hpp>
#include <boost/cstdlib.hpp>

#include "common/Program.hh"
#include "common/ProgramInfo.hh"
#include "config/Bcl2FastqOptions.hh"

#include "layout/RunInfoXml.hh"
#include "layout/Layout.hh"
#include "layout/LaneInfo.hh"
#include "conversion/Converter.hh"
#include "stats/DemuxReportGenerator.hh"


namespace detail {

static const unsigned TOP_UNKNOWN_BARCODES_MAX = 1000;

static const std::string interopStatsFileName("IndexMetricsOut.bin");
static const std::string demuxStatsFileName("DemultiplexingStats.xml");
static const std::string conversionStatsFileName("ConversionStats.xml");
static const std::string statsJsonFileName("Stats.json");
static const std::string adapterStatsFileName("AdapterTrimming.txt");
static const std::string fastqSummaryPartialFileName("FastqSummaryF1L");
static const std::string demuxSummaryPartialFileName("DemuxSummaryF1L");

/// \brief Process single lane.
/// \param layout Flowcell layout.
/// \param options Program options.
/// \param laneInfo Lane metadata.
static void processLane(
    const bcl2fastq::config::Bcl2FastqOptions &options,
    const bcl2fastq::layout::Layout &layout,
    const bcl2fastq::layout::LaneInfo &laneInfo,
    bcl2fastq::conversion::ConverterStats &stats
)
{
    bcl2fastq::conversion::Converter( options, layout, laneInfo, stats ).run();
    stats.getUnknownBarcodeStat( laneInfo ).findPopular( TOP_UNKNOWN_BARCODES_MAX, laneInfo.getNumber() );
}


} // namespace detail


static void program(
    const bcl2fastq::common::ProgramInfo &info,
    bcl2fastq::config::Bcl2FastqOptions &options
)
{
    bcl2fastq::common::adjustMaxFileHandles();

    const bcl2fastq::layout::RunInfoXml runInfoXml = bcl2fastq::layout::createRunInfoXml(options.getRunfolderDir());

    bcl2fastq::layout::Layout layout(
        options.getIntensitiesDir(),
        options.getInputDir(),
        options.getOutputDir(),
        options.getReportsDir(),
        options.getStatsDir(),
        options.getSampleSheet(),
        runInfoXml,
        options.getTilesFilterList(),
        options.getUseBasesMasks(),
        options.getMinimumTrimmedReadLength(),
        options.getAutoSetToZeroBarcodeMismatches(),
        options.getBarcodeMismatches(),
        options.getIgnoreMissingBcls(),
        options.getIgnoreMissingFilters(),
        options.getIgnoreMissingPositions(),
        options.getSampleSheet().getRead1StartCycle(),
        options.getSampleSheet().getRead2StartCycle(),
        options.getSampleSheet().getRead1EndCycle(),
        options.getSampleSheet().getRead2EndCycle(),
        options.getSampleSheet().getRead1UmiLength(),
        options.getSampleSheet().getRead2UmiLength(),
        options.getSampleSheet().getRead1UmiStartFromCycle(),
        options.getSampleSheet().getRead2UmiStartFromCycle(),
        options.getSampleSheet().trimUmi() == bcl2fastq::config::SampleSheetCsv::TriState::TRUE,
        options.getIncludeNonPfClusters()
    );

    bcl2fastq::conversion::ConverterStats allStats(
        layout,
        options.getStatsDir() / boost::filesystem::path( detail::fastqSummaryPartialFileName ),
        options.getStatsDir() /  boost::filesystem::path( detail::demuxSummaryPartialFileName ),
        options.getStatsDir() / boost::filesystem::path( detail::adapterStatsFileName ),
        options.getInteropDir() / boost::filesystem::path( detail::interopStatsFileName ),
        options.getStatsDir() / boost::filesystem::path( detail::demuxStatsFileName ),
        options.getStatsDir() / boost::filesystem::path( detail::conversionStatsFileName ),
        options.getStatsDir() / boost::filesystem::path( detail::statsJsonFileName )
    );

    std::for_each(
        layout.laneInfosBegin(),
        layout.laneInfosEnd(),
        boost::bind(
            detail::processLane,
            boost::ref(options),
            boost::ref(layout),
            _1,
            boost::ref(allStats)
        )
    );

    try {
        auto start = std::chrono::steady_clock::now();

        allStats.dumpInteropStats();
        allStats.dumpDemuxStats();
        allStats.dumpConversionStats();

        bcl2fastq::stats::DemuxReportGenerator(
            layout,
            options.getReportsDir()
        ).run(
            allStats.getDemuxStatsXmlPath(),
            allStats.getConversionStatsXmlPath()
        );

        BCL2FASTQ_LOG(bcl2fastq::common::LogLevel::DEBUG) << "Stats calculated in: " << std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - start).count() << "s" << std::endl;

    } catch( bcl2fastq::common::RuntimeError &e ) {
        BCL2FASTQ_LOG(bcl2fastq::common::LogLevel::WARNING) << "Could not compute and write statistics for this conversion: "
                                                            << e.what() << std::endl;
    }
}

int main(int argc, const char *argv[])
{
    bcl2fastq::common::run(program, argc, argv);
    return boost::exit_success;
}


