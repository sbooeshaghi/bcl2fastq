/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file DemuxReportGenerator.hh
 *
 * \brief Triggers the XSLT transformation
 *
 * \author Mauricio Varea
 */


#ifndef BCL2FASTQ_STATS_DEMUX_REPORT_GENERATOR_HH
#define BCL2FASTQ_STATS_DEMUX_REPORT_GENERATOR_HH

#include <boost/filesystem.hpp>

#include "layout/Layout.hh"
#include "common/Exceptions.hh"

namespace bcl2fastq
{
namespace stats
{

/// \brief Exception for Xslt errors.
class LibXsltError : public common::RuntimeError
{
public:
    LibXsltError(int errorNumber, const std::string &message)
    : RuntimeError(errorNumber, message)
    {}
    LibXsltError(const std::string &message)
    : RuntimeError(0, message)
    {}
    LibXsltError()
    : RuntimeError(EINVAL, "libxslt failure")
    {}

    LibXsltError(const LibXsltError &that)
    : RuntimeError(that)
    {}
public:
    virtual const char * what() const throw()
    {
        return this->getMessage().c_str();
    }
};


class DemuxReportGenerator
{
    const boost::filesystem::path outputDirectoryHtml_;

public:
    DemuxReportGenerator(
        const layout::Layout &flowcellLayout,
        const boost::filesystem::path &outputDirectory
    );

    void run( const boost::filesystem::path &basecallingStatsXmlPath,
              const boost::filesystem::path &demultiplexingStatsXmlPath );

    static void makeUniqueSampleNameMap(const layout::LaneInfo& laneInfo,
                                        std::map<std::string, std::string>& uniqueSampleNameMap);
};

} // namespace stats
} // namespace bcl2fastq

#endif // BCL2FASTQ_STATS_DEMUX_REPORT_GENERATOR_HH
