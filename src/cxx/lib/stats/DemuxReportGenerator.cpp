/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file DemuxReportGenerator.cpp
 *
 * \brief Stats
 *
 * \author Mauricio Varea
 */

#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/assign/list_of.hpp>

#include "config.h"

#include <libxslt/xslt.h>
#include <libxslt/xsltInternals.h>
#include <libxslt/transform.h>
#include <libxslt/xsltutils.h>
#include <libexslt/exslt.h>

#include "common/Debug.hh"
#include "common/FileSystem.hh"
#include "common/InstallationPath.hh"
#include "layout/Layout.hh"

#include "stats/DemuxReportGenerator.hh"

#include <set>

#ifndef WIN32
extern int xmlLoadExtDtdDefaultValue;
#endif

namespace bcl2fastq
{
namespace stats
{


DemuxReportGenerator::DemuxReportGenerator( const layout::Layout &flowcellLayout,
                                            const boost::filesystem::path &outputDirectory )
: outputDirectoryHtml_(outputDirectory/"html")
{
    std::vector<boost::filesystem::path> createList =
                boost::assign::list_of(outputDirectory)(outputDirectoryHtml_);
    const boost::filesystem::path flowcellHtmlPath = outputDirectoryHtml_ / flowcellLayout.getFlowcellInfo().getFlowcellId();
    createList.push_back(flowcellHtmlPath);
    createList.push_back(flowcellHtmlPath/"all");
    createList.push_back(flowcellHtmlPath/"all"/"all");
    createList.push_back(flowcellHtmlPath/"all"/"all"/"all");
    createList.push_back(flowcellHtmlPath/"all"/"all"/std::string(layout::Barcode::DEFAULT_BARCODE));

    BOOST_FOREACH (const layout::LaneInfo &laneInfo, std::make_pair(flowcellLayout.laneInfosBegin(), flowcellLayout.laneInfosEnd()))
    {
        std::map<std::string, std::string> uniqueSampleNameMap;
        makeUniqueSampleNameMap(laneInfo,
                                uniqueSampleNameMap);

        BOOST_FOREACH (const layout::SampleInfo &sampleInfo, std::make_pair(laneInfo.sampleInfosBegin(), laneInfo.sampleInfosEnd()))
        {
            const std::string sampleName = uniqueSampleNameMap[sampleInfo.getId()];

            /* No need to worry about repeated directories, due to more than one sample belonging to
               the same project, as createDirectories() should take care of this */
            layout::SampleInfo::Project projectName = sampleInfo.getProject();
            createList.push_back(flowcellHtmlPath/projectName);
            createList.push_back(flowcellHtmlPath/projectName/"all");
            createList.push_back(flowcellHtmlPath/projectName/"all"/"all");
            createList.push_back(flowcellHtmlPath/projectName/"all"/std::string(layout::Barcode::DEFAULT_BARCODE));
            createList.push_back(flowcellHtmlPath/projectName/sampleName);
            createList.push_back(flowcellHtmlPath/projectName/sampleName/"all");
            if ( sampleInfo.hasBarcodes() )
            {
                for (const auto& barcode : sampleInfo.getBarcodes())
                {
                    createList.push_back(flowcellHtmlPath/projectName/sampleName/barcode.toString());
                }
            } //else {
                createList.push_back(flowcellHtmlPath/projectName/sampleName/layout::Barcode::DEFAULT_BARCODE);
            //}
        }
    }
    common::createDirectories(createList);
}


void DemuxReportGenerator::run( const boost::filesystem::path &demultiplexingStatsXmlPath,
                                const boost::filesystem::path &conversionStatsXmlPath )
{
    const char *params[] = {"OUTPUT_DIRECTORY_HTML_PARAM", "''",
                            "DEMULTIPLEXING_STATS_XML_PARAM", "''",
                            "BCL2FASTQ_FULL_DATADIR_PARAM", "''",
                            0};
    std::string quotedOutputHtmlDirectory = "'" + outputDirectoryHtml_.string() + "'";

    std::string quotedDemultiplexingStatsXml= "'" + demultiplexingStatsXmlPath.string() + "'";

    const boost::filesystem::path fullDataDir(common::InstallationPath::getSingleton().expandPath(BCL2FASTQ_DATADIR));
    std::string quotedFullDataDirPath = "'" + fullDataDir.string() + "'";

    boost::filesystem::path reportGenerator(fullDataDir / "xsl" / "demux" / "GenerateReport.xsl");

    // When executing under cygwin, generic_string must be called to get the forward slashes.
    // Otherwise, xsltParseStylesheetFile will fail with cryptic error messages.
    std::string reportGeneratorStr = reportGenerator.string();

#ifdef WIN32
    std::replace(quotedOutputHtmlDirectory.begin(), quotedOutputHtmlDirectory.end(), '\\', '/');
    std::replace(quotedDemultiplexingStatsXml.begin(), quotedDemultiplexingStatsXml.end(), '\\', '/');
    std::replace(quotedFullDataDirPath.begin(), quotedFullDataDirPath.end(), '\\', '/');
    std::replace(reportGeneratorStr.begin(), reportGeneratorStr.end(), '\\', '/');
#endif

    params[1] = quotedOutputHtmlDirectory.c_str();
    params[3] = quotedDemultiplexingStatsXml.c_str();
    params[5] = quotedFullDataDirPath.c_str();

    xmlSubstituteEntitiesDefault(1);
    xmlLoadExtDtdDefaultValue = 1;

    exsltRegisterAll();

    xsltStylesheetPtr cur = xsltParseStylesheetFile((const xmlChar *) reportGeneratorStr.c_str() );
    xmlDocPtr doc = xmlParseFile(conversionStatsXmlPath.string().c_str());
    xmlDocPtr res = xsltApplyStylesheet(cur, doc, params);

    BCL2FASTQ_LOG(common::LogLevel::DEBUG) << "ReportGenerator: Applying stylesheet " << reportGenerator << std::endl;
    for(unsigned int i=0; i<( (sizeof(params)/sizeof(const char *))-2 ); i+=2)
    {
        BCL2FASTQ_LOG(common::LogLevel::DEBUG) << "ReportGenerator: --stringparam " << std::string(params[i]) << " " << std::string(params[i+1]) << std::endl;
    }
    BCL2FASTQ_LOG(common::LogLevel::DEBUG) << "ReportGenerator: Using data from " << conversionStatsXmlPath << std::endl;

    xsltFreeStylesheet(cur);
    xmlFreeDoc(res);
    xmlFreeDoc(doc);

    xsltCleanupGlobals();
    xmlCleanupParser();
    if (!res)
    {
        BOOST_THROW_EXCEPTION(LibXsltError());
    }
}

void DemuxReportGenerator::makeUniqueSampleNameMap(const layout::LaneInfo& laneInfo,
                                                   std::map<std::string, std::string>& uniqueSampleNameMap)
{
    std::set<std::string> sampleNameSet;
    std::set<std::string> duplicateNames;
    BOOST_FOREACH(const layout::SampleInfo& sampleInfo, std::make_pair(laneInfo.sampleInfosBegin(), laneInfo.sampleInfosEnd()))
    {
        if (!sampleNameSet.insert(sampleInfo.getName()).second)
        {
            duplicateNames.insert(sampleInfo.getName());
        }
    }

    BOOST_FOREACH(const layout::SampleInfo& sampleInfo, std::make_pair(laneInfo.sampleInfosBegin(), laneInfo.sampleInfosEnd()))
    {
        if (duplicateNames.count(sampleInfo.getName()) != 0)
        {
            uniqueSampleNameMap.insert(
                std::make_pair(sampleInfo.getId(),
                               sampleInfo.getName() + "_S" + boost::lexical_cast<std::string>(sampleInfo.getNumber())));
        }
        else
        {
            uniqueSampleNameMap.insert(std::make_pair(sampleInfo.getId(), sampleInfo.getName()));
        }
    }
}


} // namespace stats
} // namespace bcl2fastq
