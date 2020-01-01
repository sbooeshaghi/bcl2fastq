/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file SampleSheetCsv.cpp
 *
 * \brief Implementation of sample sheet helper.
 *
 * \author Marek Balint
 */


#include <algorithm>
#include <utility>
#include <set>

#include <boost/format.hpp>
#include <boost/regex.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

#include "common/Logger.hh"
#include "common/Debug.hh"
#include "common/Exceptions.hh"
#include "common/DirectoryValidator.hh"
#include "config/SampleSheetCsv.hh"

namespace bcl2fastq {


namespace config {

SampleSheetCsv::ExcludedTileRange::ExcludedTileRange(const std::string& flowcell,
                                                     common::LaneNumber laneNumber,
                                                     common::TileNumber firstTile,
                                                     common::TileNumber lastTile)
: flowcell_(flowcell)
, laneNumber_(laneNumber)
, firstTile_(firstTile)
, lastTile_(lastTile)
{
}

bool SampleSheetCsv::ExcludedTileRange::isExcluded(const std::string& flowcell,
                                                   common::LaneNumber laneNumber,
                                                   common::TileNumber tileNumber) const
{
    std::string lowerFlowcell = flowcell;
    boost::algorithm::to_lower(lowerFlowcell);

    if ( (flowcell_.empty() || flowcell_ == lowerFlowcell) &&
         (laneNumber_ == 0 || laneNumber_ == laneNumber))
    {
        return (tileNumber <= lastTile_ && tileNumber >= firstTile_);
    }

    return false;
}

namespace detail {


///////////////////////////////////////////////////////////////////////////////
// section metadata
///////////////////////////////////////////////////////////////////////////////

/// \brief Sample sheet section metadata.
struct SectionMetadata
{
public:

    /// \brief Section type.
    struct SectionType
    {
        enum value_type
        {
            NA = -1,  ///< N/A
            DATA = 0, ///< [Data] section.
            SETTINGS, ///< [Settings] section.
            COUNT     ///< Total number of section types.
        };
    };

public:

    /// \brief String identifiers of individual sections.
    static const char * const sectionNames_[SectionMetadata::SectionType::COUNT];

public:

    /// \brief Section type.
    SectionType::value_type type_;

    /// \brief The first line of the section.
    common::CsvGrammarAttribute::const_iterator begin_;

    /// \brief The one-past-the-end line of the section.
    common::CsvGrammarAttribute::const_iterator end_;
};

const char * const SectionMetadata::sectionNames_[SectionMetadata::SectionType::COUNT] = {
    "Data",     // SectionType::DATA
    "Settings"  // SectionType::SETTINGS
};


///////////////////////////////////////////////////////////////////////////////
// sample sheet parsing: get individual sections
///////////////////////////////////////////////////////////////////////////////

/// \brief Find all recognized sections in sample sheet data.
/// \param sampleSheetData Source data to be parsed.
/// \param sectionInserter Insert iterator for inserting section metadata.
template<typename Inserter>
static void createSections(const common::CsvGrammarAttribute &sampleSheetData, Inserter sectionInserter)
{
    // /\[((section_name_1)|(section_name_2)|...|(section_name_N)|(\w+))\]/
    std::string sectionHeaderStr;
    sectionHeaderStr.append("\\[(");
    bool needOr = false;
    for (std::size_t i = 0; i < SectionMetadata::SectionType::COUNT; ++i)
    {
        if (needOr)
        {
            sectionHeaderStr.append("|");
        }
        needOr = true;
        sectionHeaderStr.append("(");
        sectionHeaderStr.append(SectionMetadata::sectionNames_[i]);
        sectionHeaderStr.append(")");
    }
    if (needOr)
    {
        sectionHeaderStr.append("|");
    }
    sectionHeaderStr.append("(\\w+)");
    sectionHeaderStr.append(")\\]");
    // 0: whole string; 1: section name (outermost parentheses); !2!: section_name_1; ...
    std::size_t sectionHeaderIndexBase = 2; 
    const boost::regex sectionHeaderRegex(sectionHeaderStr, boost::regex::icase);

    SectionMetadata sectionMetadata;
    sectionMetadata.type_ = SectionMetadata::SectionType::NA;
    sectionMetadata.begin_ =  sampleSheetData.begin();
    sectionMetadata.end_ =  sampleSheetData.begin();
    std::vector<int> sectionPresent(SectionMetadata::SectionType::COUNT, 0);

    BOOST_FOREACH (const common::CsvGrammarAttribute::value_type &csvLine, std::make_pair(sampleSheetData.begin(), sampleSheetData.end()))
    {
        boost::cmatch res;
        if (boost::regex_match(csvLine.front().c_str(), res, sectionHeaderRegex))
        {
            if (sectionMetadata.type_ != SectionMetadata::SectionType::NA)
            {
                *sectionInserter = sectionMetadata;
                ++sectionInserter;
                sectionMetadata.type_ = SectionMetadata::SectionType::NA;
                sectionMetadata.begin_ = sectionMetadata.end_;
            }

            // magic number: +1: this is for the last catch-all component of regexp (\w+)
            BCL2FASTQ_ASSERT_MSG(res.size() == (sectionHeaderIndexBase+SectionMetadata::SectionType::COUNT+1), "Unexpected result size: did the regular expression change?");

            for (std::size_t i = sectionHeaderIndexBase; i < (sectionHeaderIndexBase+SectionMetadata::SectionType::COUNT); ++i)
            {
                if (res[i].matched)
                {
                    BCL2FASTQ_ASSERT_MSG(sectionMetadata.type_ == SectionMetadata::SectionType::NA, "Only one section name should match: did the regular expression change?");

                    sectionMetadata.type_ = static_cast<SectionMetadata::SectionType::value_type>(i-sectionHeaderIndexBase);
                    if (sectionPresent.at(sectionMetadata.type_))
                    {
                        BOOST_THROW_EXCEPTION(common::InputDataError((boost::format("Multiple '[%s]' sections in sample sheet.") % SectionMetadata::sectionNames_[sectionMetadata.type_]).str()));
                    }
                    sectionPresent.at(sectionMetadata.type_) = true;
                }
            }
        }
        ++sectionMetadata.end_;
    }
    if (sectionMetadata.type_ != SectionMetadata::SectionType::NA)
    {
        *sectionInserter = sectionMetadata;
        ++sectionInserter;
    }
}


///////////////////////////////////////////////////////////////////////////////
// section parsing: chain of responsibility base class
///////////////////////////////////////////////////////////////////////////////

/// \brief Base class for chain of responsibility for parsing individual sample sheet sections.
class SectionParser
{
public:

    /// \brief Target storage for parsed data.
    struct Data
    {
    public:

        /// \brief Constructor.
        /// \param maskAdapters Adapters to be masked.
        /// \param trimAdapters Adapters to be trim.
        /// \param samples Samples container.
        Data(
            SampleSheetCsv::AdaptersContainer &maskAdapters,
            SampleSheetCsv::AdaptersContainer &trimAdapters,
            common::SampleMetadataContainer &samples,
            SampleSheetCsv::TriState &createFastqForIndexReads,
            SampleSheetCsv::TriState &findAdaptersWithIndels,
            SampleSheetCsv::TriState &generateReverseComplementFastqs,
            common::CycleNumber &read1EndWithCycle,
            common::CycleNumber &read2EndWithCycle,
            common::CycleNumber &read1StartFromCycle,
            common::CycleNumber &read2StartFromCycle,
            common::CycleNumber &read1UmiLength,
            common::CycleNumber &read2UmiLength,
            common::CycleNumber &read1UmiStartFromCycle,
            common::CycleNumber &read2UmiStartFromCycle,
            SampleSheetCsv::TriState &trimUmi,
            SampleSheetCsv::ExcludedTilesContainer &excludedTiles
        );

    public:

        /// \brief Adapters to be masked.
        SampleSheetCsv::AdaptersContainer &maskAdapters_;

        /// \brief Adapters to be trimmed.
        SampleSheetCsv::AdaptersContainer &trimAdapters_;

        /// \brief Sample metadata.
        common::SampleMetadataContainer &samples_;

        /// \brief Create fastq for index reads.
        SampleSheetCsv::TriState &createFastqForIndexReads_;

        /// \brief Find adapters with indels.
        SampleSheetCsv::TriState &findAdaptersWithIndels_;

        /// \brief Generate reverse complement fastqs.
        SampleSheetCsv::TriState &generateReverseComplementFastqs_;

        /// \brief End cycle number for read 1
        common::CycleNumber &read1EndWithCycle_;

        /// \brief End cycle number for read 2
        common::CycleNumber &read2EndWithCycle_;

        /// \brief Start cycle number for read 1
        common::CycleNumber &read1StartFromCycle_;

        /// \brief Start cycle number for read 2
        common::CycleNumber &read2StartFromCycle_;

        /// \brief UMI length for read 1
        common::CycleNumber &read1UmiLength_;

        /// \brief UMI length for read 2
        common::CycleNumber &read2UmiLength_;

        /// \brief UMI start cycle for read 1
        common::CycleNumber &read1UmiStartFromCycle_;

        /// \brief UMI start cycle for read 2
        common::CycleNumber &read2UmiStartFromCycle_;

        /// \brief Trim UMI if true
        SampleSheetCsv::TriState &trimUmi_;

        /// \brief Excluded tiles
        SampleSheetCsv::ExcludedTilesContainer &excludedTiles_;
    };

public:

    /// \brief Pure virtual destructor.
    virtual ~SectionParser() = 0;

public:

    /// \brief Add handler to the chain.
    void addHandler(SectionParser &sectionParser);

public:

    /// \brief Parse sample sheet section.
    /// \param target Target data storage for storing result.
    /// \param source Source data to be parsed.
    /// \retval true Data has been parsed.
    /// \retval false Data has not been parsed by any of handlers.
    virtual bool parse(SectionParser::Data &target, const SectionMetadata &source) = 0;

private:

    /// \brief Next handler in the chain.
    SectionParser *next_;
};

SectionParser::~SectionParser()
{
}

void SectionParser::addHandler(SectionParser &sectionParser)
{
    if (next_)
    {
        next_->addHandler(sectionParser);
    }
    else
    {
        next_ = &sectionParser;
    }
}

bool SectionParser::parse(SectionParser::Data &target, const SectionMetadata &source)
{
    // default behavior: pass to the next handler (if any)
    if (next_)
    {
        return next_->parse(target, source);
    }
    return false;
}

SectionParser::Data::Data(
    SampleSheetCsv::AdaptersContainer &maskAdapters,
    SampleSheetCsv::AdaptersContainer &trimAdapters,
    common::SampleMetadataContainer &samples,
    SampleSheetCsv::TriState &createFastqForIndexReads,
    SampleSheetCsv::TriState &findAdaptersWithIndels,
    SampleSheetCsv::TriState &generateReverseComplementFastqs,
    common::CycleNumber &read1EndWithCycle,
    common::CycleNumber &read2EndWithCycle,
    common::CycleNumber &read1StartFromCycle,
    common::CycleNumber &read2StartFromCycle,
    common::CycleNumber &read1UmiLength,
    common::CycleNumber &read2UmiLength,
    common::CycleNumber &read1UmiStartFromCycle,
    common::CycleNumber &read2UmiStartFromCycle,
    SampleSheetCsv::TriState &trimUmi,
    SampleSheetCsv::ExcludedTilesContainer& excludedTiles
)
: maskAdapters_(maskAdapters)
, trimAdapters_(trimAdapters)
, samples_(samples)
, createFastqForIndexReads_(createFastqForIndexReads)
, findAdaptersWithIndels_(findAdaptersWithIndels)
, generateReverseComplementFastqs_(generateReverseComplementFastqs)
, read1EndWithCycle_(read1EndWithCycle)
, read2EndWithCycle_(read2EndWithCycle)
, read1StartFromCycle_(read1StartFromCycle)
, read2StartFromCycle_(read2StartFromCycle)
, read1UmiLength_(read1UmiLength)
, read2UmiLength_(read2UmiLength)
, read1UmiStartFromCycle_(read1UmiStartFromCycle)
, read2UmiStartFromCycle_(read2UmiStartFromCycle)
, trimUmi_(trimUmi)
, excludedTiles_(excludedTiles)
{
}


///////////////////////////////////////////////////////////////////////////////
// section parsing: [Data]
///////////////////////////////////////////////////////////////////////////////

/// \brief Concrete handler for chain of responsibility for parsing sample sheet data section.
class DataSectionParser : public SectionParser
{
public:

    /// \brief Column in the samples data table.
    struct Column
    {
        enum value_type
        {
            NA = -1,       ///< N/A
            SAMPLE_ID = 0, ///< Sample ID.
            SAMPLE_NAME,   ///< Sample name.
            PROJECT,       ///< Project name.
            LANE,          ///< Lane number.
            BARCODE_1,     ///< Barcode #1.
            BARCODE_2,     ///< Barcode #2.
            COUNT          ///< Total number of recognized columns.
        };
    };

public:

    /// \brief String identifiers of individual columns.
    static const std::string columnNames_[DataSectionParser::Column::COUNT][3];

public:

    /// \brief Parse column string identifier.
    /// \param columnName Column string identifier.
    /// \return Column numeric unique identifier.
    static Column::value_type identifyColumnName(std::string columnName);

public:
    virtual ~DataSectionParser();

public:

    virtual bool parse(SectionParser::Data &target, const SectionMetadata &source);

private:

    virtual void parseLanes(const std::string&      lanes,
                            common::SampleMetadata& sampleMetadata) const;

    virtual void validateCharacters(const std::string& value,
                                    const std::string& columnName) const;

    virtual void addBarcodes(common::SampleMetadata&                                 sampleMetadata,
                             std::string&                                            barcode1,
                             std::string&                                            barcode2,
                             std::set< std::pair<common::LaneNumber, std::string> >& uniqueBarcodes);
};

const std::string DataSectionParser::columnNames_[DataSectionParser::Column::COUNT][3] = {
    { "sampleid"   , "sample_id"     , "" }, // Column::SAMPLE_ID
    { "samplename" , "sample_name"   , "" }, // Column::SAMPLE_NAME
    { "project"    , "sample_project", "" }, // Column::PROJECT
    { "lane"       , "lanes"         , "" }, // Column::LANE
    { "index"      , ""              , "" }, // Column::BARCODE_1
    { "index2"     , ""              , "" }, // Column::BARCODE_2
};

DataSectionParser::Column::value_type DataSectionParser::identifyColumnName(std::string columnName)
{
    boost::algorithm::to_lower(columnName);
    for (std::size_t column = 0; column < Column::COUNT; ++column)
    {
        for (std::size_t i = 0; columnNames_[column][i][0] != '\0'; ++i)
        {
            if (columnName == columnNames_[column][i])
            {
                return static_cast<DataSectionParser::Column::value_type>(column);
            }
        }
    }
    return Column::NA;
}

DataSectionParser::~DataSectionParser()
{
}

bool DataSectionParser::parse(SectionParser::Data &target, const SectionMetadata &source)
{
    if (source.type_ != SectionMetadata::SectionType::DATA)
    {
        return SectionParser::parse(target, source);
    }

    const std::size_t columnNotPresent = -1;
    std::vector<std::size_t> colsIndicies(Column::COUNT, columnNotPresent);
    bool colsIndiciesInitialized = false;

    std::set< std::pair<common::LaneNumber, std::string> > uniqueBarcodes;

    BOOST_FOREACH (const common::CsvGrammarAttribute::value_type &line, std::make_pair(source.begin_, source.end_))
    {
        if (!colsIndiciesInitialized)
        {
            size_t colIdx = 0;
            BOOST_FOREACH (const common::CsvGrammarAttribute::value_type::value_type &value, std::make_pair(line.begin(), line.end()))
            {
                Column::value_type colId = identifyColumnName(value);
                if (colId != Column::NA)
                {
                    colsIndicies.at(colId) = colIdx;
                    colsIndiciesInitialized = true;
                }
                ++colIdx;
            }
        }
        else
        {
            std::vector<std::string> values(Column::COUNT, "");
            for (std::size_t colId = 0; colId < Column::COUNT; ++colId)
            {
                if (colsIndicies.at(colId) != columnNotPresent)
                {
                    values.at(colId) = line.at(colsIndicies.at(colId));
                    boost::algorithm::trim(values.at(colId));
                }
            }

            if (values.at(Column::SAMPLE_ID).empty() && values.at(Column::SAMPLE_NAME).empty())
            {
                // skip empty lines
                continue;
            }

            // new sample
            target.samples_.push_back(common::SampleMetadata());
            common::SampleMetadata &sampleMetadata = target.samples_.back();
            sampleMetadata.id_ = values.at(Column::SAMPLE_ID);
            sampleMetadata.name_ = values.at(Column::SAMPLE_NAME);
            sampleMetadata.project_ = values.at(Column::PROJECT);

            validateCharacters(sampleMetadata.id_, columnNames_[Column::SAMPLE_ID][0]);
            validateCharacters(sampleMetadata.name_, columnNames_[Column::SAMPLE_NAME][0]);
            validateCharacters(sampleMetadata.project_, columnNames_[Column::PROJECT][0]);

            parseLanes(values.at(Column::LANE),
                       sampleMetadata);

            std::string& barcode1 = values.at(Column::BARCODE_1);
            std::string& barcode2 = values.at(Column::BARCODE_2);
            if (!barcode1.empty() || !barcode2.empty())
            {
                addBarcodes(sampleMetadata,
                            barcode1,
                            barcode2,
                            uniqueBarcodes);
            }
        }
    }
    return true;
}

void DataSectionParser::parseLanes(const std::string&      lanes,
                                   common::SampleMetadata& sampleMetadata) const
{
    if (lanes.empty())
    {
        return;
    }

    std::vector<std::string> sampleLaneStrings;
    boost::split(sampleLaneStrings, lanes, boost::is_any_of("+"));
    for (std::string& lane : sampleLaneStrings)
    {
        boost::algorithm::trim(lane);
        try
        {
            sampleMetadata.lanes_.push_back(boost::lexical_cast<common::LaneNumber>(lane));
        }
        catch (boost::bad_lexical_cast&)
        {
            BOOST_THROW_EXCEPTION(common::InputDataError("Invalid character in lane: '" + lane + "'."));
        }
    }
}

void DataSectionParser::validateCharacters(const std::string& value,
                                           const std::string& columnName) const
{
    if (!value.empty() && !boost::regex_match(value, boost::regex("[a-zA-Z0-9_-]+")))
    {
        BOOST_THROW_EXCEPTION(common::InputDataError("Invalid characters in " + columnName +
                              ": '" + value + "'. Only alphanumeric, '-', or '_' characters are allowed."));
    }
}

void DataSectionParser::addBarcodes(common::SampleMetadata&                                 sampleMetadata,
                                    std::string&                                            barcode1,
                                    std::string&                                            barcode2,
                                    std::set< std::pair<common::LaneNumber, std::string> >& uniqueBarcodes)
{
    if (barcode1.find_first_not_of("ACGTN") != std::string::npos)
    {
        BOOST_THROW_EXCEPTION(common::InputDataError((boost::format("Error parsing barcode '%s' in sample sheet.") % barcode1).str()));
    }
    if (barcode2.find_first_not_of("ACGTN") != std::string::npos)
    {
        BOOST_THROW_EXCEPTION(common::InputDataError((boost::format("Error parsing barcode '%s' in sample sheet.") % barcode2).str()));
    }

    sampleMetadata.barcodes_.resize(sampleMetadata.barcodes_.size()+1);
    sampleMetadata.barcodes_.back().push_back(barcode1);
    sampleMetadata.barcodes_.back().push_back(barcode2);

    for (common::LaneNumber laneNumber : sampleMetadata.lanes_)
    {
        std::pair< std::set< std::pair<common::LaneNumber, std::string> >::iterator, bool > insertPos = uniqueBarcodes.insert(std::make_pair(
            laneNumber, barcode1 + ((barcode1.empty() || barcode2.empty() ? "" : "-") + barcode2) ));
    
        if (!insertPos.second)
        {        
            BOOST_THROW_EXCEPTION(common::InputDataError((boost::format("Duplicate barcode '%s' found in sample sheet.")
                % insertPos.first->second).str()));
        }
    }
}


///////////////////////////////////////////////////////////////////////////////
// section parsing: [Settings]
///////////////////////////////////////////////////////////////////////////////

/// \brief Concrete handler for chain of responsibility for parsing sample sheet settings section.
class SettingsSectionParser : public SectionParser
{
public:

    virtual ~SettingsSectionParser();

public:

    virtual bool parse(SectionParser::Data &target, const SectionMetadata &source);

private:

    virtual void handleExcludedTiles(const std::string& key,
                                     const std::string& value,
                                     SectionParser::Data &target) const;

    virtual SampleSheetCsv::TriState getBoolSetting(const std::string& key, const std::string& value) const;

    virtual size_t getIntSetting(const std::string& value) const;

    void addAdapters(SampleSheetCsv::AdaptersContainer &adapterContainer,
                     common::ReadNumber                 readNumber,
                     const std::string                 &adapters) const;
};

SettingsSectionParser::~SettingsSectionParser()
{
}

static const std::string excludeTiles("excludetiles");
static const std::string excludeTilesLane("excludetileslane");
static const std::string flowcellStr("flowcell");

bool SettingsSectionParser::parse(SectionParser::Data &target, const SectionMetadata &source)
{
    if (source.type_ != SectionMetadata::SectionType::SETTINGS)
    {
        return SectionParser::parse(target, source);
    }

    static const std::string maskAdapter("maskadapter");
    static const std::string maskAdapterRead("maskadapterread");
    static const std::string trimAdapter("trimadapter");
    static const std::string trimAdapterDeprecated("adapter");
    static const std::string trimAdapterRead("trimadapterread");
    static const std::string trimAdapterReadDeprecated("adapterread");

    static const std::string createFastqForIndexReads("createfastqforindexreads");
    static const std::string findAdapterWithIndelsDeprecated("findadapterwithindels");
    static const std::string findAdaptersWithIndels("findadapterswithindels");
    static const std::string generateReverseComplementFastqs("reversecomplement");

    static const std::string read1EndWithCycle("read1endwithcycle");
    static const std::string read2EndWithCycle("read2endwithcycle");
    static const std::string read1StartFromCycle("read1startfromcycle");
    static const std::string read2StartFromCycle("read2startfromcycle");
    static const std::string read1UmiLength("read1umilength");
    static const std::string read2UmiLength("read2umilength");
    static const std::string read1UmiStartFromCycle("read1umistartfromcycle");
    static const std::string read2UmiStartFromCycle("read2umistartfromcycle");
    static const std::string trimUmi("trimumi");

    BOOST_FOREACH (const common::CsvGrammarAttribute::value_type &line, std::make_pair(source.begin_, source.end_))
    {
        std::string key = line.at(0);
        boost::algorithm::to_lower(key);
        boost::algorithm::trim(key);

        std::string value = line.at(1);
        boost::algorithm::trim(value);

        if (key == maskAdapter)
        {
            addAdapters(target.maskAdapters_,
                        0,
                        value);
        }
        else if ((key == trimAdapter) || (key == trimAdapterDeprecated))
        {
            addAdapters(target.trimAdapters_,
                        0,
                        value);
        }
        else if ((key.find(maskAdapterRead) == 0) && (key.size() == maskAdapterRead.size()+1))
        {
            const char readKey = key.at(key.size()-1);
            if (readKey >= '1' && readKey <= '9')
            {
                addAdapters(target.maskAdapters_,
                            readKey-'0',
                            value);
            }
        }
        else if (((key.find(trimAdapterRead) == 0) && (key.size() == trimAdapterRead.size()+1)) || ((key.find(trimAdapterReadDeprecated) == 0) && (key.size() == trimAdapterReadDeprecated.size()+1)))
        {
            const char readKey = key.at(key.size()-1);
            if (readKey >= '1' && readKey <= '9')
            {
                addAdapters(target.trimAdapters_,
                            readKey-'0',
                            value);
            }
        }
        else if (key == createFastqForIndexReads)
        {
            target.createFastqForIndexReads_ = getBoolSetting(key, value);
        }
        else if (key == findAdaptersWithIndels || key == findAdapterWithIndelsDeprecated)
        {
            if (key == findAdapterWithIndelsDeprecated) {
                BCL2FASTQ_LOG(common::LogLevel::WARNING) 
                    << "Option '" << findAdapterWithIndelsDeprecated 
                    << "' is deprecated. Please use '" 
                    << findAdaptersWithIndels << "'." << std::endl;
            }
            target.findAdaptersWithIndels_ = getBoolSetting(key, value);
        }
        else if (key == generateReverseComplementFastqs)
        {
            target.generateReverseComplementFastqs_ = getBoolSetting(key, value);
        }
        else if (key == read1EndWithCycle)
        {
            target.read1EndWithCycle_ = getIntSetting(value);
        }
        else if (key == read2EndWithCycle)
        {
            target.read2EndWithCycle_ = getIntSetting(value);
        }
        else if (key == read1StartFromCycle)
        {
            target.read1StartFromCycle_ = getIntSetting(value);
        }
        else if (key == read2StartFromCycle)
        {
            target.read2StartFromCycle_ = getIntSetting(value);
        }
        else if (key == read1UmiLength)
        {
            target.read1UmiLength_ = getIntSetting(value);
        }
        else if (key == read2UmiLength)
        {
            target.read2UmiLength_ = getIntSetting(value);
        }
        else if (key == read1UmiStartFromCycle)
        {
            target.read1UmiStartFromCycle_ = getIntSetting(value);
        }
        else if (key == read2UmiStartFromCycle)
        {
            target.read2UmiStartFromCycle_ = getIntSetting(value);
        }
        else if (key == trimUmi)
        {
            target.trimUmi_ = getBoolSetting(key, value);
        }
        else if (key.substr(0, excludeTiles.size()) == excludeTiles)
        {
            handleExcludedTiles(key, value, target);
        }
    }

    return true;
}

void SettingsSectionParser::handleExcludedTiles(const std::string&   key,
                                                const std::string&   value,
                                                SectionParser::Data& target) const
{
    std::vector<std::string> tileRanges;
    boost::split(tileRanges, value, boost::is_any_of("+,"));

    common::LaneNumber laneNumber = 0;
    std::string flowcell;
    if (key.substr(0, excludeTilesLane.size()) == excludeTilesLane)
    {
        BCL2FASTQ_ASSERT_MSG(key.size() > excludeTilesLane.size(),
                             "Config parameter: ExcludedTilesLane must include the lane number. (e.g. ExcludedTilesLane5)");

        size_t flowcellPos = key.find(flowcellStr, excludeTilesLane.size());

        if (flowcellPos != std::string::npos)
        {
            BCL2FASTQ_ASSERT_MSG(key.size() > flowcellPos + flowcell.size(),
                "Config parameter ExcludeTilesLaneXFlowcell must include the flowcell id. (e.g. ExcludeTilesLane5FlowcellABCD");

            laneNumber = boost::lexical_cast<common::LaneNumber>(
                key.substr(excludeTilesLane.size(), flowcellPos-excludeTilesLane.size()));

            flowcell = key.substr(flowcellPos + flowcellStr.size());
        }
        else
        {
            laneNumber = boost::lexical_cast<common::LaneNumber>(key.substr(excludeTilesLane.size()));
        }
    }

    BOOST_FOREACH(const std::string& tileRange, tileRanges)
    {
        std::vector<std::string> tileRangeVec;
        boost::split(tileRangeVec, tileRange, boost::is_any_of("-"));

        BCL2FASTQ_ASSERT_MSG(tileRangeVec.size() == 1 || tileRangeVec.size() == 2, "Formatting error in ExcludedTiles config parameter");

        common::TileNumber firstTile = boost::lexical_cast<common::TileNumber>(tileRangeVec[0]);
        common::TileNumber lastTile = tileRangeVec.size() == 2 ?
            boost::lexical_cast<common::TileNumber>(tileRangeVec[1]) :
            firstTile;

        target.excludedTiles_.push_back(
            SampleSheetCsv::ExcludedTileRange(flowcell,
                                              laneNumber,
                                              firstTile,
                                              lastTile));
    }
}

SampleSheetCsv::TriState SettingsSectionParser::getBoolSetting(const std::string& key, const std::string& value) const
{
    static const std::string trueStrings[] = { "true", "t", "yes", "y", "1" };
    static const std::string falseStrings[] = { "false", "f", "no", "n", "0" };
    static const std::set<std::string> trueSet(trueStrings, trueStrings + sizeof(trueStrings)/sizeof(trueStrings[0]));
    static const std::set<std::string> falseSet(falseStrings, falseStrings + sizeof(falseStrings)/sizeof(falseStrings[0]));

    std::string valueLowerCase = value;
    boost::algorithm::to_lower(valueLowerCase);

    if (trueSet.find(valueLowerCase) != trueSet.end())
    {
        return SampleSheetCsv::TRUE;
    }

    if (falseSet.find(valueLowerCase) != falseSet.end())
    {
        return SampleSheetCsv::FALSE;
    }

    BOOST_THROW_EXCEPTION(common::InputDataError("Error parsing sample sheet setting: '" + key +
        "'. Expected bool (true, false, t, f, yes, no, y, n, 1, 0), got: '" + value + "'"));

    return SampleSheetCsv::INDETERMINATE;
}

size_t SettingsSectionParser::getIntSetting(const std::string& value) const
{
    return boost::lexical_cast<size_t>(value);
}

void SettingsSectionParser::addAdapters(SampleSheetCsv::AdaptersContainer &adapterContainer,
                                        common::ReadNumber                 readNumber,
                                        const std::string                 &adapters) const
{
    std::vector<std::string> adapterList;
    boost::split(adapterList, adapters, boost::is_any_of("+"));

    BOOST_FOREACH(std::string &adapter, adapterList)
    {
        boost::algorithm::trim(adapter);
        if (adapter.find_first_not_of("ACGT") != std::string::npos)
        {
            BOOST_THROW_EXCEPTION(common::InputDataError((boost::format("Error parsing adapter '%s' in sample sheet.") % adapter).str()));
        }

        if (!adapter.empty())
        {
            adapterContainer.push_back(std::make_pair(readNumber, adapter));
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
// section parsing: create chain of responsibility
///////////////////////////////////////////////////////////////////////////////


/// \brief Initialize section parser.
/// \return Section parser.
static SectionParser * createSectionParser()
{
    static DataSectionParser dataSectionParser;
    static SettingsSectionParser settingsSectionParser;

    dataSectionParser.addHandler(settingsSectionParser);
    return &dataSectionParser;
}

/// \brief Section parser.
static SectionParser * sectionParser = createSectionParser();


} // namespace detail


///////////////////////////////////////////////////////////////////////////////
// SampleSheetCsv
///////////////////////////////////////////////////////////////////////////////


SampleSheetCsv::SampleSheetCsv(const common::CsvGrammarAttribute& sampleSheetData,
                               const boost::filesystem::path&     outputDir)
: maskAdapters_()
, trimAdapters_()
, samples_()
, createFastqForIndexReads_(SampleSheetCsv::INDETERMINATE)
, findAdaptersWithIndels_(SampleSheetCsv::INDETERMINATE)
, generateReverseComplementFastqs_(SampleSheetCsv::INDETERMINATE)
, read1EndWithCycle_(0)
, read2EndWithCycle_(0)
, read1StartFromCycle_(0)
, read2StartFromCycle_(0)
, read1UmiLength_(0)
, read2UmiLength_(0)
, read1UmiStartFromCycle_(1)
, read2UmiStartFromCycle_(1)
, trimUmi_(SampleSheetCsv::INDETERMINATE)
, excludedTiles_()
{
    typedef std::vector<detail::SectionMetadata> SectionsContainer;
    SectionsContainer sectionsContainer;
    detail::createSections(sampleSheetData, std::back_inserter(sectionsContainer));

    detail::SectionParser::Data sectionDataTarget(maskAdapters_,
                                                  trimAdapters_,
                                                  samples_,
                                                  createFastqForIndexReads_,
                                                  findAdaptersWithIndels_,
                                                  generateReverseComplementFastqs_,
                                                  read1EndWithCycle_,
                                                  read2EndWithCycle_,
                                                  read1StartFromCycle_,
                                                  read2StartFromCycle_,
                                                  read1UmiLength_,
                                                  read2UmiLength_,
                                                  read1UmiStartFromCycle_,
                                                  read2UmiStartFromCycle_,
                                                  trimUmi_,
                                                  excludedTiles_);

    BOOST_FOREACH (const detail::SectionMetadata &sectionDataSource, sectionsContainer)
    {
        detail::sectionParser->parse(sectionDataTarget, sectionDataSource);

        BOOST_FOREACH(common::SampleMetadata& sampleMetadata,
                      std::make_pair(sectionDataTarget.samples_.begin(),
                                     sectionDataTarget.samples_.end()))
        {
            common::DirectoryValidator::getSingleton().validateNoClash(outputDir / sampleMetadata.project_, "Project");
        }
    }
}

SampleSheetCsv::AdaptersContainer::const_iterator SampleSheetCsv::maskAdaptersBegin() const
{
    return maskAdapters_.begin();
}

SampleSheetCsv::AdaptersContainer::const_iterator SampleSheetCsv::maskAdaptersEnd() const
{
    return maskAdapters_.end();
}

SampleSheetCsv::AdaptersContainer::const_iterator SampleSheetCsv::trimAdaptersBegin() const
{
    return trimAdapters_.begin();
}

SampleSheetCsv::AdaptersContainer::const_iterator SampleSheetCsv::trimAdaptersEnd() const
{
    return trimAdapters_.end();
}

bool SampleSheetCsv::isTileExcluded(const std::string& flowcell,
                                    common::LaneNumber laneNumber,
                                    common::TileNumber tileNumber) const
{
    BOOST_FOREACH(const SampleSheetCsv::ExcludedTileRange& excludedTile, excludedTiles_)
    {
        if (excludedTile.isExcluded(flowcell, laneNumber, tileNumber)) return true;
    }

    return false;
}

SampleSheetCsvPtr createSampleSheetCsv(const boost::filesystem::path& sampleSheet,
                                       const boost::filesystem::path& outputDir)
{
    common::CsvGrammarAttribute sampleSheetData;

    BCL2FASTQ_LOG(common::LogLevel::INFO) << "Sample sheet: '" << sampleSheet.string() << "'" << std::endl;

    if (boost::filesystem::exists(sampleSheet))
    {
        sampleSheetData = common::parseCsvFile(sampleSheet);
    }
    else
    {
        BCL2FASTQ_LOG(common::LogLevel::INFO) << "Sample sheet: NOT FOUND" << std::endl;
    }

    return SampleSheetCsvPtr( new SampleSheetCsv(sampleSheetData,
                                                 outputDir));
}


} // namespace config


} // namespace bcl2fastq


