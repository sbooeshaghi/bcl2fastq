/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file CsvGrammar.cpp
 *
 * \brief Basic CSV grammar.
 *
 * \author Marek Balint
 */


#include <algorithm>
#include <iterator>

#include <boost/format.hpp>
#include <boost/bind.hpp>
#include <boost/filesystem/fstream.hpp>

#include "common/CsvGrammar.hh"
#include "common/Exceptions.hh"


namespace bcl2fastq {


namespace common {


CsvGrammarAttribute parseCsvFile(const boost::filesystem::path& csvFile)
{
    int errnum;

    boost::filesystem::ifstream csvStream(csvFile);
    errnum = errno;
    if (!csvStream)
    {
        BOOST_THROW_EXCEPTION(common::IoError(errnum, (boost::format("Unable to open '%s' file for reading") % csvFile.string()).str()));
    }
    csvStream.unsetf(std::ios::skipws);

    typedef std::vector<char> CsvBuffer;
    CsvBuffer csvBuffer;
    errnum = 0;
    std::copy(
        std::istream_iterator<CsvBuffer::value_type>(csvStream),
        std::istream_iterator<CsvBuffer::value_type>(),
        std::back_inserter(csvBuffer)
    );
    errnum = errno;
    if (!csvStream.eof())
    {
        BOOST_THROW_EXCEPTION(common::IoError(errnum, "Failed to read " + csvFile.string() + " to the end"));
    }

    return parseCsvData(csvBuffer.begin(), csvBuffer.end());
}

CsvGrammarAttribute::value_type::size_type correctCsvColumnsCount(CsvGrammarAttribute &csvData, const std::string &defaultValue)
{
    if (csvData.empty())
    {
        return 0;
    }

    const common::CsvGrammarAttribute::value_type::size_type maxCols = std::max_element(
        csvData.begin(),
        csvData.end(), 
        boost::bind(&common::CsvGrammarAttribute::value_type::size, _1) < boost::bind(&common::CsvGrammarAttribute::value_type::size, _2)
    )->size();
    std::for_each(
        csvData.begin(),
        csvData.end(), 
        boost::bind(&common::CsvGrammarAttribute::value_type::resize, _1, maxCols, defaultValue)
    );

    return maxCols;
}

} // namespace common


} // namespace bcl2fastq


