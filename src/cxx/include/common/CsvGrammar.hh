/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file CsvGrammar.hh
 *
 * \brief Basic CSV grammar.
 *
 * \author Roman Petrovski
 * \author Marek Balint
 */


#ifndef BCL2FASTQ_COMMON_CSV_GRAMMAR_HH
#define BCL2FASTQ_COMMON_CSV_GRAMMAR_HH


#include <vector>
#include <string>

#include <boost/filesystem/path.hpp>

#include <boost/spirit/include/qi_grammar.hpp>
#include <boost/spirit/include/qi_rule.hpp>
#include <boost/spirit/home/support/unused.hpp>


namespace bcl2fastq {


namespace common {


/// \brief CSV grammar attribute.
typedef std::vector<std::vector<std::string> > CsvGrammarAttribute;

/// \brief CSV grammar.
template <typename Iterator>
struct CsvGrammar : public boost::spirit::qi::grammar<Iterator, CsvGrammarAttribute()>
{
public:

    /// \brief Default constructor.
    CsvGrammar();

public:

    /// \brief CSV file.
    boost::spirit::qi::rule<Iterator, CsvGrammarAttribute()> start_;

    /// \brief CSV line.
    boost::spirit::qi::rule<Iterator, std::vector<std::string>() > line_;

    /// \brief Csv data line.
    boost::spirit::qi::rule<Iterator, std::vector<std::string>() > dataLine_;

    /// \brief CSV comment line.
    boost::spirit::qi::rule<Iterator, boost::spirit::unused_type()> commentLine_;

    /// \brief CSV field.
    boost::spirit::qi::rule<Iterator, std::string()> field_;

    /// \brief Escaped (quoted) field value.
    boost::spirit::qi::rule<Iterator, std::string()> escaped_;

    /// \brief Non-escaped (not quoted) field value.
    boost::spirit::qi::rule<Iterator, std::string()> nonEscaped_;

    /// \brief All printable characters plus space, except dquote_ and comma_.
    boost::spirit::qi::rule<Iterator, char()> textData_;

    /// \brief CR.
    boost::spirit::qi::rule<Iterator, char()> cr_;

    /// \brief LF.
    boost::spirit::qi::rule<Iterator, char()> lf_;

    /// \brief Line ending.
    boost::spirit::qi::rule<Iterator, boost::spirit::unused_type()> crlf_;

    /// \brief Comma.
    boost::spirit::qi::rule<Iterator, char()> comma_;

    /// \brief Double-quote.
    boost::spirit::qi::rule<Iterator, char()> dquote_;
};


/// \brief Parse CSV data.
/// \param csvDataBegin Begin iterator to CSV data to be parsed.
/// \param csvDataEnd End iterator to CSV data to be parsed.
/// \return Data strucutre containing data parsed out of the CSV data.
template<typename Iterator>
CsvGrammarAttribute parseCsvData(Iterator csvDataBegin, Iterator csvDataEnd);

/// \brief Parse CSV file.
/// \param csvFile Path to CSV file to be parsed.
/// \return Data strucutre containing data parsed out of the CSV file.
CsvGrammarAttribute parseCsvFile(const boost::filesystem::path& csvFile);

/// \brief Insert missing columns into CSV data structure.
/// \param csvData CSV data structure to operate on.
/// \param defaultValue Default value for inserted fields.
/// \return Number of columns per row.
/// \post CSV data structure containing the same number of columns in each row.
CsvGrammarAttribute::value_type::size_type correctCsvColumnsCount(CsvGrammarAttribute &csvData, const std::string &defaultValue = "");


} // namespace common


} // namespace bcl2fastq


#include "common/CsvGrammar.hpp"


#endif //BCL2FASTQ_COMMON_CSV_GRAMMAR_HH


