/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file CsvGrammar.hpp
 *
 * \brief Basic CSV grammar.
 *
 * \author Roman Petrovski
 * \author Marek Balint
 */


#ifndef BCL2FASTQ_COMMON_CSVGRAMMAR_HPP
#define BCL2FASTQ_COMMON_CSVGRAMMAR_HPP


#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>

#include "common/Exceptions.hh"


namespace bcl2fastq {


namespace common {


template <typename Iterator>
CsvGrammar<Iterator>::CsvGrammar()
: CsvGrammar::base_type(start_)
{
    using boost::spirit::omit;
    using boost::spirit::qi::_val;
    using boost::spirit::qi::_1;
    using boost::spirit::ascii::char_;
    using boost::phoenix::push_back;

    // CSV ABNF mainly taken from http://www.ietf.org/rfc/rfc4180.txt
    // with some modification to make sample sheet users life easier
    cr_ = char_('\x0d');
    lf_ = char_('\x0a');
    crlf_ = +(cr_ | lf_); // unlike the real csv, sample sheets often get mixed line endings
    comma_ = char_(',');
    comma_.name(",");
    dquote_ = char_('"');

    //textData_ = char_('\x20', '\x21') | char_('\x23', '\x2b') | char_('\x2d', '\x7e');
    textData_ = char_('\x20', '\x7e') - (char_(',') | char_('"'));

    /// \todo Refactoring: Probably omit[dquote_] should be used.
    escaped_ = dquote_
        >> *(textData_[push_back(_val, _1)] |
            comma_[push_back(_val, _1)] |
            cr_[push_back(_val, _1)] |
            lf_[push_back(_val, _1)] |
            (dquote_ >> dquote_)[push_back(_val, '"')])
        >> dquote_
    ;

    nonEscaped_ = *textData_;

    field_ = escaped_ | nonEscaped_;

    commentLine_ = char_('#') >> *(char_('\x00', '\x7e') - crlf_);

    dataLine_ = field_ % comma_;

    line_ = *(commentLine_ >> +crlf_) >> dataLine_ >> *(+crlf_ >> commentLine_);

    start_ = line_ % crlf_ >> *crlf_;
}


template<typename Iterator>
CsvGrammarAttribute parseCsvData(Iterator csvDataBegin, Iterator csvDataEnd)
{
    common::CsvGrammar<Iterator> parser;

    common::CsvGrammarAttribute csvData;
    if (!boost::spirit::qi::parse(csvDataBegin, csvDataEnd, parser, csvData) || csvDataEnd != csvDataBegin)
    {
        int errnum = errno;
        const std::string message = "Could not parse the CSV stream text:\n"
             + std::string(csvDataBegin, csvDataEnd)
        ;
        BOOST_THROW_EXCEPTION(common::InputDataError(errnum, message));
    }

    correctCsvColumnsCount(csvData, "");
    return csvData;
}


} // namespace common


} // namespace bcl2fastq


#endif //BCL2FASTQ_COMMON_CSVGRAMMAR_HPP


