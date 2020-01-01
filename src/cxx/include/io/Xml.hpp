/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file XmlParser.hpp
 *
 * \brief Basic XML parser.
 *
 * \author Marek Balint
 */


#ifndef BCL2FASTQ_IO_XML_HPP
#define BCL2FASTQ_IO_XML_HPP


#include <sstream>
#include <string>

#include <boost/property_tree/xml_parser.hpp>


namespace bcl2fastq {
namespace io {


template<typename Iterator>
boost::property_tree::ptree parseXmlData(Iterator xmlDataBegin, Iterator xmlDataEnd)
{
    boost::property_tree::ptree ret;
    std::istringstream istream(std::string(xmlDataBegin, xmlDataEnd));
    boost::property_tree::read_xml<boost::property_tree::ptree>(static_cast<std::basic_istream<char> &>(istream), ret);
    return ret;
}


} // namespace io
} // namespace bcl2fastq


#endif //BCL2FASTQ_IO_XML_HPP


