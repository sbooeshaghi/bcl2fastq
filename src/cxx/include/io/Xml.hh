/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file Xml.hh
 *
 * \brief Basic XML parser.
 *
 * \author Marek Balint
 * \author Roman Petrovski
 */


#ifndef BCL2FASTQ_IO_XML_HH
#define BCL2FASTQ_IO_XML_HH


#include <boost/filesystem/path.hpp>
#include <boost/property_tree/ptree.hpp>


namespace bcl2fastq {
namespace io {

void index(std::vector<std::string> indexAttributes, boost::property_tree::ptree &tree);

std::ostream &serializeAsXml(std::ostream &os, const boost::property_tree::ptree &tree);

std::istream &parseAsXml(std::istream &is, boost::property_tree::ptree &tree);


/// \brief Parse XML data.
/// \param xmlDataBegin Begin iterator to XML data to be parsed.
/// \param xmlDataEnd End iterator to XML data to be parsed.
/// \return Data strucutre containing data parsed out of the XML data.
template<typename Iterator>
boost::property_tree::ptree parseXmlData(Iterator xmlDataBegin, Iterator xmlDataEnd);

/// \brief Parse XML file.
/// \param xmlFile Path to XML file to be parsed.
/// \return Data strucutre containing data parsed out of the XML file.
boost::property_tree::ptree parseXmlFile(const boost::filesystem::path& xmlFile);


} // namespace io
} // namespace bcl2fastq


#include "io/Xml.hpp"


#endif //BCL2FASTQ_IO_XML_HH


