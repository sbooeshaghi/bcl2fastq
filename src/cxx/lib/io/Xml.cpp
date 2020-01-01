/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file XmlParser.cpp
 *
 * \brief Basic XML parser.
 *
 * \author Roman Petrovski
 * \author Marek Balint
 */


#include <vector>
#include <string>
#include <algorithm>
#include <iterator>

#include <boost/filesystem/fstream.hpp>
#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "io/Xml.hh"
#include "common/Exceptions.hh"

#include "common/Logger.hh"

namespace bcl2fastq {
namespace io {



using namespace boost::property_tree;
static const std::string indexed("<indexed>");


/**
 * @brief Rearranges property tree according to a few basic schema customization rules.
 *
 * The original property tree would produce an xml
 * where each Read, Lane and Tile element has a number of child elements named after
 * the number of the read, lane of tile respectively. For example:
 *
 *  Tree:
 *  Summary.Read.1.Lane.4.Tile.5.value.Yield = 123
 *                                    .YieldQ30 = 12
 *                            .6.value.Yield = 123
 *                                    .YieldQ30 = 12
 *
 *  Xml:
 *  <Summary>
 *    <Read>
 *      <1>
 *        <Lane>
 *          <4>
 *            <Tile>
 *              <5>
 *              ...
 * This is inconvenient and against the rest of the iSAAC metadata xml schemas.
 * Before the tree is written out into the xml, all index nodes are collapsed into
 * xml attribute nodes so that the resulting xml looks like:
 *
 *  Tree:
 *  Summary.<indexed>Read.<number>1.<indexed>Lane.<number>4.<indexed>Tile.<number>5.value.Yield = 123
 *                                                                                      .YieldQ30 = 12
 *                                                                      .<number>6.value.Yield = 123
 *                                                                                      .YieldQ30 = 12
 *
 *  Xml:
 *  <Summary>
 *    <Read index="1">
 *      <Lane number="4">
 *        <Tile number="5">
 *          <Yield>123</Yield>
 *          <YieldQ30>12</YieldQ30>
 *        </Tile>
 *        <Tile number="6">
 *          <Yield>123</Yield>
 *          <YieldQ30>12</YieldQ30>
 *          ...
 */
void unindex(const ptree::value_type &nameElementKv, ptree &resultTree)
{
    if (indexed == nameElementKv.first.substr(0, indexed.length()))
    {
        BOOST_FOREACH(const ptree::value_type &indexElementKv, nameElementKv.second)
        {
            BOOST_ASSERT('<' == indexElementKv.first[0]);
            size_t attrNameEndPos = indexElementKv.first.find('>', 1);
            BOOST_ASSERT(std::string::npos != attrNameEndPos);
            ptree &insertedChildTree = resultTree.add_child(nameElementKv.first.substr(indexed.length()), ptree());

            std::for_each(indexElementKv.second.begin(), indexElementKv.second.end(),
                          boost::bind(unindex, _1, boost::ref(insertedChildTree)));

            if (insertedChildTree.empty())
            {
                // if child has no child elements, use its value as content
                insertedChildTree.put_value(indexElementKv.second.get_value<std::string>());
            }

            insertedChildTree.add("<xmlattr>." + indexElementKv.first.substr(1, attrNameEndPos -1),
                                  indexElementKv.first.substr(attrNameEndPos+1));
        }
    }
    else
    {
        ptree &insertedChildTree = resultTree.add_child(nameElementKv.first, ptree());
        std::for_each(nameElementKv.second.begin(), nameElementKv.second.end(),
                      boost::bind(unindex, _1, boost::ref(insertedChildTree)));
        if (insertedChildTree.empty())
        {
            //unindexing did not find any indexed nodes. just add it as is
            insertedChildTree.put_value(nameElementKv.second.get_value<std::string>());
        }
        else
        {
            //unindexing inserted the modified tree.
        }
    }
}

bool longestPathFirst(const std::string &left, const std::string &right)
{
    return std::count(left.begin(), left.end(), '.') > std::count(right.begin(), right.end(), '.');
}

void index(std::vector<std::string > indexNodes, ptree &tree)
{
    // indexing of longer paths should happen first for the shorter paths to stay valid until they get indexed
    std::sort(indexNodes.begin(), indexNodes.end(), longestPathFirst);

    BOOST_FOREACH(const std::string &idxAttrPath, indexNodes)
    {
        size_t attrNameDotPos(idxAttrPath.rfind('.'));
        BOOST_ASSERT(std::string::npos != attrNameDotPos && "Could not parse index attribute name");
        const std::string indexedNodeKey(idxAttrPath.substr(0, attrNameDotPos));
        const std::string indexAttributeName(idxAttrPath.substr(attrNameDotPos+1));

        size_t indexedNodeNameDotPos(indexedNodeKey.rfind('.'));
        BOOST_ASSERT(std::string::npos != indexedNodeNameDotPos && "Could not parse indexed node name");
        const std::string indexedNodeParentKey(indexedNodeKey.substr(0, indexedNodeNameDotPos));
        const std::string indexedNodeName(indexedNodeKey.substr(indexedNodeNameDotPos+1));

        ptree &indexedNodeParent(tree.get_child(indexedNodeParentKey));
        BOOST_FOREACH(ptree::value_type &potentiallyIndexedNode, indexedNodeParent)
        {
            if (indexedNodeName == potentiallyIndexedNode.first)
            {
                std::string indexedNodeKey("<indexed>" + indexedNodeName + ".<" + indexAttributeName + ">"
                                           + potentiallyIndexedNode.second.get<std::string>("<xmlattr>." + indexAttributeName));
                BOOST_ASSERT(potentiallyIndexedNode.second.get_child("<xmlattr>").erase(indexAttributeName));
                indexedNodeParent.add_child(indexedNodeKey, potentiallyIndexedNode.second);
            }
        }
        indexedNodeParent.erase(indexedNodeName);
    }
}

std::ostream &serializeAsXml(std::ostream &os, const boost::property_tree::ptree &tree)
{
    boost::property_tree::ptree treeWithIndexAttributes;
    if (!tree.empty())
    {
        unindex(*tree.begin(), treeWithIndexAttributes);
#ifndef WIN32
        boost::property_tree::write_xml(os, treeWithIndexAttributes, boost::property_tree::xml_writer_make_settings(' ', 2));
#else
        boost::property_tree::write_xml(os, treeWithIndexAttributes, boost::property_tree::xml_writer_make_settings<std::string>(' ', 2));
#endif
    }
    else
    {
#ifndef WIN32
        boost::property_tree::write_xml(os, tree, boost::property_tree::xml_writer_make_settings(' ', 2));
#else
        boost::property_tree::write_xml(os, tree, boost::property_tree::xml_writer_make_settings<std::string>(' ', 2));
#endif
    }
    return os;
}

std::istream &parseAsXml(std::istream &is, boost::property_tree::ptree &tree)
{
    boost::property_tree::read_xml<boost::property_tree::ptree>(is, tree);
    return is;
}

boost::property_tree::ptree parseXmlFile(const boost::filesystem::path& xmlFile)
{
    boost::filesystem::ifstream xmlStream(xmlFile);
    xmlStream.unsetf(std::ios::skipws);
    int errnum = errno;
    if (!xmlStream)
    {
        BOOST_THROW_EXCEPTION(common::IoError(errnum, (boost::format("Unable to open '%s' file for reading") % xmlFile.string()).str()));
    }

    typedef std::vector<char> XmlBuffer;
    XmlBuffer xmlBuffer;
    errnum = 0;
    std::copy(
        std::istream_iterator<XmlBuffer::value_type>(xmlStream),
        std::istream_iterator<XmlBuffer::value_type>(),
        std::back_inserter(xmlBuffer)
    );

    return parseXmlData(xmlBuffer.begin(), xmlBuffer.end());
}


} // namespace io
} // namespace bcl2fastq


