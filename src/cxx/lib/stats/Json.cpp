/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file JsonStats.cpp
 *
 * \brief Json creation
 *
 * \author Aaron Day
 */

#include "stats/Json.hh"

namespace bcl2fastq
{
namespace stats
{

void JsonType::streamIndent(std::ostream& os, int depth) const
{
    for (int i = 0; i < depth; ++i)
    {
        os << " ";
    }
}

void JsonString::stream(std::ostream& os, int depth /*= 0*/) const
{
    os << "\"" << value_ << "\"";
}

void JsonNumber::stream(std::ostream& os, int depth /*= 0*/) const
{
    if (value_ == std::numeric_limits<int64_t>::max())
    {
        os << "null";
    }
    else 
    {
        os << value_;
    }
}

void JsonBool::stream(std::ostream& os, int depth /*= 0*/) const
{
    os << (value_ ? "true" : "false");
}


void JsonArray::stream(std::ostream& os, int depth /*= 0*/) const
{
    os << "[" << std::endl;

    bool first = true;
    for (auto& jsonValue : values_)
    {
        if (first)
        {
            first = false;
        }
        else
        {
            os << "," << std::endl;
        }

        streamIndent(os, depth+1);
        jsonValue->stream(os, depth+1);
    }
    os << std::endl;
    streamIndent(os, depth);
    os << "]";
}

void JsonObject::stream(std::ostream& os, int depth /*= 0*/) const
{
    os << "{" << std::endl;
    int nextDepth = depth + 1;

    bool first = true;
    for (auto& jsonValue : attributes_)
    {
        if (first)
        {
            first = false;
        }
        else
        {
           os << "," << std::endl;
        }

        streamIndent(os, nextDepth);
        os << "\"" << jsonValue.first << "\" : ";
        jsonValue.second->stream(os, nextDepth);
    }

    os << std::endl;
    streamIndent(os, depth);
    os << "}";
}

} //namespace stats
} //namespace bcl2fastq

