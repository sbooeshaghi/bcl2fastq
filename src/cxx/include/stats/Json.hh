/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file Json.hh
 *
 * \brief Writes a json file
 *
 * \author Aaron Day
 */


#ifndef BCL2FASTQ_STATS_JSON_HH
#define BCL2FASTQ_STATS_JSON_HH

#include "common/Debug.hh"

#include <boost/noncopyable.hpp>

#include <limits>
#include <string>
#include <vector>
#include <memory>
#include <ostream>

namespace bcl2fastq
{
namespace stats
{

class JsonType
{
public:
    JsonType() { }
    virtual ~JsonType() { }

    virtual void stream(std::ostream& os, int depth = 0) const = 0;

protected:
    void streamIndent(std::ostream& os, int depth) const;
};

class JsonString : public JsonType
{
public:
    JsonString() : JsonType() { }
    JsonString(const std::string& value) : JsonType(), value_(value) { }

    virtual void stream(std::ostream& os, int depth = 0) const;

private:
    const std::string value_;
};

class JsonNumber : public JsonType
{
public:
    JsonNumber() : JsonType(), value_(std::numeric_limits<int64_t>::max()) { }
    JsonNumber(int64_t number) : JsonType(), value_(number) { }

    void add(int64_t number) { BCL2FASTQ_ASSERT_MSG(value_ != std::numeric_limits<int64_t>::max(), "Can't increment a null JsonNumber"); value_ += number; }
    virtual void stream(std::ostream& os, int depth = 0) const;

    int64_t getValue() const { return value_; }

private:
    int64_t value_;
};

class JsonBool : public JsonType
{
public:
    JsonBool(bool value) : JsonType(), value_(value) { }

    virtual void stream(std::ostream& os, int depth = 0) const;

private:
    bool value_;
};

class JsonArray : public JsonType
{
public:
    JsonArray() : JsonType(), values_() { }
    virtual ~JsonArray() { }

    void add(const std::shared_ptr<JsonType>& jsonValue) { values_.push_back(jsonValue); }

    void stream(std::ostream& os, int depth = 0) const;

private:
    std::vector<std::shared_ptr<JsonType>> values_;
};

class JsonObject : public JsonType
{
public:
    JsonObject() : JsonType(), attributes_() { }
    virtual ~JsonObject() { }

    void add(const std::string& key, const std::shared_ptr<JsonType>& jsonValue)
    {
        attributes_.push_back(std::make_pair(key, jsonValue));
    }

    void stream(std::ostream& os, int depth = 0) const;

private:
    std::vector<std::pair<std::string, std::shared_ptr<JsonType>>> attributes_;
};

} // namespace stats
} // namespace bcl2fastq

#endif // BCL2FASTQ_STATS_JSON_HH
