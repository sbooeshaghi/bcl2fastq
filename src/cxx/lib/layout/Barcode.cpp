/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file Barcode.cpp
 *
 * \brief Implementation of barcode.
 *
 * \author Marek Balint
 */


#include <boost/format.hpp>
#include <boost/foreach.hpp>

#include "common/Debug.hh"
#include "common/Exceptions.hh"
#include "layout/Barcode.hh"


namespace bcl2fastq {


namespace layout {


const char Barcode::COMPONENT_SEPARATOR = '+';

Barcode::Component::Component(const std::string &value)
: value_(value)
{
    if (value.empty())
    {
        BOOST_THROW_EXCEPTION(common::LogicError("Barcode component empty"));
    }
}

const char Barcode::DEFAULT_BARCODE[] = "unknown";

void Barcode::reset(size_t numIndexReads)
{
    // Do not clear the vector. We want to avoid reallocations for the strings.
    components_.resize(numIndexReads);

    for (auto& component : components_)
    {
        component.reset();
    }
}

void Barcode::output(std::vector<char>& buffer) const
{
    size_t bufferSize = buffer.size();
    const std::string& firstBarcodeString = components_.front().getString();
    buffer.resize(bufferSize + firstBarcodeString.size());
    std::copy(firstBarcodeString.begin(), firstBarcodeString.end(), buffer.begin() + bufferSize);

    for (size_t i = 1; i < components_.size(); ++i)
    {
        buffer.push_back(COMPONENT_SEPARATOR);

        bufferSize = buffer.size();
        buffer.resize(bufferSize + components_[i].getString().size());
        std::copy(components_[i].getString().begin(), components_[i].getString().end(), buffer.begin() + bufferSize);
    }
}

void Barcode::mask(common::ReadNumber indexNumber)
{
    BCL2FASTQ_ASSERT_MSG(indexNumber <= components_.size(), "Attempting to mask a barcode that doesn't exist.");
    components_.erase(components_.begin()+indexNumber-1);
}

std::ostream& operator<<(std::ostream& os, const Barcode& barcode)
{
    if (barcode.components_.empty())
    {
        os << Barcode::DEFAULT_BARCODE;
    }
    else
    {
        os << barcode.components_[0].getString();
        for (size_t i = 1; i < barcode.components_.size(); ++i)
        {
            os << Barcode::COMPONENT_SEPARATOR;
            os << barcode.components_[i].getString();
        }
    }

    return os;
}

// The conversion operator was abused, and led to terrible memory allocations.
// Use this method sparingly.
std::string Barcode::toString() const
{
    std::string ret;
    if( components_.empty() )
    {
        return DEFAULT_BARCODE;
    }
    BOOST_FOREACH (const Component &component, components_)
    {
        ret.append(component.getString());
        ret += COMPONENT_SEPARATOR;
    }
    ret.resize(ret.size()-1);
    return ret;
}

bool Barcode::operator<(const Barcode &rhs) const
{
    /// \todo Bug: Fix this function to have the same result as string comparison.
    const ComponentsContainer::size_type lhsSize = components_.size();
    const ComponentsContainer::size_type rhsSize = rhs.components_.size();
    if (lhsSize != rhsSize)
    {
        return lhsSize < rhsSize;
    }

    for (ComponentsContainer::size_type i = 0; i < lhsSize; ++i)
    {
        if (components_[i].less(rhs.components_[i]))
        {
            return true;
        }
        if (rhs.components_[i].less(components_[i]))
        {
            return false;
        }
    }
    return false;
}

bool Barcode::operator==(const Barcode &rhs) const
{
    return !(*this != rhs);
}

bool Barcode::operator!=(const Barcode &rhs) const
{
    const ComponentsContainer::size_type lhsSize = components_.size();
    const ComponentsContainer::size_type rhsSize = rhs.components_.size();
    if (lhsSize != rhsSize)
    {
        return true;
    }

    for (ComponentsContainer::size_type i = 0; i < lhsSize; ++i)
    {
        if (components_[i] != rhs.components_[i])
        {
            return true;
        }
    }
    return false;
}

} // namespace layout


} // namespace bcl2fastq


