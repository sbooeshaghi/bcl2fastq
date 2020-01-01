/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file Barcode.hpp
 *
 * \brief Implementation of barcode.
 *
 * \author Marek Balint
 */


#ifndef BCL2FASTQ_LAYOUT_BARCODE_HPP
#define BCL2FASTQ_LAYOUT_BARCODE_HPP


#include <vector>
#include <iterator>

#include <boost/foreach.hpp>

#include "common/Debug.hh"
#include "common/Exceptions.hh"


namespace bcl2fastq {


namespace layout {

template<typename Inserter>
void Barcode::Component::generateMismatches(
    Component::SizeType position,
    Inserter componentInserter
) const
{
    BCL2FASTQ_ASSERT_MSG(position < value_.size(), "Position " << int(position) << " out of range [0," << value_.size() << ")");

    static const char *mismatchesTable[] = {
        "CGTN", // (A)
        "AGTN", // (C)
        "ACGN", // (T)
        "ACGT", // (N)
        "ACTN", // (G)
    };

    size_t mismatchIndex = 0;
    switch (value_[position])
    {
    case 'A':
        mismatchIndex = 0;
        break;
    case 'C':
        mismatchIndex = 1;
        break;
    case 'T':
        mismatchIndex = 2;
        break;
    case 'N':
        mismatchIndex = 3;
        break;
    case 'G':
        mismatchIndex = 4;
        break;
    default:
        BOOST_ASSERT(false);
    };

    // magic number 3: 3 bits per base
    // magic number 0x07: 00000111
    BOOST_FOREACH (const char mismatchedBase, mismatchesTable[mismatchIndex])
    {
        std::string mismatchedValue(value_);
        mismatchedValue.at(position) = mismatchedBase;
        const Component component(mismatchedValue);
        *(*componentInserter) = component;
        ++(*componentInserter);
    }
}


template<typename Iterator>
Barcode::Barcode(Iterator componentsBegin, Iterator componentsEnd)
: components_(componentsBegin, componentsEnd)
{
    if (0) //(componentsBegin == componentsEnd)
    {
        BOOST_THROW_EXCEPTION(common::LogicError("Barcode must have at least one component"));
    }
}

template<typename Inserter>
void Barcode::generateMismatches(
    ComponentsContainer::const_iterator component,
    Component::SizeType position,
    Inserter barcodeInserter
) const
{
    BCL2FASTQ_ASSERT_MSG(component >= components_.begin(), "Invalid component");
    BCL2FASTQ_ASSERT_MSG(component < components_.end(), "Invalid component");

    // magic number 4: There are 5 distinct values ('A', 'C', 'G', 'T', 'N'), therefore there are 4 distinct mismatches
    std::vector<Component> mismatchedComponents;
    mismatchedComponents.reserve(4);
    component->generateMismatches(position, std::back_inserter(mismatchedComponents));

    BOOST_FOREACH (const Component &mismatchedComponent, mismatchedComponents)
    {
        Barcode barcode(*this);
        barcode.components_.at(component-components_.begin()) = mismatchedComponent;
        *(*barcodeInserter) = barcode;
        ++(*barcodeInserter);
    }
}


} // namespace layout


} // namespace bcl2fastq


#endif // BCL2FASTQ_LAYOUT_BARCODE_HPP


