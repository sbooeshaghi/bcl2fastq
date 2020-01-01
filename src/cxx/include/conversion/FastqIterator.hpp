/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file TaskQueue.hpp
 *
 * \brief Declaration of task queue.
 *
 * \author Aaron Day
 */


#ifndef BCL2FASTQ_CONVERSION_FASTQITERATOR_HPP
#define BCL2FASTQ_CONVERSION_FASTQITERATOR_HPP

#include "conversion/FastqIterator.hh"

namespace bcl2fastq {
namespace conversion {

// Silence a g++ warning where the template class function warns the const qualifier is not 
// // needed BUT it is for the class function. I think it looks at this as a non-class function. 
#pragma GCC diagnostic ignored "-Wignored-qualifiers"

template<>
inline typename FastqIteratorTemplate<true, common::NumBasesPerByte::ONE>::reference FastqIteratorTemplate<true, common::NumBasesPerByte::ONE>::dereference() const
{
    return iterator_->bcls_[offset_];
}
template<>
inline typename FastqIteratorTemplate<false, common::NumBasesPerByte::ONE>::reference FastqIteratorTemplate<false, common::NumBasesPerByte::ONE>::dereference() const
{
    return iterator_->bcls_[offset_];
}


template<>
inline typename FastqIteratorTemplate<true, common::NumBasesPerByte::TWO>::reference FastqIteratorTemplate<true, common::NumBasesPerByte::TWO>::dereference() const
{
    auto cycleByte = iterator_->includeNonPf_ ? iterator_->bcls_[offset_] : iterator_->bcls_[offsetCbclFiltered_];
    auto subByteIndex = iterator_->includeNonPf_ ? offsetSubByte_ : offsetCbclFilteredSubByte_;
    return (cycleByte >> 4*subByteIndex) & 0x0f;
}
template<>
inline typename FastqIteratorTemplate<false, common::NumBasesPerByte::TWO>::reference FastqIteratorTemplate<false, common::NumBasesPerByte::TWO>::dereference() const
{
    auto cycleByte = iterator_->includeNonPf_ ? iterator_->bcls_[offset_] : iterator_->bcls_[offsetCbclFiltered_];
    auto subByteIndex = iterator_->includeNonPf_ ? offsetSubByte_ : offsetCbclFilteredSubByte_;
    return (cycleByte >> 4*subByteIndex) & 0x0f;
}

template<>
inline typename FastqIteratorTemplate<true, common::NumBasesPerByte::FOUR>::reference FastqIteratorTemplate<true, common::NumBasesPerByte::FOUR>::dereference() const
{
    auto cycleByte = iterator_->includeNonPf_ ? iterator_->bcls_[offset_] : iterator_->bcls_[offsetCbclFiltered_];
    auto subByteIndex = iterator_->includeNonPf_ ? offsetSubByte_ : offsetCbclFilteredSubByte_;
    return (cycleByte >> 2*subByteIndex) & 0x03;
}
template<>
inline typename FastqIteratorTemplate<false, common::NumBasesPerByte::FOUR>::reference FastqIteratorTemplate<false, common::NumBasesPerByte::FOUR>::dereference() const
{
    auto cycleByte = iterator_->includeNonPf_ ? iterator_->bcls_[offset_] : iterator_->bcls_[offsetCbclFiltered_];
    auto subByteIndex = iterator_->includeNonPf_ ? offsetSubByte_ : offsetCbclFilteredSubByte_;
    return (cycleByte >> 2*subByteIndex) & 0x03;
}

} // namespace conversion
} // namespace bcl2fastq

#endif // BCL2FASTQ_CONVERSION_FASTQITERATOR_HPP

