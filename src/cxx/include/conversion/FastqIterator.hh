/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file FastqIterator.hh
 *
 * \brief Declaration of FASTQ iterator.
 *
 * \author Marek Balint
 */


#ifndef BCL2FASTQ_CONVERSION_FASTQITERATOR_HH
#define BCL2FASTQ_CONVERSION_FASTQITERATOR_HH


#include <vector>
#include <iterator>

#include <boost/iterator/iterator_facade.hpp>
#include <boost/iterator/reverse_iterator.hpp>
#include <boost/mpl/if.hpp>

#include "data/BclBuffer.hh"
#include "data/CbclFile.hh"

namespace bcl2fastq {
namespace conversion {

// VERY IMPORTANT: This class is templated on common::NumBasesPerByte because
// the dereference method must be very fast, since it is called on every base.
// Even a simple switch statement in dereference add 30% to the entire run time.

// NOTE: This class is fundamentally flawed in that it does not allow for the
// use of auto-vecotrization. This is a core piece of the code, and so it should
// be rewritten in such a way that the compiler can vectorize the instructions.

/// \brief FASTQ iterator.
/// \note As Boost::Iterator documentation notes, the superior way to implement
/// this should be by using @c boost::iterator_adaptor. However, that class
/// is flawed, because it allows interoperability of completely unrelated
/// iterators (@c boost::is_convertible is @c true for @c T and @c const @c T,
/// but also for e.g. @c char and @c int).
template<bool IsConstant, common::NumBasesPerByte numBasesPerByte>
class FastqIteratorTemplate : public boost::iterator_facade<
    FastqIteratorTemplate<IsConstant, numBasesPerByte>, // Derived
    typename boost::mpl::if_c<                          // Value
    IsConstant,
    const char,
    char
    >::type,
    boost::random_access_traversal_tag,                 // CategoryOrTraversal
    typename boost::mpl::if_c<                          // Reference
    IsConstant,
    const char,
    char
    >::type
>
{
public:

    /// \brief Constructor.
    /// \param iterator BCL chunk iterator.
    /// \param offset Offset from beginning of the BCL chunk.
    FastqIteratorTemplate(
        const data::BclBuffer& bclBuffer,
        data::PerCycleData::BclsContainer::size_type offset
        )
        : iterator_(bclBuffer.bcls_.begin())
        , offset_(offset / static_cast<uint32_t>(numBasesPerByte))
        , offsetSubByte_(offset % static_cast<uint32_t>(numBasesPerByte))
        , offsetCbclFiltered_(bclBuffer.cbclFilterOffsets_.empty() ? 0 : bclBuffer.cbclFilterOffsets_.at(offset) / static_cast<uint32_t>(numBasesPerByte))
        , offsetCbclFilteredSubByte_(bclBuffer.cbclFilterOffsets_.empty() ? 0 : bclBuffer.cbclFilterOffsets_.at(offset) % static_cast<uint32_t>(numBasesPerByte))
    {
    }

private:

    friend class boost::iterator_core_access;

    /// \brief Dereference core operation.
    /// \return Reference to element being pointed to.
    /// \pre Iterator is dereferencable.
    typename FastqIteratorTemplate<IsConstant, numBasesPerByte>::reference dereference() const;
/*    {
        auto cycleByte = iterator_->includeNonPf_ ? iterator_->bcls_[offset_] : iterator_->bcls_[offsetCbclFiltered_];
        auto subByteIndex = iterator_->includeNonPf_ ? offsetSubByte_ : offsetCbclFilteredSubByte_;
        switch (numBasesPerByte)
        {
        case common::NumBasesPerByte::ONE:
            return cycleByte;
            break;
        case common::NumBasesPerByte::TWO:
            return (cycleByte >> 4*subByteIndex) & 0x0f;
            break;
        case common::NumBasesPerByte::FOUR:
            return (cycleByte >> 2*subByteIndex) & 0x03;
            break;
        }

        BCL2FASTQ_ASSERT_MSG(false, "Invalid number of bases per byte");

        return 0;
    }
*/
    /// \brief Equality comparison core operation.
    /// \param rhs Right-hand-side parameter.
    /// \retval true Left-hand-side and right-hand-side parameters are equal.
    /// \retval false Left-hand-side and right-hand-side parameters are not equal.
    template<bool IsRhsConstant>
    bool equal(const FastqIteratorTemplate<IsRhsConstant, numBasesPerByte> &rhs) const
    {
        return (iterator_ == rhs.iterator_) &&
               (offset_ == rhs.offset_);
    }

    /// \brief Increment core operation.
    /// \pre Iterator is dereferencable.
    /// \post Iterator points to the next element.
    void increment()
    {
        ++iterator_;
    }

    /// \brief Decrement core operation.
    /// \pre Iterator is dereferencable.
    /// \post Iterator points to the next element.
    void decrement()
    {
        --iterator_;
    }

    /// \brief Advance core operation.
    /// \param diff Distance by which the iterator is to be advanced.
    /// \pre Iterator is dereferencable.
    /// \post Iterator has been advanced by specified distance.
    void advance(typename FastqIteratorTemplate<IsConstant, numBasesPerByte>::difference_type diff)
    {
        iterator_ += diff;
    }

    /// \brief Distance core operation.
    /// \param rhs Right-hand-side parameter.
    /// \return Distance between this iterator and right-hand-side iterator.
    template<bool IsRhsConstant>
    typename FastqIteratorTemplate<IsConstant, numBasesPerByte>::difference_type distance_to(const FastqIteratorTemplate<IsRhsConstant, numBasesPerByte> &rhs) const
    {
        return rhs.iterator_ - iterator_;
    }

private:

    /// \brief BCL chunk iterator.
    typename boost::mpl::if_c<IsConstant, data::BclBuffer::BclCycleContainer::const_iterator, data::BclBuffer::BclCycleContainer::iterator>::type iterator_;

    /// \brief Offset from beginning of the BCL chunk.
    data::PerCycleData::BclsContainer::size_type offset_;

    data::PerCycleData::BclsContainer::size_type offsetSubByte_;

    /// \brief Offset from beginning of the CBCL chunk (for cycles which don't include non-pf clusters).
    data::PerCycleData::BclsContainer::size_type offsetCbclFiltered_;

    data::PerCycleData::BclsContainer::size_type offsetCbclFilteredSubByte_;
};

/// \brief FASTQ iterator.
//typedef FastqIteratorTemplate<false, NumBasesPerByte::ONE> FastqIterator;
template<common::NumBasesPerByte numBasesPerByte>
using FastqIterator = FastqIteratorTemplate<false, numBasesPerByte>;

/// \brief FASTQ reverse iterator.
//typedef boost::reverse_iterator<FastqIterator> FastqReverseIterator;
template<common::NumBasesPerByte numBasesPerByte>
using FastqReverseIterator = boost::reverse_iterator<FastqIterator<numBasesPerByte>>;

/// \brief FASTQ constant iterator.
//typedef FastqIteratorTemplate<true, NumBasesPerByte::ONE> FastqConstIterator;
template<common::NumBasesPerByte numBasesPerByte>
using FastqConstIterator = FastqIteratorTemplate<true, numBasesPerByte>;

/// \brief FASTQ constant reverse iterator.
//typedef boost::reverse_iterator<FastqConstIterator> FastqConstReverseIterator;
template<common::NumBasesPerByte numBasesPerByte>
using FastqConstReverseIterator = boost::reverse_iterator<FastqConstIterator<numBasesPerByte>>;


} // namespace conversion
} // namespace bcl2fastq

#include "conversion/FastqIterator.hpp"

#endif // BCL2FASTQ_CONVERSION_FASTQITERATOR_HH
