/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file SampleIndex.hh
 *
 * \brief Declaration of sample index.
 *
 * \author Marek Balint
 * \author Mauricio Varea
 */


#ifndef BCL2FASTQ_CONVERSION_SAMPLEINDEX_HH
#define BCL2FASTQ_CONVERSION_SAMPLEINDEX_HH


#include <boost/noncopyable.hpp>

#include "layout/LaneInfo.hh"
#include "data/BclBuffer.hh"


namespace bcl2fastq {
namespace conversion {


/// \brief Index of FASTQ offsets and their assignments to samples.
class SampleIndex : private boost::noncopyable
{
public:
    /// \brief FASTQ offsets container type definition.
    typedef std::vector<std::pair<data::BclBufferVec::size_type, common::RawDataBuffer::size_type>> FastqOffsetsContainer;

    static const FastqOffsetsContainer::size_type defaultFastqStep_ = 1;

public:

    /// \brief Constructor.
    /// \param samplesCount Total number of samples.
    /// \param maxFastqsCount Maximum total number of indexed FASTQ sequences.
    SampleIndex(
        layout::LaneInfo::SampleInfosContainer::size_type samplesCount,
        SampleIndex::FastqOffsetsContainer::size_type maxFastqsCount
    );

    /// \brief Reset to initial state.
    /// \param bclBuffers buffers to index
    void reset(const data::BclBufferVec& bclBuffers);

private:

    /// \brief Get the total number of clusters.
    /// \param bclBuffers bcl buffers.
    /// \return Iterator to end of offsets.
    size_t getTotalClusters(const data::BclBufferVec& bclBuffers) const;

    /// \brief Reset to initial state.
    /// \param fastqsCount Actual total number of indexed FASTQ sequences.
    /// \pre Value passed in @c fastqsCount parameter must not be larger than
    /// value passed to constructor's @c maxFastqsCount paramteter.
    void reset(SampleIndex::FastqOffsetsContainer::size_type fastqsCount);

    /// \brief Increment count of FASTQ offsets for given sample.
    /// \param sample Sample index.
    /// \param count Number of FASTQ offsets to be counted in.
    /// \pre Neither @c reserve nor @c finalize has been called since
    /// last call to @c reset, which must have been called at least once.
    void incrementFastqCount(
        layout::LaneInfo::SampleInfosContainer::size_type sample,
        SampleIndex::FastqOffsetsContainer::size_type count = defaultFastqStep_
    );

    /// \brief Create data structures for specified amounts of FASTQ offsets in individual samples.
    /// \pre Neither @c reserve nor @c finalize has been called since
    /// last call to @c reset, which must have been called at least once.
    void reserve();

    /// \brief Add FASTQ offset for given sample.
    /// \param sample Sample index.
    /// \param offset FASTQ offset to be added.
    /// \pre Neither @c reset nor @c finalize has been called since
    /// last call to @c reserve, which must have been called at least once.
    void addOffset(
        layout::LaneInfo::SampleInfosContainer::size_type sample,
        const FastqOffsetsContainer::value_type& offset
    );

    /// \brief Finalize index.
    /// \pre Neither @c finalize nor @c reset has been called since
    /// last call to @c reserve, which must have been called at least once.
    void finalize();

public:

    /// \brief Get beginning of offsets for given sample.
    /// \param sample Sample index.
    /// \return Iterator to beginning of offsets.
    /// \pre Neither @c reset nor @c reserve has been called since
    /// last call to @c finalize, which must have been called at least once.
    SampleIndex::FastqOffsetsContainer::const_iterator offsetsBegin(layout::LaneInfo::SampleInfosContainer::size_type sample) const;

    /// \brief Get end of offsets for given sample.
    /// \param sample Sample index.
    /// \return Iterator to end of offsets.
    /// \pre Neither @c reset nor @c reserve has been called since
    /// last call to @c finalize, which must have been called at least once.
    SampleIndex::FastqOffsetsContainer::const_iterator offsetsEnd(layout::LaneInfo::SampleInfosContainer::size_type sample) const;

    /// \brief Get beginning of offsets for given sample.
    /// \param sample Sample index.
    /// \return Iterator to beginning of offsets.
    /// \pre Neither @c reset nor @c reserve has been called since
    /// last call to @c finalize, which must have been called at least once.
    SampleIndex::FastqOffsetsContainer::iterator offsetsBegin(layout::LaneInfo::SampleInfosContainer::size_type sample);

    /// \brief Get end of offsets for given sample.
    /// \param sample Sample index.
    /// \return Iterator to end of offsets.
    /// \pre Neither @c reset nor @c reserve has been called since
    /// last call to @c finalize, which must have been called at least once.
    SampleIndex::FastqOffsetsContainer::iterator offsetsEnd(layout::LaneInfo::SampleInfosContainer::size_type sample);

private:

    /// \brief Internal state.
    struct State
    {
        enum value_type
        {
            INITIALIZING,      ///< Object has not been fully initialized, yet.
            CALCULATING_SIZES, ///< Write-only state: accumulating sample sizes.
            BUILDING_INDEX,    ///< Write-only state: acquiring declared numbers of FASTQ offsets.
            INDEX_BUILT        ///< Read-only state: index is ready for use.
        };
    };

    /// \brief Sample meta-data.
    struct SampleMetaData
    {
        /// \brief Beginning of offsets list.
        FastqOffsetsContainer::iterator begin;

        /// \brief End of offsets list.
        FastqOffsetsContainer::iterator end;

        /// \brief Next offset record to be filled in.
        FastqOffsetsContainer::iterator next;
    };

private:

    /// \brief Helper for reserve method.
    /// \param offset Reference to offset variable, being used inside the algorithm.
    /// \param sampleMetaData Sample meta data structure.
    static void reserveHelper(
        std::vector<SampleIndex::SampleMetaData>::size_type &offset,
        SampleIndex::SampleMetaData &sampleMetaData
    );

    /// \brief Helper predicate for finalize method.
    /// \param offsetsBegin Beginning of FASTQ offset array.
    /// \param offsetsEnd End of FASTQ offset array.
    /// \param sampleMetaData Sample meta data structure.
    /// \retval true Sample meta data structure *IS* *NOT* invalid.
    /// \retval false Sample meta data structure *IS* valid.
    static bool finalizePredicate(
        SampleIndex::FastqOffsetsContainer::const_iterator offsetsBegin,
        SampleIndex::FastqOffsetsContainer::const_iterator offsetsEnd,
        const SampleIndex::SampleMetaData &sampleMetaData
    );

private:
public:
    /// \brief Internal state.
    State::value_type state_;

    /// \brief FASTQ offset array.
    FastqOffsetsContainer offsets_;

    /// \brief Samples metadata.
    std::vector<SampleMetaData> samplesMetaData_;

    /// \brief Number of FASTQ sequences (actual number as opposed to maximum number).
    FastqOffsetsContainer::size_type fastqsCount_;
};

struct DemuxBuffer
{
    DemuxBuffer()
        : bclBuffers_(), sampleIndex_(), groupNumber_(0)
    {
    }

    DemuxBuffer(layout::LaneInfo::SampleInfosContainer::size_type samplesCount,
                SampleIndex::FastqOffsetsContainer::size_type maxFastqsCount,
                size_t numCycles)
        : bclBuffers_(),
          sampleIndex_(std::make_shared<SampleIndex>(samplesCount, maxFastqsCount)),
          groupNumber_(0)
    {
    }

    void init(layout::LaneInfo::SampleInfosContainer::size_type samplesCount,
              SampleIndex::FastqOffsetsContainer::size_type maxFastqsCount)
    {
        sampleIndex_ = std::make_shared<SampleIndex>(samplesCount, maxFastqsCount);
    }

    void calculateIndex()
    {
        sampleIndex_->reset(bclBuffers_);
    }

    bool hasCycleData() const 
    {
        for (const data::BclBuffer &bclBuffer : bclBuffers_) {
            if (bclBuffer.hasCycleData()) 
                 return true;
        }
        return false;
    }

    uint32_t getGroupNumber() const { return groupNumber_; }
    void setGroupNumber(uint32_t groupNumber) { groupNumber_ = groupNumber; }

    data::BclBufferVec bclBuffers_;
    std::shared_ptr<SampleIndex> sampleIndex_;
    uint32_t groupNumber_;
};


} // namespace conversion
} // namespace bcl2fastq


#endif // BCL2FASTQ_CONVERSION_SAMPLEINDEX_HH


