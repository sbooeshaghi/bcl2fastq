/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file SampleIndex.cpp
 *
 * \brief Implementation of sample index.
 *
 * \author Marek Balint
 * \author Mauricio Varea
 */


#include <algorithm>

#include "common/Debug.hh"
#include "conversion/SampleIndex.hh"


namespace bcl2fastq {


namespace conversion {


SampleIndex::SampleIndex(
    layout::LaneInfo::SampleInfosContainer::size_type samplesCount,
    SampleIndex::FastqOffsetsContainer::size_type maxFastqsCount
)
: state_(State::INITIALIZING)
, offsets_(maxFastqsCount)
, samplesMetaData_(samplesCount, SampleIndex::SampleMetaData())
, fastqsCount_()
{
}

void SampleIndex::reset(const data::BclBufferVec& bclBuffers)
{
    reset(getTotalClusters(bclBuffers));

    for (const auto& inputBuffer : bclBuffers)
    {
         // find out how many reads in each sample
         size_t bufferSize = inputBuffer.samples_.size();
         for (data::BclBuffer::BclCycleContainer::size_type offset = 0; offset < bufferSize; ++offset)
         {
             incrementFastqCount(inputBuffer.samples_[offset].first.sampleIndex_);
         }
    }

    reserve();

    // initialize offsets
    size_t bufferIndex = 0;
    for (const auto& inputBuffer : bclBuffers)
    {
        size_t bufferSize = inputBuffer.samples_.size();
        for (data::BclBuffer::BclCycleContainer::size_type offset = 0; offset < bufferSize; ++offset)
        {
            addOffset(inputBuffer.samples_[offset].first.sampleIndex_,
                SampleIndex::FastqOffsetsContainer::value_type(bufferIndex, offset));
        }

        ++bufferIndex;
    }

    finalize();
}

size_t SampleIndex::getTotalClusters(const data::BclBufferVec& bclBuffers) const
{
    size_t dataSize = 0;

    for (const auto& inputBuffer : bclBuffers)
    {
        dataSize += inputBuffer.samples_.size();
    }

    return dataSize;
}


void SampleIndex::reset(SampleIndex::FastqOffsetsContainer::size_type fastqsCount)
{
    fastqsCount_ = fastqsCount;
    if( fastqsCount > offsets_.size() )
    {
        BCL2FASTQ_LOG(common::LogLevel::DEBUG) << "Resizing offsets table to cope with " << fastqsCount << " clusters/tile" << std::endl;
        offsets_.resize(fastqsCount + defaultFastqStep_);
    }
    SampleMetaData smd;
    smd.begin = offsets_.begin();
    smd.end = offsets_.begin();
    smd.next = offsets_.begin();
    std::fill(samplesMetaData_.begin(), samplesMetaData_.end(), smd);

    state_ = State::CALCULATING_SIZES;
}

void SampleIndex::incrementFastqCount(
    layout::LaneInfo::SampleInfosContainer::size_type sample,
    SampleIndex::FastqOffsetsContainer::size_type count
)
{
    BCL2FASTQ_ASSERT_MSG( sample < samplesMetaData_.size(),
                          "Invalid sample index: " << sample << "/" << (samplesMetaData_.size()-1));
    BCL2FASTQ_ASSERT_MSG( offsets_.size() > static_cast<size_t>(samplesMetaData_[sample].end - samplesMetaData_[sample].begin) + count,
                          "Sample index #" << sample << " tries to increment FASTQ count past the end of the offsets table"
                          " (size: " << offsets_.size() << ")" );

    std::advance(samplesMetaData_[sample].end, count);
}

void SampleIndex::reserveHelper(
    std::vector<SampleIndex::SampleMetaData>::size_type &offset,
    SampleIndex::SampleMetaData &sampleMetaData
)
{
    sampleMetaData.begin += offset;
    sampleMetaData.end += offset;
    sampleMetaData.next = sampleMetaData.begin;

    offset += sampleMetaData.end - sampleMetaData.begin;
}

void SampleIndex::reserve()
{
    BCL2FASTQ_ASSERT_MSG(
        state_ == State::CALCULATING_SIZES,
        "Sample index must be in state CALCULATING_SIZES: " << state_
    );

    std::vector<SampleIndex::SampleMetaData>::size_type offset = 0;
    std::for_each(
        samplesMetaData_.begin(),
        samplesMetaData_.end(),
        boost::bind(
            reserveHelper,
            boost::ref(offset),
            _1
        )
    );

    BCL2FASTQ_ASSERT_MSG(
        samplesMetaData_.front().begin == offsets_.begin(),
        "The first offset of the first sample must be at the first position of offsets array"
    );
    BCL2FASTQ_ASSERT_MSG(
        samplesMetaData_.back().end == (offsets_.begin() + fastqsCount_),
        "The last offset of the last sample must be at the position of last FASTQ offset in offsets array"
    );

    state_ = State::BUILDING_INDEX;
}

void SampleIndex::addOffset(
    layout::LaneInfo::SampleInfosContainer::size_type sample,
    const FastqOffsetsContainer::value_type& offset
)
{
    BCL2FASTQ_ASSERT_MSG(
        sample < samplesMetaData_.size(),
        "Invalid sample index: " << sample << "/" << (samplesMetaData_.size()-1)
    );
    BCL2FASTQ_ASSERT_MSG(
        state_ == State::BUILDING_INDEX,
        "Sample index must be in state BUILDING_INDEX: " << state_
    );
    BCL2FASTQ_ASSERT_MSG(
        samplesMetaData_[sample].next < samplesMetaData_[sample].end,
        "Too many offsets for sample: " << sample << ":" << (samplesMetaData_[sample].end - samplesMetaData_[sample].begin)
    );

    *(samplesMetaData_[sample].next++) = offset;
}

bool SampleIndex::finalizePredicate(
    SampleIndex::FastqOffsetsContainer::const_iterator offsetsBegin,
    SampleIndex::FastqOffsetsContainer::const_iterator offsetsEnd,
    const SampleIndex::SampleMetaData &sampleMetaData
)
{
    // mind reversed logic of this predicate
    return !(
        sampleMetaData.begin >= offsetsBegin
        &&
        sampleMetaData.begin <= offsetsEnd
        &&
        sampleMetaData.end >= offsetsBegin
        &&
        sampleMetaData.end <= offsetsEnd
        &&
        sampleMetaData.next >= offsetsBegin
        &&
        sampleMetaData.next <= offsetsEnd
        &&
        sampleMetaData.begin <= sampleMetaData.end
        &&
        sampleMetaData.next == sampleMetaData.end
    );
}

void SampleIndex::finalize()
{
    BCL2FASTQ_ASSERT_MSG(
        state_ == State::BUILDING_INDEX,
        "Sample index must be in state BUILDING_INDEX: " << state_
    );

    {
        std::vector<SampleMetaData>::const_iterator invalidMetaData = find_if(
            samplesMetaData_.begin(),
            samplesMetaData_.end(),
            boost::bind(
                SampleIndex::finalizePredicate,
                offsets_.begin(),
                offsets_.end(),
                _1
            )
        );

#if 0
        BCL2FASTQ_LOG(common::LogLevel::TRACE) << "SAMPLE INDEX: " << std::endl;
        BOOST_FOREACH( const SampleIndex::SampleMetaData &sm, std::make_pair( samplesMetaData_.begin(), samplesMetaData_.end() ))
        {
            BclBuffer::BclsContainer::value_type::size_type &begin = *(sm.begin);
            BclBuffer::BclsContainer::value_type::size_type &end   = *(sm.end);
            BclBuffer::BclsContainer::value_type::size_type &next  = *(sm.next);

            BCL2FASTQ_LOG(common::LogLevel::TRACE) << "\t{" << begin
                                                   << ", " << end
                                                   << ", " << next << "}" << std::endl;
        }
#endif

        BCL2FASTQ_ASSERT_MSG(
            invalidMetaData == samplesMetaData_.end(),
            "Invalid index: " << (invalidMetaData - samplesMetaData_.begin())
        );
    }

    state_ = State::INDEX_BUILT;
}

SampleIndex::FastqOffsetsContainer::const_iterator SampleIndex::offsetsBegin(layout::LaneInfo::SampleInfosContainer::size_type sample) const
{
    BCL2FASTQ_ASSERT_MSG(
        sample < samplesMetaData_.size(),
        "Invalid sample index: " << sample << "/" << (samplesMetaData_.size()-1)
    );
    BCL2FASTQ_ASSERT_MSG(
        state_ == State::INDEX_BUILT,
        "Sample index must be in state INDEX_BUILT: " << state_
    );

    return samplesMetaData_[sample].begin;
}

SampleIndex::FastqOffsetsContainer::const_iterator SampleIndex::offsetsEnd(layout::LaneInfo::SampleInfosContainer::size_type sample) const
{
    BCL2FASTQ_ASSERT_MSG(
        sample < samplesMetaData_.size(),
        "Invalid sample ID: " << sample << "/" << (samplesMetaData_.size()-1)
    );
    BCL2FASTQ_ASSERT_MSG(
        state_ == State::INDEX_BUILT,
        "Sample index must be in state INDEX_BUILT: " << state_
    );

    return samplesMetaData_[sample].end;
}

SampleIndex::FastqOffsetsContainer::iterator SampleIndex::offsetsBegin(layout::LaneInfo::SampleInfosContainer::size_type sample)
{
    BCL2FASTQ_ASSERT_MSG(
        sample < samplesMetaData_.size(),
        "Invalid sample index: " << sample << "/" << (samplesMetaData_.size()-1)
    );
    BCL2FASTQ_ASSERT_MSG(
        state_ == State::INDEX_BUILT,
        "Sample index must be in state INDEX_BUILT: " << state_
    );

    return samplesMetaData_[sample].begin;
}

SampleIndex::FastqOffsetsContainer::iterator SampleIndex::offsetsEnd(layout::LaneInfo::SampleInfosContainer::size_type sample)
{
    BCL2FASTQ_ASSERT_MSG(
        sample < samplesMetaData_.size(),
        "Invalid sample ID: " << sample << "/" << (samplesMetaData_.size()-1)
    );
    BCL2FASTQ_ASSERT_MSG(
        state_ == State::INDEX_BUILT,
        "Sample index must be in state INDEX_BUILT: " << state_
    );

    return samplesMetaData_[sample].end;
}


} // namespace conversion


} // namespace bcl2fastq


