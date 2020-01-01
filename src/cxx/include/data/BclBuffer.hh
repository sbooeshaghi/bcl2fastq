/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file BclBuffer.hh
 *
 * \brief Declaration of BCL buffer.
 *
 * \author Marek Balint
 */


#ifndef BCL2FASTQ_DATA_BCLBUFFER_HH
#define BCL2FASTQ_DATA_BCLBUFFER_HH


#include <vector>

#include "data/PositionsFile.hh"
#include "data/FilterFile.hh"
#include "data/ControlFile.hh"
#include "data/RawBclBuffer.hh"
#include "layout/Layout.hh"
#include "layout/LaneInfo.hh"
#include "layout/BarcodeTranslationTable.hh"

namespace bcl2fastq {
namespace data {

// Read 5 million clusters from a bcl file at once.
static const size_t ClustersPerReadBuffer = 5000000;

static const size_t ClustersPerTask = 16384;
 
/// \brief BCL buffer.
struct BclBuffer
{
public:

    /// \brief BCL cycles data container type definition.
    typedef std::vector<PerCycleData> BclCycleContainer;

    /// \brief Cluster to sample assignment data container type definition.
    typedef std::vector< std::pair<layout::BarcodeTranslationTable::SampleMetadata, layout::Barcode> > SamplesContainer;

    /// \brief Position data container type definition.
    typedef std::vector<data::PositionsFile::Record> PositionsContainer;

    /// \brief Filter data container type definition.
    typedef std::vector<data::FilterFile::Record> FiltersContainer;

    /// \brief Control data container type definition.
    typedef std::vector<data::ControlFile::Record> ControlsContainer;

    class PatternedPositionsContainer : public PositionsContainer
    {
    public:
        PatternedPositionsContainer() : PositionsContainer(), mut_(), cv_(), isReady_(false) { }

        virtual ~PatternedPositionsContainer() { }

        void setReady() { isReady_ = true; cv_.notify_all(); }

        void waitForReady() { if (!isReady_) { std::unique_lock<std::mutex> lock(mut_); cv_.wait(lock, [this] { return isReady_.load(); }); } }

    private:
        std::mutex mut_;
        std::condition_variable cv_;
        std::atomic<bool> isReady_;
    };

    bool hasCycleData() const {
        for (const PerCycleData & perCycleData : bcls_) {
            if (!(perCycleData.bcls_.empty())) return true;
        }
        return false;
    }

public:

    /// \brief Buffers for cycles data.
    BclCycleContainer bcls_;

    /// \brief Cluster to sample assignments.
    SamplesContainer samples_;

    /// \brief CBCL filter offsets. Used for the offsets for the last cycles which don't include non-pf data.
    std::vector<PerCycleData::BclsContainer::size_type> cbclFilterOffsets_;

    /// \brief Tile meta data.
    layout::LaneInfo::TileInfosContainer::const_iterator tileInfo_;

    /// \brief Buffer for position data.
    std::shared_ptr<PositionsContainer> positions_;

    /// \brief Buffer for filter data.
    FiltersContainer filters_;

    /// \brief Buffer for control data.
    ControlsContainer controls_;

public:

    /// \brief Number of BCL clusters to be processed in single execution iteration.
    /// \todo Refactoring: Put this into CMake configuration.
    static const size_t clustersPerChunk_;
};

typedef std::vector<BclBuffer> BclBufferVec;

/// \brief Swap function for BCL buffer.
/// \param lhs Left-hand-side parameter.
/// \param rhs Right-hand-side parameter.
void swap(BclBuffer &lhs, BclBuffer &rhs);


} // namespace data
} // namespace bcl2fastq

#endif // BCL2FASTQ_DATA_BCLBUFFER_HH

