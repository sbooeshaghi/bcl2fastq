/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file BarcodeStats.hh
 *
 * \brief Declaration of Barcode Stats.
 *
 * \author Mauricio Varea
 */


#ifndef BCL2FASTQ_STATS_BARCODESTATS_HH
#define BCL2FASTQ_STATS_BARCODESTATS_HH


#include "common/Types.hh"
#include "layout/BarcodeTranslationTable.hh"

#include <boost/noncopyable.hpp>
#include <numeric>

namespace bcl2fastq
{
namespace stats
{


/// \brief Demultiplexing statistics for single barcode.
class BarcodeStats
{
public:
    // layout::BarcodeTranslationTable::SampleMetadata index;
    typedef std::pair< layout::BarcodeTranslationTable::SampleMetadata,
                       layout::LaneInfo::TileInfosContainer::size_type > Index;

    Index index;

    common::TotalClustersCount getBarcodeCount() const {
        return static_cast<common::TotalClustersCount>(
            std::accumulate(
                barcodeMismatchCounts_.begin(),
                barcodeMismatchCounts_.end(),
                static_cast<common::TotalClustersCount>(0)));
    }

    // This method is used only if there was no demultiplexing.
    void setBarcodeCount(common::TotalClustersCount barcodeCount)
    {
        BOOST_ASSERT(getBarcodeCount() == 0);

        barcodeMismatchCounts_[0] = barcodeCount;
    }

    void incrementBarcodeCount(common::TotalClustersCount numMismatches)
    {
        BOOST_ASSERT(numMismatches <= MAX_MISMATCHES);
        ++barcodeMismatchCounts_[numMismatches];
    }

    static const size_t MAX_MISMATCHES = 5;

    // The index into the array is the number of mismatches
    const std::array<common::TotalClustersCount, MAX_MISMATCHES+1>& getBarcodeMismatchCounts() const { return barcodeMismatchCounts_; }

private:
    std::array<common::TotalClustersCount, MAX_MISMATCHES+1> barcodeMismatchCounts_;

    static Index defaultIndex() { return std::make_pair( layout::BarcodeTranslationTable::SampleMetadata(), 0 ); }

public:
    BarcodeStats( const layout::BarcodeTranslationTable::SampleMetadata &sampleMetadata,
                  layout::LaneInfo::TileInfosContainer::size_type tileIndex )
    : index( std::make_pair(sampleMetadata,tileIndex) ), barcodeMismatchCounts_() {}
    BarcodeStats(const Index &idx)
    : index(idx), barcodeMismatchCounts_() {}
    BarcodeStats()
    : index( defaultIndex() ), barcodeMismatchCounts_() {}

    ~BarcodeStats() { /*BCL2FASTQ_LOG(common::LogLevel::DEBUG) << "Stats going away" << std::endl;*/ }
    BarcodeStats &operator+=(const BarcodeStats &rhs)
    {
        for (size_t i = 0; i < MAX_MISMATCHES+1; ++i)
        {
            barcodeMismatchCounts_[i] += rhs.barcodeMismatchCounts_[i];
        }
        return *this;
    }

    bool operator<( const BarcodeStats &rhs ) const
    {
        return index < rhs.index;
    }

    void reset()
    {
        barcodeMismatchCounts_.fill(0);
    }
};


} // namespace stats
} // namespace bcl2fastq


#endif // BCL2FASTQ_STATS_BARCODESTATS_HH


