/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file TileStats.hpp
 *
 * \brief Declaration of Tile Stats.
 *
 * \author Mauricio Varea
 */


#ifndef BCL2FASTQ_STATS_TILESTATS_HPP
#define BCL2FASTQ_STATS_TILESTATS_HPP


#include "common/Types.hh"
#include "layout/BarcodeTranslationTable.hh"
#include "layout/Layout.hh"



namespace bcl2fastq
{
namespace stats
{


/// \brief Across-all-reads statistics for a single tile.
struct AllReadsStats
{
    common::TotalClustersCount clustersCount;

    typedef std::vector<size_t> TrimLengthCount;
    typedef std::vector<TrimLengthCount> TrimLengthCountForReads;

    TrimLengthCountForReads trimLengthCountForReads;

    AllReadsStats()
    : clustersCount(0)
    , trimLengthCountForReads()
    {}

    /*AllReadsStats(size_t readLength)
    : clustersCount(0)
    , trimLengthCount(readLength+1)
    {}*/

    void init(size_t readLength,
              size_t numNonIndexReads)
    {
        BOOST_ASSERT(readLength != 0);

        trimLengthCountForReads.resize(numNonIndexReads);

        BOOST_FOREACH(TrimLengthCount& trimLengthCount, trimLengthCountForReads)
        {
            trimLengthCount.resize(readLength+1);
        }
    }

    void RecordAdapterStats(size_t trimmedBases,
                            size_t readNumber)
    {
        BOOST_ASSERT(readNumber <= trimLengthCountForReads.size());
        BOOST_ASSERT(trimmedBases < trimLengthCountForReads[readNumber-1].size());
        
        trimLengthCountForReads[readNumber - 1][trimmedBases] += 1;
    }

    AllReadsStats operator+=(const AllReadsStats &rhs)
    {
        clustersCount += rhs.clustersCount;

        size_t numReads = rhs.trimLengthCountForReads.size();
        trimLengthCountForReads.resize(numReads);

        for (size_t i = 0; i < numReads; ++i)
        {
            size_t readLength = rhs.trimLengthCountForReads[i].size();
            trimLengthCountForReads[i].resize(readLength);

            for (size_t j = 0; j < readLength; ++j)
            {
                trimLengthCountForReads[i][j] += rhs.trimLengthCountForReads[i][j];
            }
        }

        return *this;
    }

    void reset()
    {
        clustersCount = 0;
        BOOST_FOREACH(std::vector<size_t>& trimLengthCount, std::make_pair(trimLengthCountForReads.begin(), trimLengthCountForReads.end()))
        {
            BOOST_FOREACH (size_t& count, std::make_pair(trimLengthCount.begin(), trimLengthCount.end()))
            {
                count = 0;
            }
        }
    }
};

/// \brief Per-read statistics for a single tile.
struct ReadStats
{
    common::Yield yield;
    common::Yield yieldQ30;
    common::QualityScore qualityScoreSum;

    ReadStats()
    : yield(0)
    , yieldQ30(0)
    , qualityScoreSum(0)
    {}

    ReadStats &operator+=(const ReadStats &rhs)
    {
        yield += rhs.yield;
        yieldQ30 += rhs.yieldQ30;
        qualityScoreSum += rhs.qualityScoreSum;
        return *this;
    }

    const ReadStats operator+(const ReadStats &rhs) const
    {
       return ReadStats(*this) += rhs;
    }

    void reset()
    {
        yield = 0;
        yieldQ30 = 0;
        qualityScoreSum = 0;
    }

};

namespace detail
{
    static const size_t TileStatsCount = 2;
} // namespace detail

/// \brief Raw and PF stats for a single tile-barcode
template< typename T >
class TileStats : public boost::array<T, detail::TileStatsCount>
{
public:
    enum type {RAW=0,PF};

    static const char *str_type[detail::TileStatsCount];

public:
    TileStats( const boost::array<T, detail::TileStatsCount> &a )
    : boost::array<T, detail::TileStatsCount>(a) {}
    TileStats()
    : boost::array<T, detail::TileStatsCount>()  {}
    TileStats( const T &ts )
    : boost::array<T, detail::TileStatsCount>() { this->fill(ts); }

    virtual ~TileStats() {}

    TileStats<T> &operator+=( const TileStats<T> &rhs )
    {
        size_t length = this->size();
        for (size_t i = 0; i < length; ++i)
        {
            (*this)[i] += rhs[i];
        }
        return *this;
    }

    void reset()
    {
        for( typename boost::array< T, detail::TileStatsCount>::size_type i=0;
             i < boost::array<T, detail::TileStatsCount>::static_size;
             i++ )
        {
            (*this)[i].reset();
        }
    }
};


/// \brief stats for a single tile-barcode
struct TileBarcodeStats : public TileStats<AllReadsStats>
{
    TileBarcodeStats( const TileStats<AllReadsStats> &ts )          : TileStats<AllReadsStats>(ts)  {}
    TileBarcodeStats( common::ReadNumber readIndex = 0 )            : TileStats<AllReadsStats>()    {}
    TileBarcodeStats( const AllReadsStats &ts )                     : TileStats<AllReadsStats>(ts)  {}

    TileBarcodeStats( const TileBarcodeStats &rbs )                 : TileStats<AllReadsStats>(rbs) {}

    TileBarcodeStats &operator+=( const TileBarcodeStats &rhs )
    {
        TileStats<AllReadsStats>::operator+=(rhs);
        return *this;
    }
};


/// \brief stats for a single tile-barcode, annotated by read number
struct ReadBarcodeStats : public TileStats<ReadStats>
{
    layout::ReadInfosContainer::size_type read;

    ReadBarcodeStats( const TileStats<ReadStats> &tbs,
                      layout::ReadInfosContainer::size_type readIndex = 0 )
    : TileStats<ReadStats>(tbs), read(readIndex)  {}
    ReadBarcodeStats( common::ReadNumber readIndex = 0 )
    : TileStats<ReadStats>(), read(readIndex)     {}
    ReadBarcodeStats( const ReadStats &ts,
                      layout::ReadInfosContainer::size_type readIndex = 0 )
    : TileStats<ReadStats>(ts), read(readIndex)   {}

    ReadBarcodeStats(const ReadBarcodeStats &rbs) : TileStats<ReadStats>(rbs), read(rbs.read) {}

    ReadBarcodeStats &operator+=( const ReadBarcodeStats &rhs )
    {
        TileStats<ReadStats>::operator+=(rhs);
        return *this;
    }
};


} // namespace stats
} // namespace bcl2fastq


#endif // BCL2FASTQ_STATS_TILESTATS_HPP


