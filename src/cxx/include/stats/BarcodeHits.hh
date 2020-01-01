/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file BarcodeHits.hh
 *
 * \brief Declaration of Barcode Hits.
 *
 * \author Mauricio Varea
 */


#ifndef BCL2FASTQ_STATS_BARCODEHITS_HH
#define BCL2FASTQ_STATS_BARCODEHITS_HH

#include <algorithm>
#include <numeric>

#include <boost/foreach.hpp>

#include "common/Types.hh"
#include "common/Debug.hh"


namespace bcl2fastq
{
namespace stats
{


/// \brief Keeping track of barcodes
class BarcodeHits
{
    typedef std::pair< layout::Barcode, uint64_t >      Hits;
    typedef std::map< Hits::first_type, Hits::second_type >  Collector;
public:
    typedef std::vector< Hits >                              Popular;

public:
    BarcodeHits() : raw_(), popular_(), numBarcodes_(0) {}

    void record(const Collector::key_type &sequence)
    {
        static const uint64_t UNKNOWN_BARCODE_SAMPLE_FREQ = 20;

        ++numBarcodes_;
        if (numBarcodes_ % UNKNOWN_BARCODE_SAMPLE_FREQ == 0)
        {
            raw_[sequence] += UNKNOWN_BARCODE_SAMPLE_FREQ;
        }
    }

    void merge( const BarcodeHits &rhs )
    {
        BOOST_FOREACH( const Collector::value_type &barcodeHit, rhs.raw_ )
        {
            std::pair< Collector::iterator, bool > ret;
            ret = raw_.insert ( barcodeHit );
            if (ret.second==false)
            {
                ret.first->second += barcodeHit.second;
            }
        }

        numBarcodes_ += rhs.numBarcodes_;
    }

    void reset()
    {
        raw_.clear();
        Collector().swap(raw_);
        numBarcodes_ = 0;
    }

    void findPopular( size_t size, common::LaneNumber laneNumber )
    {
        popular_.reserve(size);
        BOOST_FOREACH( const Collector::value_type &barcodeHit, raw_ )
        {
            Popular::iterator insertIt = std::lower_bound( popular_.begin(),
                                                           popular_.end(),
                                                           std::make_pair(layout::Barcode(), barcodeHit.second),
                                                           BarcodeHits::orderByHits );

            if (popular_.size() < size)
            {
                popular_.insert(insertIt, barcodeHit);
            }
            else if (popular_.end() != insertIt)
            {
                popular_.insert(insertIt, barcodeHit);
                popular_.pop_back();
            }
        }
        const std::string unknownBarcodeCount( popular_.size()
                                             ? ("the " + boost::lexical_cast<std::string>(popular_.size()) + " most popular")
                                             : "no");
        BCL2FASTQ_LOG(common::LogLevel::DEBUG) << "Taking " << unknownBarcodeCount << " unknown barcodes "
                                               << "(out of " << raw_.size() << ") "
                                               << "in lane '" << laneNumber << "'" << std::endl;
        // destroying raw_ once popular_ has been populated, as this takes too much space!
        reset();
    }

    const Popular &getPopularBarcodes() const
    {
        //BCL2FASTQ_ASSERT_MSG(!popular_.empty(), "Cannot retrieve popular barcodes");
        return popular_;
    }

    Hits::second_type sum() const
    {
        return std::accumulate( raw_.begin(), raw_.end(), Hits::second_type(0), BarcodeHits::sumHits );
    }

    static bool orderBySequence( const Hits &lhs, const Hits &rhs )
    {
        return lhs.first < rhs.first;
    }
    static bool orderByHits( const Hits &lhs, const Hits &rhs )
    {
        // we want the Greatest Hits on top!
        return lhs.second > rhs.second;
    }

private:
    static Hits::second_type sumHits( Hits::second_type lhs, const Hits &rhs )
    {
        // we want the Greatest Hits on top!
        return lhs + rhs.second;
    }

private:
    /// \brief count of all barcodes
    Collector raw_;

    /// \brief subset with the most popular barcodes
    Popular   popular_;

    /// \brief Total number of barcodes recorded
    uint64_t numBarcodes_;
};


} // namespace stats
} // namespace bcl2fastq


#endif // BCL2FASTQ_STATS_BARCODEHITS_HH


