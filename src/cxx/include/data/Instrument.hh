/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file Instrument.hh
 *
 * \brief Declaration of Instrument file.
 *
 * \author Daniel Berard 
 */

#ifndef BCL2FASTQ_DATA_INSTRUMENT_HH
#define BCL2FASTQ_DATA_INSTRUMENT_HH

#include "data/BclFileReader.hh"
#include "common/Types.hh"

// include/common/Types.hh:enum class TileAggregationMode { NON_AGGREGATED, AGGREGATED, CBCL };

namespace bcl2fastq {

namespace layout {
class TileInfo;
}

namespace data {

class Instrument
{
public:

    Instrument(common::TileAggregationMode tileAggregationMode)
    : tileAggregationMode_(tileAggregationMode)
    {
    }

    virtual ~Instrument() { }

    static bool isFilterFileAggregated(
        common::TileAggregationMode mode) 
    {
        return (mode == common::TileAggregationMode::AGGREGATED);
    }

    bool isFilterFileAggregated() const
    {
        return isFilterFileAggregated(tileAggregationMode_); 
    }

    static bool isPositionFileAggregated(common::TileAggregationMode mode)
    {
        return (mode != common::TileAggregationMode::NON_AGGREGATED);
    }

    bool isPositionFileAggregated() const 
    {
        return isPositionFileAggregated(tileAggregationMode_); 
    }

private:

   common::TileAggregationMode tileAggregationMode_;

}; // end class Instrument

} // namespace data
} // namespace bcl2fastq

#endif // BCL2FASTQ_DATA_INSTRUMENT_HH


