/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file TileStats.cpp
 *
 * \brief Tile and Tile-Barcode Stats.
 *
 * \author Mauricio Varea
 */

#include "stats/TileStats.hpp"

namespace bcl2fastq
{
namespace stats
{

// Do not be tempted to use [detail::TileStatsCount]
//    keeping the size hardcoded forces us to come back to this part of the code every time we do a change
//    in the header and forget to change it here as well
template<> const char *TileStats<AllReadsStats>::str_type[detail::TileStatsCount] = {"Raw","Pf"};
template<> const char *TileStats<ReadStats>::str_type[detail::TileStatsCount]     = {"Raw","Pf"};


} //namespace stats
} //namespace bcl2fastq

