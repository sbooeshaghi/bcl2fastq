/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file Types.hh
 *
 * \brief Declaration of common data types
 *
 * \author Marek Balint
 */


#ifndef BCL2FASTQ_COMMON_TYPES_HH
#define BCL2FASTQ_COMMON_TYPES_HH

#include "common/Debug.hh"
#include <boost/filesystem.hpp>
#include <cstddef>
#include <stdint.h>
#include <ostream>
#include <map>
#include <vector>

namespace bcl2fastq {
namespace common {

enum class NumBasesPerByte : uint32_t { ONE=1, TWO=2, FOUR=4 };

/// \brief Read type definition.
enum ReadType
{
    DATA,
    INDEX,
    UMI
};


/// \brief Raw data buffer type definition.
class RawDataBuffer : public std::vector<char>
{
public:
    RawDataBuffer() : std::vector<char>() { }
    RawDataBuffer(size_t sz) : std::vector<char>(sz) { }

    virtual ~RawDataBuffer() { }

    boost::filesystem::path path_;
};

template<typename T>
void resizeAndRealloc(T& container, size_t size)
{
/*    if (container.capacity() < size)
    {
        // We need to allocate memory. Reserve 20% more than we need to reduce the
        // liklihood we will need to do this again next time.
        container.reserve(size*1.2);
    }
*/
    container.resize(size);
}

template<typename T, typename TYPE>
void resizeAndRealloc(T& container, size_t size, TYPE defaultValue)
{
/*    if (container.capacity() < size)
    {
        // We need to allocate memory. Reserve 20% more than we need to reduce the
        // liklihood we will need to do this again next time.
        container.reserve(size*1.2);
    }
*/
    container.resize(size, defaultValue);
}

/// \brief Lane number type definition.
typedef std::size_t LaneNumber;

/// \brief Sample number type definition.
typedef std::size_t SampleNumber;

/// \brief Tile number type definition.
typedef std::size_t TileNumber;

/// \brief Read number type definition.
typedef std::size_t ReadNumber;

/// \brief Cycle number type definition.
typedef std::size_t CycleNumber;

/// \brief Cycle number range.definition.
typedef std::pair<CycleNumber, CycleNumber> CycleRange;

/// \brief Number of clusters per file/tile.
typedef uint32_t ClustersCount;

/// \brief Number of clusters
typedef uint64_t TotalClustersCount;

/// \brief Number of bases; throughput.
typedef uint64_t Yield;

/// \brief Number of bases; throughput.
typedef uint64_t QualityScore;

/// \brief Control flag type definition.
typedef uint16_t ControlFlag;

/// \brief TileAggregationMode type definition.
enum class TileAggregationMode { NON_AGGREGATED, AGGREGATED, CBCL };

/// brief type definition for map from lane and tile number to file name
typedef std::vector<std::map<common::TileNumber, boost::filesystem::path>> TileFileMap;

inline std::ostream& operator<<(std::ostream& os, ReadType readType)
{
    switch (readType)
    {
        case DATA:
            os << "data read";
            break;
        case INDEX:
            os << "index read";
            break;
        case UMI:
            os << "umi read";
            break;
        default:
            BCL2FASTQ_ASSERT_MSG(false, "Unrecognized ReadType");
        break;
    }

    return os;
}

} // namespace common
} // namespace bcl2fastq


#endif // BCL2FASTQ_DATA_TYPES_HH


