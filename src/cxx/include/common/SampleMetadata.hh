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
 * \author Aaron Day
 */


#ifndef BCL2FASTQ_COMMON_SAMPLE_METADATA_HH
#define BCL2FASTQ_COMMON_SAMPLE_METADATA_HH


#include "common/Types.hh"

#include <string>
#include <vector>

namespace bcl2fastq {


namespace common {

/// \brief Sample meta data.
struct SampleMetadata
{
public:

    /// \brief Sample barcodes container type definition.
    typedef std::vector<std::vector<std::string> > BarcodesContainer;

public:

    /// \brief Sample ID.
    std::string id_;

    /// \brief Sample name.
    std::string name_;

    /// \brief Project name.
    std::string project_;

    /// \brief Lane number.
    std::vector<common::LaneNumber> lanes_;

    /// \brief Sample barcodes.
    BarcodesContainer barcodes_;
};

/// \brief Sample metadta container type definition.
typedef std::vector<SampleMetadata> SampleMetadataContainer;

} // namespace common


} // namespace bcl2fastq


#endif // BCL2FASTQ_COMMON_SAMPLE_METADATA_HH


