/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file BarcodeCollisionDetector.hh
 *
 * \brief Declaration of BarcodeCollisionDetector
 *
 * \author Aaron Day
 */


#ifndef BCL2FASTQ_LAYOUT_BARCODE_COLLISION_DETECTOR_HH
#define BCL2FASTQ_LAYOUT_BARCODE_COLLISION_DETECTOR_HH

#include "layout/Barcode.hh"
#include "common/Exceptions.hh"
#include "config/Bcl2FastqOptions.hh"

#include <boost/noncopyable.hpp>
#include <vector>

namespace bcl2fastq
{
namespace layout
{

class BarcodeCollisionError : public common::RuntimeError
{
public:
    BarcodeCollisionError(const std::string& message);
    BarcodeCollisionError(const BarcodeCollisionError& other);
};

class BarcodeCollisionDetector : private boost::noncopyable
{
public:
    BarcodeCollisionDetector(config::BarcodeMismatchCountsContainer& componentMaxMismatches,
                             bool autoSetToZeroMismatches = false);

    void validateBarcode(BarcodesContainer barcodes);

private:
    void validateNewBarcodeSizes(const BarcodesContainer& barcodes,
                                 const std::vector<std::size_t>& barcodeComponentSizes) const;

    void validateNewBarcodeSizesAgainstExisting(const std::vector<std::size_t>& barcodeComponentSizes) const;

    void validateBarcodeSizes(const BarcodesContainer& barcodes);

    void validateBarcodes(const BarcodesContainer::value_type& barcode1,
                          const BarcodesContainer::value_type& barcode2);

    void handleCollision(const BarcodesContainer::value_type& barcode1,
                         const BarcodesContainer::value_type& barcode2);

    size_t getNumMismatches(const Barcode::Component& barcode1,
                            const Barcode::Component& barcode2) const;

    bool                                    autoSetToZeroMismatches_;
    bool                                    triggerZeroMismatches_;
    config::BarcodeMismatchCountsContainer& componentMaxMismatches_;
    BarcodesContainer                       barcodes_;
};

} // namespace layout
} // namespace bcl2fastq

#endif // BCL2FASTQ_LAYOUT_BARCODE_COLLISION_DETECTOR_HH
