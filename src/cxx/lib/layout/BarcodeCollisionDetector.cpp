/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file BarcodeCollisionDetector.cpp
 *
 * \brief Implementation of BarcodeCollisionDetector.
 *
 * \author Aaron Day
 */

#include "layout/BarcodeCollisionDetector.hh"
#include "layout/Barcode.hh"
#include "common/Debug.hh"

namespace bcl2fastq
{
namespace layout
{

BarcodeCollisionError::BarcodeCollisionError(const std::string& message)
: common::RuntimeError(0, message)
{
}

BarcodeCollisionError::BarcodeCollisionError(const BarcodeCollisionError& other)
: common::RuntimeError(other)
{
}


BarcodeCollisionDetector::BarcodeCollisionDetector(config::BarcodeMismatchCountsContainer& componentMaxMismatches,
                                                   bool autoSetToZeroMismatches /*=false*/)
 : autoSetToZeroMismatches_(autoSetToZeroMismatches)
 , componentMaxMismatches_(componentMaxMismatches)
 , barcodes_()
{
}

void BarcodeCollisionDetector::validateBarcode(BarcodesContainer barcodes)
{
    validateBarcodeSizes(barcodes);

    // Compare each of the barcodes for a single sample against all barcodes for all other samples
    for (const auto& barcode2 : barcodes)
    {
        for (const auto& barcode1 : barcodes_)
        {
            validateBarcodes(barcode1, barcode2);
        }

        // There were no collisions, add this sample's barcodes to the list of all barcodes
        barcodes_.push_back(barcode2);
    }
}

void BarcodeCollisionDetector::validateBarcodeSizes(const BarcodesContainer& barcodes)
{
    BCL2FASTQ_ASSERT_MSG(!barcodes.empty(), "No barcodes for sample.");
    BCL2FASTQ_ASSERT_MSG(!barcodes.front().getComponents().empty(), "Barcode has 0 components");

    std::vector<std::size_t> barcodeComponentSizes;
    for (const auto& component : barcodes.front().getComponents())
    {
        barcodeComponentSizes.push_back(component.size());
    }

    validateNewBarcodeSizes(barcodes,
                            barcodeComponentSizes);

    validateNewBarcodeSizesAgainstExisting(barcodeComponentSizes);
}

void BarcodeCollisionDetector::validateNewBarcodeSizes(const BarcodesContainer& barcodes,
                                                       const std::vector<std::size_t>& barcodeComponentSizes) const
{
    size_t numComponents = barcodes.front().getComponents().size();

    for (const auto& barcode : barcodes)
    {
        if (numComponents != barcode.getComponents().size())
        {
            BOOST_THROW_EXCEPTION(BarcodeCollisionError("Barcodes have an unequal number of components."));
        }

        std::vector<std::size_t>::const_iterator barcodeComponentSize = barcodeComponentSizes.begin();
        for (const auto& component : barcode.getComponents())
        {
            if (component.size() != *barcodeComponentSize)
            {
                BOOST_THROW_EXCEPTION(BarcodeCollisionError("Barcodes components must have the same size."));
            }
            ++barcodeComponentSize;
        }
    }
}

void BarcodeCollisionDetector::validateNewBarcodeSizesAgainstExisting(const std::vector<std::size_t>& barcodeComponentSizes) const
{
    if (!barcodes_.empty())
    {
        if (barcodeComponentSizes.size() != barcodes_.back().getComponents().size())
        {
            BOOST_THROW_EXCEPTION(BarcodeCollisionError("Barcodes have an unequal number of components."));
        }

        std::vector<std::size_t>::const_iterator barcodeComponentSize = barcodeComponentSizes.begin();
        for (const auto& component : barcodes_.back().getComponents())
        {
            if (component.size() != *barcodeComponentSize)
            {
                BOOST_THROW_EXCEPTION(BarcodeCollisionError("Barcodes components must have the same size."));
            }
            ++barcodeComponentSize;
        }
    }
}

void BarcodeCollisionDetector::validateBarcodes(const BarcodesContainer::value_type& barcode1,
                                                const BarcodesContainer::value_type& barcode2)
{
    size_t numComponentsWithCollision = 0;
    std::vector<std::size_t>::const_iterator componentMaxMismatchesIter = componentMaxMismatches_.begin();
    auto barcode2Iter = barcode2.componentsBegin();
    for (const auto& component : barcode1.getComponents())
    {
        // Multiply by 2 because both barcodes can have mismatches
        if (2* (*componentMaxMismatchesIter) >= getNumMismatches(component, *barcode2Iter))
        {
            ++numComponentsWithCollision;
        }

        ++barcode2Iter;
        if (componentMaxMismatchesIter + 1 != componentMaxMismatches_.end())
        {
            ++componentMaxMismatchesIter;
        }
    }

    if (numComponentsWithCollision == barcode1.getComponents().size())
    {
        handleCollision(barcode1, barcode2);
    }
}

void BarcodeCollisionDetector::handleCollision(const BarcodesContainer::value_type& barcode1,
                                               const BarcodesContainer::value_type& barcode2)
{
    std::stringstream errStr;

    if (barcode1 == barcode2)
    {
        errStr << "Identical barcode found for multiple samples in the sample sheet: " << barcode1;
    }
    else
    {
        errStr << "Barcode collision for barcodes: "
               << barcode1
               << ", "
               << barcode2;

        if (autoSetToZeroMismatches_)
        {
            for (auto& maxMismatches : componentMaxMismatches_)
            {
                maxMismatches = 0;
            }
            BCL2FASTQ_LOG(common::LogLevel::WARNING) << errStr.str()
                << ". Setting the number of allowed mismatches to 0." << std::endl;

            return;
        }
        else
        {
            errStr << std::endl << "By default, bcl2fastq allows 1 mismatch in each barcode."
                                   " Barcodes with too few mismatches are ambiguous ( less than 2 times the number of mismatches plus 1)."
                                   " To reduce the number of allowed mismatches, use the command line option: '--barcode-mismatches'."
                                   " Note that particularly for barcodes with only 1 mismatch, there is the danger that some reads will"
                                   " be written to the wrong sample due to errors in the barcode sequence.";
        }
    }

    BOOST_THROW_EXCEPTION(BarcodeCollisionError(errStr.str()));
}

size_t BarcodeCollisionDetector::getNumMismatches(const Barcode::Component& barcodeComponent1,
                                                  const Barcode::Component& barcodeComponent2) const
{
    size_t numMismatches = 0;

    size_t barcodeComponentSize = barcodeComponent1.size();
    for (size_t i = 0; i < barcodeComponentSize; ++i)
    {
        // We've already verified that the sizes are equal
        if (barcodeComponent1.getString()[i] != barcodeComponent2.getString()[i])
        {
            ++numMismatches;
        }
    }

    return numMismatches;
}

} // namespace bcl2fastq
} // namespace layout

