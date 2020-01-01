/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file BarcodeTranslationTable.cpp
 *
 * \brief Implementation of barcode translation table.
 *
 * \author Roman Petrovski
 * \author Marek Balint
 */


#include <algorithm>
#include <utility>

#include <boost/assign/list_of.hpp>

#include "common/Debug.hh"
#include "common/Exceptions.hh"
#include "layout/BarcodeTranslationTable.hh"


namespace bcl2fastq
{
namespace layout
{


BarcodeTranslationTable::BarcodeTranslationTable(
    const layout::LaneInfo::SampleInfosContainer::const_iterator samplesBegin,
    const layout::LaneInfo::SampleInfosContainer::const_iterator samplesEnd,
    const config::BarcodeMismatchCountsContainer& componentMaxMismatches
)
: barcodeSampleMappings_(this->initBarcodeSampleMappings(samplesBegin, samplesEnd, componentMaxMismatches))
{
}

BarcodeTranslationTable::SampleMetadata BarcodeTranslationTable::translateBarcode(
    const layout::Barcode &barcode,
    unsigned int &mismatches
) const
{
    const BarcodeSampleMapping mappingToTranslate(barcode, SampleMetadata());
    const BarcodeSampleMappingsContainer::const_iterator mapping = std::lower_bound(
        barcodeSampleMappings_.begin(),
        barcodeSampleMappings_.end(),
        mappingToTranslate,
        &BarcodeTranslationTable::compareBarcodeSampleMappings
    );

    if (mapping == barcodeSampleMappings_.end() || mapping->barcode != mappingToTranslate.barcode)
    {
        return SampleMetadata();
    }
    else
    {
        //BCL2FASTQ_LOG(common::LogLevel::TRACE) << "LOOKUP:" << *(mapping.first) << std::endl;
        mismatches = mapping->mismatches;
        return mapping->sample;
    }
}

const BarcodeTranslationTable::BarcodeSampleMappingsContainer BarcodeTranslationTable::initBarcodeSampleMappings(
    const layout::LaneInfo::SampleInfosContainer::const_iterator samplesBegin,
    const layout::LaneInfo::SampleInfosContainer::const_iterator samplesEnd,
    const config::BarcodeMismatchCountsContainer& componentMaxMismatches
)
{
    BarcodeSampleMappingsContainer barcodeMappings;
    for (layout::LaneInfo::SampleInfosContainer::const_iterator samplesIter = samplesBegin; samplesIter != samplesEnd; ++samplesIter)
    {
        const auto& barcodes = samplesIter->getBarcodes();
        for (size_t i = 0; i < barcodes.size(); ++i)
        {
            SampleMetadata sample;
            sample.sampleIndex_ = (samplesIter - samplesBegin);
            sample.barcodeIndex_ = i;
            barcodeMappings.push_back(BarcodeSampleMappingsContainer::value_type(barcodes[i], sample));
        }
    }

    BarcodeSampleMappingsContainer ret;

    BOOST_FOREACH(const BarcodeSampleMappingsContainer::value_type &barcodeMapping, barcodeMappings)
    {
        /// \todo Features: Special handling for sample No.0. will be needed in case, when sample No.0 will have some value other than '********-********', which is currently hardcoded in LaneInfo.cpp.

        BarcodeTranslationTable::insertBarcodeMapping(
            barcodeMapping,
            barcodeMapping.barcode.componentsBegin(),
            0,
            componentMaxMismatches.front(),
            componentMaxMismatches.begin(),
            componentMaxMismatches.end(),
            ret
        );
    }

    std::sort(
        ret.begin(),
        ret.end(),
        &BarcodeTranslationTable::compareBarcodeSampleMappings
    );

    return ret;
}

void BarcodeTranslationTable::insertBarcodeMapping(
    const BarcodeSampleMapping &barcodeMapping,
    layout::Barcode::ComponentsContainer::const_iterator component,
    layout::Barcode::Component::SizeType position,
    unsigned int allowedMismatchesLeft,
    std::vector<std::size_t>::const_iterator componentMaxMismatchesIter,
    std::vector<std::size_t>::const_iterator componentMaxMismatchesEnd,
    BarcodeSampleMappingsContainer &container
)
{
    BCL2FASTQ_LOG(common::LogLevel::TRACE) << "\t[barcode table] " << barcodeMapping << std::endl;
    if(component != barcodeMapping.barcode.componentsEnd())
    {
        // simulate mismatches
        if (allowedMismatchesLeft > 0) {
            const layout::Barcode::Component::SizeType componentLength = component->getLength();
            for (layout::Barcode::Component::SizeType mismatchedPosition = position; mismatchedPosition < componentLength; ++mismatchedPosition)
            {
                // magic number 4: There are 5 distinct values ('A', 'C', 'G', 'T', 'N'), therefore there are 4 distinct mismatches
                std::vector<layout::Barcode> mismatchedBarcodes;
                mismatchedBarcodes.reserve(4);
                barcodeMapping.barcode.generateMismatches(component, mismatchedPosition, std::back_inserter(mismatchedBarcodes));

                BOOST_FOREACH (const layout::Barcode &mismatchedBarcode, mismatchedBarcodes)
                {
                    const BarcodeSampleMapping mismatchedBarcodeMapping(mismatchedBarcode,
                                                                        barcodeMapping.sample,
                                                                        *componentMaxMismatchesIter - allowedMismatchesLeft + 1);
                    BarcodeTranslationTable::insertBarcodeMapping(
                        mismatchedBarcodeMapping,
                        mismatchedBarcodeMapping.barcode.componentsBegin() + (component - barcodeMapping.barcode.componentsBegin()),
                        mismatchedPosition + 1,
                        allowedMismatchesLeft - 1,
                        componentMaxMismatchesIter,
                        componentMaxMismatchesEnd,
                        container
                    );
                }
            }
        }

        // move to next component (do not change current component)
        std::vector<std::size_t>::const_iterator nextMismatchesIter = (componentMaxMismatchesIter + 1 != componentMaxMismatchesEnd
            ?
            componentMaxMismatchesIter + 1
            : 
            componentMaxMismatchesIter
        );
        BarcodeTranslationTable::insertBarcodeMapping(
            barcodeMapping,
            component + 1,
            0,
            *nextMismatchesIter,
            nextMismatchesIter,
            componentMaxMismatchesEnd,
            container
        );
    }
    else
    {
        // actual insertion
        container.push_back(barcodeMapping);
    }
}

bool BarcodeTranslationTable::compareBarcodeSampleMappings(const BarcodeSampleMapping &lhs, const BarcodeSampleMapping &rhs)
{
    return lhs.barcode < rhs.barcode;
}


std::ostream& operator<<(std::ostream& os, const BarcodeTranslationTable::BarcodeSampleMapping &bsm)
{
    os << bsm.barcode;
    if( 0 != bsm.sample.sampleIndex_ )
    {
        os << std::string(" => sample #") << bsm.sample.sampleIndex_;
    } else {
        os << std::string(" : *");
    }
    os << " : " << bsm.mismatches;
    return os;
}


} // namespace conversion
} // namespace bcl2fastq


