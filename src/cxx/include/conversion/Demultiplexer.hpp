/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file Demultiplexer.hpp
 *
 * \brief Declaration of Demultiplexer.
 *
 * \author Aaron Day
 */


#ifndef BCL2FASTQ_CONVERSION_DEMULTIPLEXER_HPP
#define BCL2FASTQ_CONVERSION_DEMULTIPLEXER_HPP


#include "conversion/Demultiplexer.hh"

namespace bcl2fastq {
namespace conversion {

template<common::NumBasesPerByte numBasesPerByte>
bool DemultiplexTask::execute(int32_t threadNum)
{
    for (data::BclBuffer::BclCycleContainer::size_type offset = offsetBegin_; offset < offsetEnd_; ++offset)
    {
        const bool filterFlag = inputBuffer_.filters_.at(offset).data_;

        layout::Barcode& barcode = outputBuffer_.at(offset).second;
        // Clear the barcode. This does not free memory, which is good.
        // We don't want to reallocate the string every time, which
        // is particularly bad since this is multithreaded and the memory allocation will block.
        barcode.reset(numIndexReads_);

        if (filterFlag || includeNonPfClusters_)
        {
            size_t i = 0;
            layout::ReadInfo::CycleInfosContainer::const_iterator::difference_type previousReadsLength = 0;
            for (const auto& readInfo : readInfos_)
            {
                size_t currentReadLength = readInfo.cycleInfos().size();
                if (readInfo.isIndexRead() && currentReadLength != 0)
                {
                    const FastqConstIterator<numBasesPerByte> begin = FastqConstIterator<numBasesPerByte>(inputBuffer_, offset) + previousReadsLength;
                    const FastqConstIterator<numBasesPerByte> end = FastqConstIterator<numBasesPerByte>(inputBuffer_, offset) + previousReadsLength + currentReadLength;

                    size_t dist = end-begin;
                    std::string& barcodeString = barcode.components_[i].getString();
                    barcodeString.resize(dist);
                    for (size_t j = 0; j < dist; ++j)
                    {
                        barcodeString[j] = convertBcl2FastqBase(begin[j]);
                    }
                    ++i;
                }
                previousReadsLength += readInfo.cyclesToLoad().size();
            }
            unsigned int mismatches(0);
            layout::BarcodeTranslationTable::SampleMetadata sampleMetadata = barcodeTranslationTable_.translateBarcode(barcode,mismatches);

            outputBuffer_.at(offset).first = sampleMetadata;
            demuxStats_[threadNum]->getStats().at(sampleMetadata.sampleIndex_)[sampleMetadata.barcodeIndex_][inputBuffer_.tileInfo_ - tileInfosBegin_].incrementBarcodeCount(mismatches);
        }
        else
        {
            outputBuffer_.at(offset).first = layout::BarcodeTranslationTable::SampleMetadata();
        }
    }
    return true;
}

} // namespace conversion
} // namespace bcl2fastq

#endif // BCL2FASTQ_CONVERSION_DEMULTIPLEXER_HPP

