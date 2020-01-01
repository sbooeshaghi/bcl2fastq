/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file BarcodeTranslationTable.hh
 *
 * \brief Declaration of barcode translation table.
 *
 * \author Roman Petrovski
 * \author Marek Balint
 * \author Mauricio Varea
 */


#ifndef BCL2FASTQ_LAYOUT_BARCODETRANSLATIONTABLE_HH
#define BCL2FASTQ_LAYOUT_BARCODETRANSLATIONTABLE_HH


#include <vector>
#include <string>

#include <boost/filesystem.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include "common/Exceptions.hh"
#include "config/Bcl2FastqOptions.hh"
#include "layout/Barcode.hh"
#include "layout/SampleInfo.hh"
#include "layout/LaneInfo.hh"


namespace bcl2fastq
{
namespace layout
{


/// \brief Barcode translation table.
/// \todo Refactoring: This class is just a pretty-printed copy-paste form old converter, forced into this code.
/// \todo Features: Make number of allowed mismatches per component configurable (currently hardcoded to 2).
class BarcodeTranslationTable : private boost::noncopyable
{
public:

    /// \brief Sample meta data.
    struct SampleMetadata
    {
        /// \brief Sample index.
        layout::LaneInfo::SampleInfosContainer::size_type sampleIndex_;

        /// \brief Barcode index.
        BarcodesContainer::size_type barcodeIndex_;

        SampleMetadata( layout::LaneInfo::SampleInfosContainer::size_type sampleIndex,
                        BarcodesContainer::size_type barcodeIndex = 0)
        : sampleIndex_(sampleIndex), barcodeIndex_(barcodeIndex) {}

        SampleMetadata()  ///< \brief Default sample
        : sampleIndex_(0), barcodeIndex_(0) {}

        /// \param lhs Left-hand-side argument.
        /// \param rhs Right-hand-side argument.
        /// \retval true Left-hand-side argument is less than right-hand-side argument.
        /// \retval false Left-hand-side argument is not less than right-hand-side argument.
        bool operator<( const SampleMetadata &rhs ) const
        {
            if (this->sampleIndex_ != rhs.sampleIndex_)
            {
                return this->sampleIndex_ < rhs.sampleIndex_;
            }
            return this->barcodeIndex_ < rhs.barcodeIndex_;
        }

        bool operator==( const SampleMetadata &rhs ) const
        {
            return (this->sampleIndex_ == rhs.sampleIndex_) && (this->barcodeIndex_ == rhs.barcodeIndex_);
        }

        bool operator!=( const SampleMetadata &rhs ) const  { return !(*this == rhs); }
    };

public:

    /// \brief Constructor.
    /// \param samplesBegin Beginning of samples.
    /// \param samplesEnd End of samples.
    /// \param componentMaxMismatches Maximum number of allowed mismatches per barcode component.
    BarcodeTranslationTable(
        const layout::LaneInfo::SampleInfosContainer::const_iterator samplesBegin,
        const layout::LaneInfo::SampleInfosContainer::const_iterator samplesEnd,
        const config::BarcodeMismatchCountsContainer& componentMaxMismatches
    );

public:

    /// \brief Translate barcode to sample.
    /// \param barcode Barcode to be translated.
    /// \param mismatches Return number of mismatches w.r.t. perfect barcode
    /// \return Sample associated with the barcode.
    BarcodeTranslationTable::SampleMetadata translateBarcode(
        const layout::Barcode &barcode,
        unsigned int &mismatches
    ) const;

private:

    /// \brief Barcode to sample mapping type definition.
    struct BarcodeSampleMapping
    {
        layout::Barcode barcode;
        SampleMetadata sample;
        unsigned int mismatches;

        BarcodeSampleMapping(const layout::Barcode &b, const SampleMetadata &s, unsigned int m=0)
        : barcode(b), sample(s), mismatches(m) {}
    };

    /// \brief Barcode to sample mappings container type definition.
    typedef std::vector<BarcodeSampleMapping> BarcodeSampleMappingsContainer;

    /// \brief Barcode to sample mappings range definition.
    typedef std::pair<BarcodeSampleMappingsContainer::const_iterator,
                      BarcodeSampleMappingsContainer::const_iterator> BarcodeSampleMappingsRange;
private:

    /// \brief Initialize barcode to sample mappings.
    /// \param samplesBegin Beginning of samples.
    /// \param samplesEnd End of samples.
    /// \param componentMaxMismatches Maximum number of allowed mismatches per barcode component.
    /// \return Initialized barcode to sample mappings container.
    static const BarcodeTranslationTable::BarcodeSampleMappingsContainer initBarcodeSampleMappings(
        const layout::LaneInfo::SampleInfosContainer::const_iterator samplesBegin,
        const layout::LaneInfo::SampleInfosContainer::const_iterator samplesEnd,
        const config::BarcodeMismatchCountsContainer& componentMaxMismatches
    );

    /// \brief Insert barcode to sample mapping into container.
    /// \param barcodeMapping Barcode to sample mapping.
    /// \param component Component being currently processed.
    /// \param position Position of component, where mismatches are to be simulated.
    /// \param allowedMismatchesLeft Number of mismatches left for current component.
    /// \param componentMaxMismatchesIter Iterator to maximum number of allowed mismatches for current component.
    /// \param componentMaxMismatchesEnd End iterator to maximum number of allowed mismatches per barcode component.
    /// \param container Container to insert to.
    static void insertBarcodeMapping(
        const BarcodeSampleMapping &barcodeMapping,
        layout::Barcode::ComponentsContainer::const_iterator component,
        layout::Barcode::Component::SizeType position,
        unsigned int allowedMismatchesLeft,
        std::vector<std::size_t>::const_iterator componentMaxMismatchesIter,
        std::vector<std::size_t>::const_iterator componentMaxMismatchesEnd,
        BarcodeSampleMappingsContainer &container
    );

    /// \brief Compare barcode to sample mappings.
    /// \param lhs Left-hand-side argument.
    /// \param rhs Right-hand-side argument.
    /// \retval true Left-hand-side argument is less than right-hand-side argument.
    /// \retval false Left-hand-side argument is not less than right-hand-side argument.
    static bool compareBarcodeSampleMappings(const BarcodeSampleMapping &lhs, const BarcodeSampleMapping &rhs);

private:

    /// \brief Barcode to sample mappings.
    const BarcodeSampleMappingsContainer barcodeSampleMappings_;

    friend std::ostream& operator<<(std::ostream& os, const BarcodeTranslationTable::BarcodeSampleMapping &bsm);
};


} // namespace conversion
} // namespace bcl2fastq


#endif // BCL2FASTQ_LAYOUT_BARCODETRANSLATIONTABLE_HH


