/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file SampleInfo.hh
 *
 * \brief Declaration of sample metadata.
 *
 * \author Marek Balint
 */


#ifndef BCL2FASTQ_LAYOUT_SAMPLEINFO_HH
#define BCL2FASTQ_LAYOUT_SAMPLEINFO_HH


#include <vector>
#include <string>

#include "common/Types.hh"
#include "layout/Barcode.hh"


namespace bcl2fastq
{
namespace layout
{


/// \brief Sample metadata.
class SampleInfo
{
public:

    /// \brief Sample ID type definition.
    typedef std::string SampleId;

    /// \brief Sample name type definition.
    typedef std::string SampleName;

    /// \brief Project type definition.
    typedef std::string Project;

public:

    /// \brief Label to use when no project name has been provided.
    static const std::string defaultProject;
    /// \brief Label to use when no sample id has been provided.
    static const std::string defaultId;
    /// \brief Label to use when no sample name has been provided.
    static const std::string defaultName;

    /// \brief Constructor.
    /// \param number Sample number.
    /// \param sampleId Sample ID.
    /// \param sampleName Sample name.
    /// \param project Project name.
    SampleInfo(
        common::SampleNumber number     = 0,
        SampleId             sampleId   = SampleInfo::defaultId,
        SampleName           sampleName = SampleInfo::defaultName,
        Project              project    = SampleInfo::defaultProject
    );

public:

    /// \brief Get sample number.
    /// \return Sample number.
    common::SampleNumber getNumber() const;

    /// \brief Set sample number.
    void setNumber(common::SampleNumber number) { number_ = number; }

    /// \brief Get sample ID.
    /// \return Sample ID.
    const SampleInfo::SampleId& getId() const { return sampleId_; }

    /// \brief Get sample name.
    /// \return Sample name.
    const SampleInfo::SampleName& getName() const { return sampleName_; }

    /// \brief Get Project.
    /// \return Project.
    const SampleInfo::Project& getProject() const { return project_; }

    /// \brief Get Barcode.
    /// \return The Barcode.
    const Barcode& getBarcode(BarcodesContainer::size_type idx) const;

    /// \return whether Sample Id is defined or not.
    bool hasId() const;

    /// \return whether Sample Name is defined or not.
    bool hasName() const;

    /// \return whether Project is defined or not.
    bool hasProject() const;

    /// \return whether sample has any barcodes or not.
    bool hasBarcodes() const;

    /// \brief Mask the barcode for the given index number
    /// \param index number.
    void maskBarcode(common::ReadNumber indexNumber);

public:

    const BarcodesContainer& getBarcodes() const { return barcodes_; }

public:

    /// \brief Add barcode.
    /// \param barcode Barcode to be added to sample.
    void addBarcode(const Barcode &barcode);

private:

    /// \brief Sample number.
    ///
    /// Source: SampleSheet.csv (0,1..n)
    common::SampleNumber number_;

    /// \brief Sample ID.
    ///
    /// Source: SampleSheet.csv
    SampleId sampleId_;

    /// \brief Sample name.
    ///
    /// Source: SampleSheet.csv
    SampleName sampleName_;

    /// \brief Project name.
    ///
    /// Source: SampleSheet.csv
    Project project_;

    /// \brief Barcodes.
    ///
    /// Source: SampleSheet.csv
    BarcodesContainer barcodes_;

private:
    friend std::ostream& operator<<(std::ostream& os, const SampleInfo &sampleInfo);

};


} // namespace layout
} // namespace bcl2fastq


#endif // BCL2FASTQ_LAYOUT_SAMPLEINFO_HH


