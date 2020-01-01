/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file LaneInfo.hh
 *
 * \brief Declaration of lane metadata.
 *
 * \author Marek Balint
 */


#ifndef BCL2FASTQ_LAYOUT_LANEINFO_HH
#define BCL2FASTQ_LAYOUT_LANEINFO_HH


#include <vector>

#include <boost/filesystem/path.hpp>
#include <boost/format.hpp>

#include "common/Types.hh"
#include "layout/SampleInfo.hh"
#include "layout/TileInfo.hh"
#include "layout/ReadInfo.hh"


namespace bcl2fastq {


namespace layout {


/// \brief Lane metadata.
class LaneInfo
{
public:

    /// \brief Sample metadata container type definition.
    typedef std::vector<SampleInfo> SampleInfosContainer;

    /// \brief Tile metadata container type definition.
    typedef std::vector<TileInfo> TileInfosContainer;

public:

    /// \brief Constructor.
    /// \param number Lane number.
    explicit LaneInfo(common::LaneNumber number);

public:

    /// \brief Get lane number.
    /// \return Lane number.
    common::LaneNumber getNumber() const { return number_; }

    /// \brief Get lane directory name.
    /// \return Lane directory name.
    boost::filesystem::path getDirName() const { return boost::filesystem::path((boost::format("L%03d") % number_).str()); }

public:

    /// \brief Remove the default sample.
    void removeDefaultSample();

    /// \brief Mask the barcode for the given index number
    /// \param index number.
    void maskBarcode(common::ReadNumber indexNumber);

    /// \brief Get beginning of lanes.
    /// \return Iterator to beginning of lanes.
    LaneInfo::SampleInfosContainer::const_iterator sampleInfosBegin() const { return sampleInfos_.begin(); }

    /// \brief Get end of lanes.
    /// \return Iterator to end of lanes.
    LaneInfo::SampleInfosContainer::const_iterator sampleInfosEnd() const { return sampleInfos_.end(); }

    /// \brief Get the sample info container.
    /// \return Sample info container.
    const SampleInfosContainer& getSampleInfos() const { return sampleInfos_; }

    /// \brief Sort sampleInfos by sample number.
    void sortSampleInfos();

    /// \brief Get the tile infos
    /// \return Container of tile infos
    const TileInfosContainer& getTileInfos() const { return tileInfos_; }

    /// \brief Get the read infos
    /// \return Read infos container
    const ReadInfosContainer& readInfos() const { return readInfos_; }

    /// \brief Get the read infos. (The functions returning iterators are silly.)
    /// \return Read infos.
    const ReadInfosContainer& getReadInfos() const { return readInfos_; }

    /// \brief Get the minimum trimmed read length.
    /// \return Minimum trimmed read length.
    size_t getMinimumTrimmedReadLength() const { return minimumTrimmedReadLength_; }

    /// \brief Get the maximum number of barcode mismatches.
    /// \return a vector of the maximum number of barcode mismatches. 1 entry for each index.
    const std::vector<size_t>& getMaxBarcodeMismatches() const { return maxBarcodeMismatches_; }

    /// \brief Get number of cycles to load from BCL files.
    /// \return Number of cycles.
    common::CycleNumber getNumCyclesToLoad() const;

public:

    /// \brief Get clusters count.
    /// \return Number of clusters on all tiles on this lane.
    /// \pre <tt>this->haveClustersCount() == true</tt>
    common::TotalClustersCount getClustersCount() const;

    /// \brief Check whether clusters count is known.
    /// \retval true Clusters count is known.
    /// \retval false Clusters count is not known.
    bool haveClustersCount() const;

public:

    /// \brief Add sample metadata.
    /// \param sampleInfo Sample metadata to be added to read.
    void addSample(const SampleInfo &sampleInfo) { sampleInfos_.push_back(sampleInfo); }

    /// \brief Add tile metadata.
    /// \param tileInfo Tile metadata to be added to read.
    void addTile(const TileInfo &tileInfo) { tileInfos_.push_back(tileInfo); }

    /// \brief Add read metadata.
    //  \param readInfo Read metadata to be added.
    void addRead(const ReadInfo& readInfo) { readInfos_.push_back(readInfo); }

    /// \brief Set the minimum trimmed read length.
    /// \param minimumTrimmedReadLength The minimum trimmed read length.
    void setMinimumTrimmedReadLength(size_t minimumTrimmedReadLength) { minimumTrimmedReadLength_ = minimumTrimmedReadLength; }

    /// \brief Set the maximum number of barcode mismatches.
    /// \param Vector of maximum number of barcode mismatches. 1 entry for each index.
    void setMaxBarcodeMismtaches(const std::vector<std::size_t>& maxBarcodeMismatches) { maxBarcodeMismatches_ = maxBarcodeMismatches; }

    /// \brief Overloaded comparison operator.
    bool operator<(const LaneInfo &laneInfo) const  { return number_ < laneInfo.number_; }

private:

    /// \brief Lane number.
    ///
    /// Source: RunInfo.xml (1..n)
    common::LaneNumber number_;

    /// \brief Samples metadata container.
    SampleInfosContainer sampleInfos_;

    /// \brief Tile metadata container.
    TileInfosContainer tileInfos_;

    /// \brief Reads metadata container.
    ReadInfosContainer readInfos_;

    /// \brief Minimum trimmed read length.
    size_t minimumTrimmedReadLength_;

    /// \brief Maximum number of mismatches allowed in a barcode.
    std::vector<std::size_t> maxBarcodeMismatches_;
};

/// \brief Lane metadata container type definition.
typedef std::vector<LaneInfo> LaneInfosContainer;


} // namespace layout


} // namespace bcl2fastq


#endif // BCL2FASTQ_LAYOUT_LANEINFO_HH


