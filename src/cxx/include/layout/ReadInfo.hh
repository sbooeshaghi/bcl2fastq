/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file ReadInfo.hh
 *
 * \brief Declaration of read metadata.
 *
 * \author Marek Balint
 */


#ifndef BCL2FASTQ_LAYOUT_READINFO_HH
#define BCL2FASTQ_LAYOUT_READINFO_HH


#include <vector>
#include <string>
#include <set>

#include "common/Types.hh"
#include "layout/CycleInfo.hh"


namespace bcl2fastq {


namespace layout {


/// \brief Read metadata.
class ReadInfo
{
public:

    /// \brief Cycle metadata container type definition.
    typedef std::vector<CycleInfo> CycleInfosContainer;

    /// \brief Adapters container type definition.
    typedef std::vector<std::string> AdaptersContainer;

public:

    /// \brief Constructor.
    /// \param number Read number.
    /// \param readType Type of read (DATA, INDEX).
    /// \param cyclesToUse cycles to use (cycles to convert for data reads, index bases used in demultiplexing for index reads).
    /// \param unmaskedCycles All cycles found in RunInfo.xml, including those that were subsequently masked.
    /// \param umiCycles Cycles to use for the UMI, which is placed in the headers.
    ReadInfo(common::ReadNumber        number,
             common::ReadType          readType,
             const common::CycleRange& cyclesToUse,
             const common::CycleRange& unmaskedCycles,
             const common::CycleRange& umiCycles);

public:

    /// \brief Get read number.
    /// \return Read number.
    common::ReadNumber getNumber() const { return number_; }

    /// \brief Determine if this is an index read.
    /// \retval true This is index read.
    /// \retval false This is not index read.
    bool isIndexRead() const { return readType_ == common::ReadType::INDEX; }

    /// \brief Determine if this is a data read.
    /// \retval true This is a data read.
    /// \retval false This is not a data read.
    bool isDataRead() const { return readType_ == common::ReadType::DATA; }

    /// \brief Get the read type
    /// \retval Read type.
    common::ReadType getReadType() const { return readType_; }

public:

    /// \brief Get cycles
    /// \return Unmasked cycles
    const ReadInfo::CycleInfosContainer& cycleInfos() const { return cycleInfos_; }

    /// \brief Get unmasked cycles
    /// \return Unmasked cycles
    const ReadInfo::CycleInfosContainer& unmaskedCycleInfos() const { return unmaskedCycleInfos_; }

    /// \brief Get cycles to load from BCL files
    /// \return Cycles to load
    const std::set<CycleInfo>& cyclesToLoad() const { return cyclesToLoad_; }

    /// \brief Get beginning of adapters to be masked.
    /// \return Iterator to beginning of adapters to be masked.
    ReadInfo::AdaptersContainer::const_iterator maskAdaptersBegin() const { return maskAdapters_.begin(); }

    /// \brief Get end of adapters to be masked.
    /// \return Iterator to end of adapters to be trimmed.
    ReadInfo::AdaptersContainer::const_iterator maskAdaptersEnd() const { return maskAdapters_.end(); }

    /// \brief Get beginning of adapters to be trimmed.
    /// \return Iterator to beginning of adapters to be trimmed.
    ReadInfo::AdaptersContainer::const_iterator trimAdaptersBegin() const { return trimAdapters_.begin(); }

    /// \brief Get end of adapters to be trimmed.
    /// \return Iterator to end of adapters to be trimmed.
    ReadInfo::AdaptersContainer::const_iterator trimAdaptersEnd() const { return trimAdapters_.end(); }

public:

    /// \brief Get the number of cycles to load from the BCL files. Note that we skip cycles we won't use.
    /// \return Number of cycles to load
    size_t getNumCyclesToLoad() const { return cyclesToLoad_.size(); }

    /// \brief Get the BCL buffer offset for this read
    /// \return BCL buffer offset
    common::CycleNumber getBclBufferOffset() const { return bclBufferOffset_; }

    /// \brief Get the first cycle to load from the BCL file.
    /// \return The cycle number of the first cycle.
    common::CycleNumber getFirstLoadedCycle() const { return cyclesToLoad_.begin()->getNumber(); }

    /// \brief Get the cycle numbers of the UMI cycles.
    /// \return The range of cycles to use for UMI.
    const common::CycleRange& getUmiCycles() const { return umiCycleRange_; }

    /// \brief Add adapter to be masked.
    /// \param maskAdapter Adapter to be masked.
    void addMaskAdapter(const std::string &maskAdapter) { maskAdapters_.push_back(maskAdapter); }

    /// \brief Add adapter to be trimmed.
    /// \param trimAdapter Adapter to be trimmed.
    void addTrimAdapter(const std::string &trimAdapter) { trimAdapters_.push_back(trimAdapter); }

    /// \brief Return true if a FASTQ file should be created for this read.
    /// \return True if a FASTQ file should be created.
    bool shouldCreateFastqForRead(bool createFastqsForIndexReads) const { return !cycleInfos_.empty() && (isDataRead() || (createFastqsForIndexReads && isIndexRead())); }

private:

    /// \brief Add the UMI cycles to the read.
    /// \param The range of UMI cycles to add.
    void addUmiBases(common::CycleRange umiCycleRange);

    /// \brief Add cycle metadata.
    /// \param cycleInfo Cycle metadata to be added to read.
    void addCycle(const CycleInfo &cycleInfo) { cycleInfos_.push_back(cycleInfo); cyclesToLoad_.insert(cycleInfo); }

    /// \brief Add unmasked cycle metadata.
    /// \param cycleInfo Cycle metadata to be added to read.
    void addUnmaskedCycle(const CycleInfo &cycleInfo);

    /// \brief Read number.
    ///
    /// Source: RunInfo.xml (1..n)
    common::ReadNumber number_;

    /// \brief Type of read (DATA, INDEX, or UMI).
    ///
    /// Source: RunInfo.xml, SampleSheet.csv, command line
    common::ReadType readType_;

    /// \brief Cycles metadata container. (May drop masked cycles)
    CycleInfosContainer cycleInfos_;

    /// \brief Unmasked cycles.
    CycleInfosContainer unmaskedCycleInfos_;

    /// \brief Cycle to load from BCL files.
    ///
    /// \Source: RunInfo.xml, --use-bases-mask, SampleSheet.csv.
    std::set<CycleInfo> cyclesToLoad_;

    /// \brief Adapters to be masked.
    ///
    /// Source: SampleSheet.csv
    AdaptersContainer maskAdapters_;

    /// \brief Adapters to be trimmed.
    ///
    /// Source: SampleSheet.csv
    AdaptersContainer trimAdapters_;

    /// \brief UMI cycle range
    ///
    /// Source: SampleSheet.csv
    common::CycleRange umiCycleRange_;

    /// \brief BCL buffer offset
    ///
    /// Source: SampleSheet.csv, --use-bases-mask, RunInfo.xml
    common::CycleNumber bclBufferOffset_;
};

/// \brief Read metadata container type definition.
typedef std::vector<ReadInfo> ReadInfosContainer;

} // namespace layout


} // namespace bcl2fastq


#endif // BCL2FASTQ_LAYOUT_READINFO_HH


