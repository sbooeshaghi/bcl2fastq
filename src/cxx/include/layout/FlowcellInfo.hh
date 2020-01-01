/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file FlowcellInfo.hh
 *
 * \brief Declaration of flowcell metadata.
 *
 * \author Marek Balint
 * \author Mauricio Varea
 */


#ifndef BCL2FASTQ_LAYOUT_FLOWCELLINFO_HH
#define BCL2FASTQ_LAYOUT_FLOWCELLINFO_HH


#include "layout/RunInfoXml.hh"


namespace bcl2fastq {


namespace layout {


/// \brief Flowcell metadata.
class FlowcellInfo
{
public:

    /// \brief Get instrument.
    /// \return Instrument.
    const std::string& getInstrument() const { return instrument_; }

    /// \brief Set instrument.
    /// \param instrument Value to be set.
    void setInstrument(const std::string &instrument);

    /// \brief Get run number.
    /// \return Run number.
    const std::string& getRunNumber() const { return runNumber_; }

    /// \brief Set run number.
    /// \param runNumber Value to be set.
    void setRunNumber(const std::string &runNumber);

    /// \brief Get flowcell ID.
    /// \return Flowcell ID.
    const std::string& getFlowcellId() const { return flowcellId_; }

    /// \brief Set flowcell ID.
    /// \param flowcellId Value to be set.
    void setFlowcellId(const std::string &flowcellId);

    /// \brief Get aggregate tiles flag.
    /// \retval true Data for all tiles are aggregated into single file.
    /// \retval false There is separate file for each tile.
    common::TileAggregationMode getAggregateTilesMode() const;

    /// \brief Set value of aggregate tiles flag.
    /// \param aggregateTilesFlag Value to be set.
    void setAggregateTilesFlag(common::TileAggregationMode aggregateTilesFlag);

    /// \brief Get patterned flowcell flag.
    /// \retval true The clusters positions is known up-front.
    /// \retval false The clusters positions is not known up-front.
    bool isPatternedFlowcell() const;

    /// \brief Set value of patterned flowcell flag to true.
    void setPatternedFlowcell();

    /// \brief Reset value of patterned flowcell flag to false.
    void resetPatternedFlowcell();

    /// \brief Set the run id.
    /// \param Run Id to set.
    void setRunId(const std::string& runId) { runId_ = runId; }

    /// \brief Get the run id.
    /// \retval Run Id.
    const std::string& getRunId() const { return runId_; }

private:

    /// \brief Instrument.
    ///
    /// Source: RunInfo.xml
    std::string instrument_;

    /// \brief Run number.
    ///
    /// Source: RunInfo.xml
    std::string runNumber_;

    /// \brief Flowcell ID.
    ///
    /// Source: RunInfo.xml
    std::string flowcellId_;

    /// \brief Aggregate tiles flag.
    ///
    /// Source: command line argument --aggregated-tiles
    common::TileAggregationMode aggregateTilesFlag_;

    /// \brief Patterned Flowcell flag.
    ///
    /// Source: whether s.locs exist or not
    bool patternedFlowcell_;

    /// \brief Run ID.
    ///
    /// Source: RunInfo.xml
    std::string runId_;
};


} // namespace layout


} // namespace bcl2fastq


#endif // BCL2FASTQ_LAYOUT_FLOWCELLINFO_HH


