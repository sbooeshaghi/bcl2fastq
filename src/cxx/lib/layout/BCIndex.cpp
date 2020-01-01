/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file BCIndex.cpp
 *
 * \brief Implementation of BCI helper.
 *
 * \author Marek Balint
 */


#include <sstream>
#include <utility>

#include <boost/filesystem/fstream.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>

#include "common/Debug.hh"
#include "common/Exceptions.hh"
#include "layout/BCIndex.hh"


namespace bcl2fastq {


namespace layout {


namespace detail {


/// \brief Generate path to BCI file.
/// \param inputDir Path to input directory.
/// \param laneNumber Lane number.
/// \return Path to BCI file.
boost::filesystem::path generateBciFilePath(const boost::filesystem::path& inputDir, common::LaneNumber laneNumber)
{
    std::ostringstream laneDirName;
    laneDirName << boost::format("L%03d") % laneNumber;
   
    std::ostringstream bciFileName;
    bciFileName << boost::format("s_%d.bci") % laneNumber;

    return boost::filesystem::path(inputDir / laneDirName.str() / bciFileName.str());
}


/// \brief BCI file record type definition.
#pragma pack(push, 1)
struct Record
{
public:

    /// \brief Tile number.
    uint32_t tileNumber_;

    /// \brief Number of clusters on the tile.
    uint32_t clustersCount_;

};
#pragma pack(pop)


} // namespace detail


BCIndex::BCIndex(const RawData &bciData)
: tiles_()
{
    this->initTileMetadata(bciData);
}

BCIndex::TileMetadataContainer::const_iterator BCIndex::tileMetadataBegin() const
{
    return tiles_.begin();
}

BCIndex::TileMetadataContainer::const_iterator BCIndex::tileMetadataEnd() const
{
    return tiles_.end();
}

void BCIndex::initTileMetadata(const RawData &bciData)
{
    if (bciData.empty())
    {
        return;
    }

    BCL2FASTQ_ASSERT_MSG(bciData.size() % sizeof(detail::Record) == 0, "Invalid BCI data: data_size=" << bciData.size() << " record_size=" << sizeof(detail::Record));
    BOOST_FOREACH (const detail::Record &record, std::make_pair(reinterpret_cast<const detail::Record *>(&bciData[0]), reinterpret_cast<const detail::Record *>(&bciData[0] + bciData.size())))
    {
        tiles_.push_back(TileMetadata());
        TileMetadata &tileMetadata = tiles_.back();
        tileMetadata.tileNumber_ = record.tileNumber_;
        tileMetadata.clustersCount_ = static_cast<common::ClustersCount>(record.clustersCount_);
    }
}

BCIndex::RawData parseBciFile(const boost::filesystem::path& bciFile)
{
    int errnum;

    boost::filesystem::ifstream bciStream(bciFile, std::ios_base::in | std::ios_base::binary);
    errnum = errno;
    if (!bciStream)
    {
        BOOST_THROW_EXCEPTION(common::IoError(errnum, (boost::format("Unable to open '%s' file for reading") % bciFile.string()).str()));
    }

    BCIndex::RawData bciBuffer;
    errnum = 0;
    std::copy(
        std::istreambuf_iterator<BCIndex::RawData::value_type>(bciStream),
        std::istreambuf_iterator<BCIndex::RawData::value_type>(),
        std::back_inserter(bciBuffer)
    );

    return bciBuffer;
}

bool checkBciFile(boost::filesystem::path inputDir, common::LaneNumber laneNumber)
{
    return boost::filesystem::exists(detail::generateBciFilePath(inputDir, laneNumber));
}

BCIndex createBCIndex(boost::filesystem::path inputDir, common::LaneNumber laneNumber)
{
    return BCIndex(parseBciFile(detail::generateBciFilePath(inputDir, laneNumber)));
}


} // namespace layout


} // namespace bcl2fastq


