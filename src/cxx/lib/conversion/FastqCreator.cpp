/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file FastqCreator.cpp
 *
 * \brief Implementation of FASTQ creator.
 *
 * \author Marek Balint
 * \author Mauricio Varea
 */


#include <string>
#include <algorithm>
#include <iterator>
#include <utility>

#include "conversion/FastqCreator.hh"

namespace bcl2fastq
{
namespace conversion
{

std::atomic<uint32_t> FastqCreateTaskManager::numTaskManagers_(0);
std::condition_variable FastqCreateTaskManager::cvAllTaskManagersDone_;
std::mutex FastqCreateTaskManager::mut_;

FastqCreateTaskManager::FastqCreateTaskManager(std::shared_ptr<DemuxBuffer>& inputBuffer,
                                               ThreadSafeQueue<std::shared_ptr<DemuxBuffer>>& inputBuffersToRecycle,
                                               std::shared_ptr<FastqBuffer>& outputBuffer,
                                               ThreadSafeMap<uint32_t, std::shared_ptr<FastqBuffer>>& outputBufferMap)
: TaskManager()
, inputBuffer_(inputBuffer)
, inputBuffersToRecycle_(inputBuffersToRecycle)
, outputBuffer_(outputBuffer)
, outputBufferMap_(outputBufferMap)
{
    ++numTaskManagers_;
}

FastqCreateTaskManager::~FastqCreateTaskManager()
{
    static int sequentialId = 0; ++sequentialId;

    outputBuffer_->setGroupNumber(inputBuffer_->getGroupNumber());

    // We're done with the inputBuffer; Recycle it.
    inputBuffersToRecycle_.addData(inputBuffer_);

    BCL2FASTQ_LOG(common::LogLevel::DEBUG) << "Q fastq-writer: " << outputBuffer_->getGroupNumber() << std::endl;

    // TODO: This can block. Do this on a separate thread. Be careful about ensuring the thread completes.
    outputBufferMap_.addData(outputBuffer_->getGroupNumber(), outputBuffer_);
    --numTaskManagers_;

    if (numTaskManagers_ == 0)
    {
        // This could be called more than once on different threads. I think
        // this is ok.
        cvAllTaskManagersDone_.notify_all();
    }
}

void FastqCreateTaskManager::waitForAllTasksToFinish()
{
    std::unique_lock<std::mutex> lock(mut_);
    cvAllTaskManagersDone_.wait(lock, [] { return FastqCreateTaskManager::numTaskManagers_ == 0; });
}

AllConversionStats::AllConversionStats(ConversionStatsForSampleTileBarcode& summaryConversionStats,
                                       stats::BarcodeHits& summaryUnknownBarcodes)
: summaryConversionStats_(summaryConversionStats)
, summaryUnknownBarcodes_(summaryUnknownBarcodes)
, conversionStats_(summaryConversionStats)
, unknownBarcodes_(summaryUnknownBarcodes)
{
}

AllConversionStats::~AllConversionStats()
{
// TODO: Ensure all of these objects are deleted on a single thread

    size_t sampleIndex = 0;
    for (const auto& statsForSample : conversionStats_)
    {
        size_t tileIndex = 0;
        for (const auto& statsForTile : statsForSample)
        {
            size_t barcodeIndex = 0;
            for (const auto& statsForBarcode : statsForTile)
            {
                summaryConversionStats_.at(sampleIndex).at(tileIndex).at(barcodeIndex).first += statsForBarcode.first;
                size_t readIndex = 0;
                for (const auto& statsForRead : statsForBarcode.second)
                {
                    summaryConversionStats_.at(sampleIndex).at(tileIndex).at(barcodeIndex).second.at(readIndex) += statsForRead;
                    ++readIndex;
                }

                ++barcodeIndex;
            }

            ++tileIndex;
        }

        ++sampleIndex;
    }

    summaryUnknownBarcodes_.merge( unknownBarcodes_ );
}


} // namespace conversion
} // namespace bcl2fastq


