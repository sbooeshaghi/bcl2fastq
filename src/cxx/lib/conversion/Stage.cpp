/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file Stage.cpp
 *
 * \brief Implementation of a processing stage.
 *
 * \author Marek Balint
 */


//#define WITH_CPU_TIMER


#include <boost/thread/locks.hpp>
#include <boost/format.hpp>

#include <future>
#include "common/Timer.hh"
#include "common/Logger.hh"
#include "conversion/Stage.hh"

#include <chrono>
#include <thread>

namespace bcl2fastq {
namespace conversion {


Stage::Stage(const std::string& stageName,
             TaskQueue&         taskQueue)
: taskQueue_(taskQueue)
, stageName_(stageName)
{
}

Stage::~Stage()
{
}

void Stage::run()
{
    try
    {
        BCL2FASTQ_LOG(common::LogLevel::DEBUG) << "Starting " << stageName_ << " stage" << std::endl;

        while (this->startNewTasks()) { }
    }
    catch (...)
    {
        this->terminate();
        throw;
    }

    BCL2FASTQ_LOG(common::LogLevel::DEBUG) << stageName_ << " stage done" << std::endl;
}

TaskQueue & Stage::getTaskQueue()
{
    return taskQueue_;
}

} // namespace task
} // namespace bcl2fastq

