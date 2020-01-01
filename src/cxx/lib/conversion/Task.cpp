/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file Task.cpp
 *
 * \brief Implementation of a workload item a.k.a. task.
 *
 * \author Marek Balint
 */


#include "conversion/Task.hh"
#include "conversion/Stage.hh"

#include "common/Debug.hh"

namespace bcl2fastq {
namespace conversion {

Task::Task(std::shared_ptr<TaskManager> taskManager)
 : taskManager_(taskManager)
{
}

Task::~Task()
{
}

//
// Class TaskManager::GroupSerializer
//

TaskManager::GroupSerializer::GroupSerializer()
: nextGroupNumber_(0)
, cvGroupNumber_()
, groupMutex_()
{
}

void TaskManager::GroupSerializer::waitForNextGroupNumber(uint32_t groupNumber)
{
    std::unique_lock<std::mutex> lock(groupMutex_);
    cvGroupNumber_.wait(lock, [this, groupNumber] { return nextGroupNumber_ == groupNumber; });
}

void TaskManager::GroupSerializer::signalProceed()
{
    {
        std::unique_lock<std::mutex> lock(groupMutex_);
        ++nextGroupNumber_;
    }

    cvGroupNumber_.notify_all();
}

uint32_t TaskManager::GroupSerializer::GetNextGroupNumber() const
{
    return nextGroupNumber_;
}

} // namespace task
} // namespace bcl2fastq


