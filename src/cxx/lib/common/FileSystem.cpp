/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file FileSystem.cpp
 *
 * file system helper utilities.
 *
 * \author Roman Petrovski
 * \author Mauricio Varea
 */
#include <boost/foreach.hpp>
#include <boost/filesystem.hpp>

#include "common/Exceptions.hh"
#include "common/FileSystem.hh"

namespace bcl2fastq
{
namespace common
{

void createDirectories(std::vector<boost::filesystem::path> createList)
{
    std::sort(createList.begin(), createList.end());
    createList.erase(std::unique(createList.begin(), createList.end()), createList.end());
    BOOST_FOREACH(const boost::filesystem::path &directory, createList)
    {
        if (!boost::filesystem::exists(directory))
        {
            BCL2FASTQ_LOG(common::LogLevel::INFO) << "creating directory " << directory << std::endl;
            boost::system::error_code errorCode;
            if(!boost::filesystem::create_directories(directory, errorCode) &&!exists(directory))
            {
                using common::IoError;
                BOOST_THROW_EXCEPTION(IoError(errorCode.value(), "Failed to create directory " + directory.string()));
            }
        }
    }
}


char makeDirectorySeparatorChar()
{
    boost::filesystem::path slash("/");
    slash.make_preferred();
    BCL2FASTQ_ASSERT_MSG(1 == slash.string().size(), "Unexpected directory separator char length greater than 1: " << slash)
    return *slash.string().begin();
}

static const char DIRECTORY_SEPARATOR_CHAR = makeDirectorySeparatorChar();

char getDirectorySeparatorChar()
{
    return  DIRECTORY_SEPARATOR_CHAR;
}


} // namespace common
} // namespace bcl2fastq
