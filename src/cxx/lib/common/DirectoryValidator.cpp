/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file DirectoryValidator.cpp
 *
 * \brief Implementation of directory validator.
 *
 * \author Aaron Day
 */


#include "common/DirectoryValidator.hh"

#include <boost/foreach.hpp>
#include <boost/filesystem/operations.hpp>

namespace bcl2fastq
{
namespace common
{

DirectoryValidationError::DirectoryValidationError(const boost::filesystem::path& uniquePath,
                                                   const std::string&             uniquePathDescription,
                                                   const boost::filesystem::path& invalidPath,
                                                   const std::string&             invalidPathDescription)
: common::RuntimeError(0, createErrorMessage(uniquePath,
                                             uniquePathDescription,
                                             invalidPath,
                                             invalidPathDescription))
{
}

DirectoryValidationError::DirectoryValidationError(const DirectoryValidationError& other)
: common::RuntimeError(other)
{
}

std::string DirectoryValidationError::createErrorMessage(const boost::filesystem::path& uniquePath,
                                                         const std::string&             uniquePathDescription,
                                                         const boost::filesystem::path& invalidPath,
                                                         const std::string&             invalidPathDescription)
{
    return invalidPathDescription + " path: '" + invalidPath.string() +
           "' clashes with " + uniquePathDescription + " path: '" + uniquePath.string() + "'";
}


DirectoryValidator::DirectoryValidator() : uniquePaths_()
{
}

DirectoryValidator& DirectoryValidator::getSingleton()
{
    static DirectoryValidator instance;
    return instance;
}

void DirectoryValidator::addUniquePath(const boost::filesystem::path& path,
                                       const std::string& description)
{
    uniquePaths_.push_back(std::make_pair(boost::filesystem::absolute(path), description));
}

void DirectoryValidator::validateNoClash(const boost::filesystem::path& path,
                                         const std::string&             description) const
{
    static const auto pathSeparator = boost::filesystem::path("/").make_preferred().native();

    std::string absolutePath = boost::filesystem::absolute(path / pathSeparator).string();

    typedef std::pair<boost::filesystem::path, std::string> PathDescriptionPair;
    BOOST_FOREACH(const PathDescriptionPair& uniquePath,
                  std::make_pair(uniquePaths_.begin(), uniquePaths_.end()))
    {
        if (absolutePath.size() >= uniquePath.first.string().size() &&
            absolutePath.substr(0, uniquePath.first.string().size()) == uniquePath.first.string())
        {
            BOOST_THROW_EXCEPTION(DirectoryValidationError(uniquePath.first,
                                                           uniquePath.second,
                                                           absolutePath,
                                                           description));
        }
    }
}

} // namespace common
} // namespace bcl2fastq


