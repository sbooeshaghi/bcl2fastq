/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file DirectoryValidator.hh
 *
 * \brief Declaration of directory validator.
 *
 * \author Aaron Day
 */


#ifndef BCL2FASTQ_COMMON_DIRECTORY_VALIDATOR_HH
#define BCL2FASTQ_COMMON_DIRECTORY_VALIDATOR_HH

#include "common/Exceptions.hh"

#include <boost/noncopyable.hpp>
#include <boost/filesystem/path.hpp>
#include <vector>
#include <utility>

namespace bcl2fastq
{
namespace common
{

class DirectoryValidationError : public common::RuntimeError
{
public:
    DirectoryValidationError(const boost::filesystem::path& uniquePath,
                             const std::string&             uniquePathDescription,
                             const boost::filesystem::path& invalidPath,
                             const std::string&             invalidPathDescription);

    DirectoryValidationError(const DirectoryValidationError& other);

    std::string createErrorMessage(const boost::filesystem::path& uniquePath,
                                   const std::string&             uniquePathDescription,
                                   const boost::filesystem::path& invalidPath,
                                   const std::string&             invalidPathDescription);
};

/// \brief Directory validator.
class DirectoryValidator : private boost::noncopyable
{
public:
    static DirectoryValidator& getSingleton();

    /// \brief Add a unique path.
    /// \param Unique path
    /// \param description of path (used in error message)
    void addUniquePath(const boost::filesystem::path& path,
                       const std::string&             description);

    /// \brief Validate that there is no clash with the unique paths.
    /// \param path to validate
    /// \param description of path (used in error message)
    void validateNoClash(const boost::filesystem::path& path,
                         const std::string&             description) const;

private:

    /// \brief Constructor.
    DirectoryValidator();

    /// \brief uniquePaths_.
    ///
    /// Source: All paths which must be unique.
    std::vector< std::pair<boost::filesystem::path, std::string> > uniquePaths_;
};


} // namespace common
} // namespace bcl2fastq


#endif // BCL2FASTQ_COMMON_DIRECTORY_VALIDATOR_HH


