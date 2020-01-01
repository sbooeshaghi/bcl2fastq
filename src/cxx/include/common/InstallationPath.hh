/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file InstallatinPath.hh
 *
 * \brief Definition of InstallationPath.
 *
 * \author Aaron Day
 */

#ifndef BCL2FASTQ_COMMON_INSTALLATION_PATH_HH
#define BCL2FASTQ_COMMON_INSTALLATION_PATH_HH

#include <boost/filesystem/path.hpp>

namespace bcl2fastq
{
namespace common
{

class InstallationPath
{
private:
    InstallationPath();

    boost::filesystem::path installationRoot_;

public:
    static InstallationPath& getSingleton();

    boost::filesystem::path expandPath(const boost::filesystem::path& path) const;
};

}
}

#endif // BCL2FASTQ_COMMON_INSTALLATION_PATH_HH
