/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file InstallationPath.cpp
 *
 * file Locates the installation path.
 *
 * \author Aaron Day
 */

#ifdef WIN32
#include <windows.h>
#endif

#include "common/InstallationPath.hh"

#include "common/Debug.hh"

namespace bcl2fastq
{
namespace common
{

InstallationPath& InstallationPath::getSingleton()
{
    static InstallationPath instance;
    return instance;
}


InstallationPath::InstallationPath()
    : installationRoot_()
{
#ifndef WIN32
    char szBuffer[10240];
    BCL2FASTQ_ASSERT_MSG(-1 != readlink("/proc/self/exe", szBuffer, sizeof(szBuffer)), "TODO: handle the readlink error: " << errno);

    installationRoot_ = boost::filesystem::path(szBuffer).parent_path().parent_path();
#else
    wchar_t buffer[MAX_PATH];
    HMODULE hModule = GetModuleHandle(NULL);
    if (hModule != NULL)
    {
        // When passing NULL to GetModuleHandle, it returns handle of exe itself
        GetModuleFileName(hModule,buffer, (sizeof(buffer)));
    }

    installationRoot_ = boost::filesystem::path(buffer).parent_path().parent_path();
#endif
}

boost::filesystem::path InstallationPath::expandPath(const boost::filesystem::path& path) const
{
    if (installationRoot_.empty())
    {
        return path;
    }

    return installationRoot_ / path;
}

} // namespace common
} // namespace bcl2fastq

