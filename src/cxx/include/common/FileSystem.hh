/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file FileSystem.hh
 *
 * file system helper utilities.
 *
 * \author Roman Petrovski
 * \author Mauricio Varea
 */

#ifndef BCL2FASTQ_COMMON_FILE_SYSTEM_HH
#define BCL2FASTQ_COMMON_FILE_SYSTEM_HH

#include <vector>
#include <boost/filesystem/path.hpp>

#include "common/Debug.hh"

namespace bcl2fastq
{
namespace common
{

void createDirectories(std::vector<boost::filesystem::path> createList);


inline bool isDotGzPath(const boost::filesystem::path& path)
{
    static const char dotGz[] = {'.', 'g', 'z'};
    return path.string().length() > sizeof(dotGz) &&
        0 == path.string().compare(path.string().size() - sizeof(dotGz), sizeof(dotGz), dotGz);
}

char getDirectorySeparatorChar();

} //namespace common
} //namespace bcl2fastq

#endif // BCL2FASTQ_COMMON_FILE_SYSTEM_HH
