/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file RegistryName.hh
 *
 * \brief Management of the registry names for the cppunit tests.
 *
 * \author Come Raczy
 */


#ifndef BCL2FASTQ_UNIT_TEST_REGISTRY_NAME
#define BCL2FASTQ_UNIT_TEST_REGISTRY_NAME

#include <stdexcept>
#include <string>
#include <vector>

const std::vector<std::string> &getRegistryNameList();
std::string registryName(const std::string &name) throw (std::invalid_argument);

#endif // #ifndef BCL2FASTQ_UNIT_TEST_REGISTRY_NAME


