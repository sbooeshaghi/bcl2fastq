################################################################################
##
## BCL to FASTQ file converter
## Copyright (c) 2007-2017 Illumina, Inc.
##
## This software is covered by the accompanying EULA
## and certain third party copyright/licenses, and any user of this
## source file is bound by the terms therein.
##
################################################################################
##
## file CMakeLists.txt
##
## Configuration file for the cxx/unittest subfolder
##
## author Come Raczy
##
################################################################################


include_directories(${CPPUNIT_INCLUDE_DIR})
add_library(bcl2fastq_cppunit cppunitTest.cpp RegistryName.cpp)


