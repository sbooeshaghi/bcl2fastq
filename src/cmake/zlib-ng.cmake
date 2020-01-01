################################################################################
##
## BCL to FASTQ file converter
## Copyright (c) 2007-2016 Illumina, Inc.
##
## This software is covered by the accompanying EULA
## and certain third party copyright/licenses, and any user of this
## source file is bound by the terms therein.
##
################################################################################
##
## file CMakeLists.txt
##
## Configuration file for zlib-ng installation
##
## author Aaron Day
##
################################################################################

set (ENV{CFLAGS} "-mavx -O3 -I${CMAKE_BINARY_DIR}/bootstrap/include")
set(REINSTDIR ${CMAKE_BINARY_DIR}/bootstrap)

set(CMAKE_C_FLAGS_CPY "${CMAKE_C_FLAGS}")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mavx -O3")
redist_package(ZLIB-NG ${BCL2FASTQ_ZLIB-NG_VERSION}
               "--zlib-compat;--prefix=${REINSTDIR}")
set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS_CPY})

find_library_redist(ZLIB-NG ${REINSTDIR} zlib.h z)
include_directories(BEFORE SYSTEM ${ZLIB_INCLUDE_DIR})
