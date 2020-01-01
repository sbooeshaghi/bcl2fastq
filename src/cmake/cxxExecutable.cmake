################################################################################
##
## BCL to FASTQ file converter
## Copyright (c) 2007-2015 Illumina, Inc.
##
## This software is covered by the accompanying EULA
## and certain third party copyright/licenses, and any user of this
## source file is bound by the terms therein.
##
################################################################################
##
## file cxxExecutable.cmake
##
## CMake configuration file for all the c++ executables.
##
## author Roman Petrovski
##
################################################################################


include (${BCL2FASTQ_GLOBALS_CMAKE})

# Support for static linking. Notice this is done here and not in cxxConfigure to 
# allow dynamic linking for cppunit tests as some platforms lack static libcppunit.
# Note that this implies that all libraries must be found with the
# exact file name (libXXX.a or libXXX.so)
if    (BCL2FASTQ_FORCE_STATIC_LINK)
    message(STATUS "All libraries will be statically linked")
    set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "-static")
    set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "-static")
    # ensure that even if cmake decides to allow for dynamic libs resolution,
    # this gets overriden into static...
    set(CMAKE_EXE_LINK_DYNAMIC_CXX_FLAGS ${CMAKE_EXE_LINK_STATIC_CXX_FLAGS})
    set(BCL2FASTQ_LIBRARY_PREFIX ${CMAKE_STATIC_LIBRARY_PREFIX})
    set(BCL2FASTQ_LIBRARY_SUFFIX ${CMAKE_STATIC_LIBRARY_SUFFIX})
else  (BCL2FASTQ_FORCE_STATIC_LINK)
    set(BCL2FASTQ_LIBRARY_PREFIX "")
    set(BCL2FASTQ_LIBRARY_SUFFIX "")
endif (BCL2FASTQ_FORCE_STATIC_LINK)

get_filename_component(BCL2FASTQ_CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
message (STATUS "Adding the cpp    program subdirectory: ${BCL2FASTQ_CURRENT_DIR_NAME}")
include_directories (${BCL2FASTQ_CXX_ALL_INCLUDES})
include_directories (${CMAKE_CURRENT_BINARY_DIR})
include_directories (${CMAKE_CURRENT_SOURCE_DIR})
include_directories (${BCL2FASTQ_CXX_CONFIG_H_DIR})


