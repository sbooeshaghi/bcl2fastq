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
## Configuration file for the cxx subfolder
##
## author Come Raczy
##
################################################################################


include (${BCL2FASTQ_GLOBALS_CMAKE})

set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set (BCL2FASTQ_CXX_EXECUTABLE_CMAKE "${CMAKE_SOURCE_DIR}/cmake/cxxExecutable.cmake")
set (BCL2FASTQ_CXX_LIBRARY_CMAKE "${CMAKE_SOURCE_DIR}/cmake/cxxLibrary.cmake")
set (BCL2FASTQ_CXX_CONFIGURE_CMAKE "${CMAKE_SOURCE_DIR}/cmake/cxxConfigure.cmake")

set(BCL2FASTQ_CXX_CONFIG_H_DIR ${CMAKE_CURRENT_BINARY_DIR}/common)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")
include ("${BCL2FASTQ_CXX_CONFIGURE_CMAKE}")

if (HAVE_CPPUNIT AND BCL2FASTQ_UNIT_TESTS)
    set (BCL2FASTQ_CPPUNIT_CMAKE "${CMAKE_SOURCE_DIR}/cmake/cppunit.cmake")
    add_subdirectory (unittest)
endif (HAVE_CPPUNIT AND BCL2FASTQ_UNIT_TESTS)


##
## The include directories
##
set (BCL2FASTQ_CXX_ALL_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/include" "${CMAKE_CURRENT_BINARY_DIR}")

##
## Build all the libraries for the project
##
add_subdirectory (lib)

##
## build all the applications for the project
##

add_subdirectory (bin)
add_subdirectory (libexec)

##
## build all the internal applications for the project
##

#add_subdirectory (libexec)

##
## build the documentation when available
##
include  (FindDoxygen)
message (STATUS "Doxygen: ${DOXYGEN_EXECUTABLE}. Dot: ${DOXYGEN_DOT_EXECUTABLE}.")
if (DOXYGEN_FOUND)
    set (DOXYFILE ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
    message (STATUS "Creating Doxygen config file: ${DOXYFILE}")
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${DOXYFILE} @ONLY IMMEDIATE)
    add_custom_target(doc ${DOXYGEN_EXECUTABLE} ${DOXYFILE})
endif (DOXYGEN_FOUND)


