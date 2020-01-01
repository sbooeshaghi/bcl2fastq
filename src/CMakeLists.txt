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
## Top level cmake configuration file.
##
## author Come Raczy
##
################################################################################


cmake_minimum_required(VERSION 2.8.4)
enable_testing()

if (DEFINED CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel.")
else()
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel.")
endif()

project (bcl2fastq2)
set(CPACK_RPM_SPEC_MORE_DEFINE "Prefix: ${BCL2FASTQ_PREFIX}")

# Version information
# these are just some defaults, real values should come from configure script via command line arguments
if (NOT DEFINED BCL2FASTQ_NAME_SHORT)
    set (BCL2FASTQ_NAME_SHORT "${PROJECT_NAME}")
endif()
if (NOT DEFINED BCL2FASTQ_NAME_LONG)
    set (BCL2FASTQ_NAME_LONG "${PROJECT_NAME}")
endif()
if (NOT DEFINED BCL2FASTQ_COPYRIGHT)
    set (BCL2FASTQ_COPYRIGHT "Copyright (c) Illumina, Inc.")
endif()
if (NOT DEFINED BCL2FASTQ_VERSION_MAJOR)
    set (BCL2FASTQ_VERSION_MAJOR "2")
endif()
if (NOT DEFINED BCL2FASTQ_VERSION_MINOR)
    set (BCL2FASTQ_VERSION_MINOR "0")
endif()
if (NOT DEFINED BCL2FASTQ_VERSION_PATCH)
    set (BCL2FASTQ_VERSION_PATCH "0")
endif()
if (NOT DEFINED BCL2FASTQ_VERSION_BUILD)
    set (BCL2FASTQ_VERSION_BUILD "0")
endif()

if (CPACK_GENERATOR)
    message (STATUS "Configuring to produce the following package types: ${CPACK_GENERATOR}")
    SET(CPACK_PACKAGE_VENDOR "Illumina")
    SET(CPACK_PACKAGE_VERSION_MAJOR "${BCL2FASTQ_VERSION_MAJOR}")
    SET(CPACK_PACKAGE_VERSION_MINOR "${BCL2FASTQ_VERSION_MINOR}")
    SET(CPACK_PACKAGE_VERSION_PATCH "${BCL2FASTQ_VERSION_PATCH}.${BCL2FASTQ_VERSION_BUILD}")
    SET(CPACK_PACKAGE_VERSION "v${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
    SET(CPACK_SET_DESTDIR ON)
    SET(CPACK_PACKAGE_RELOCATABLE OFF)
    SET(CPACK_DEBIAN_PACKAGE_OWNER root)
    SET(CPACK_DEBIAN_PACKAGE_GROUP root)
    INCLUDE(CPack)
endif (CPACK_GENERATOR)

if (NOT CMAKE_PARALLEL)
    set (CMAKE_PARALLEL "1")
endif (NOT CMAKE_PARALLEL)

# Installation directories
if    (NOT BCL2FASTQ_PREFIX)
    set(BCL2FASTQ_PREFIX "/usr/local")
endif (NOT BCL2FASTQ_PREFIX)

if    (NOT CMAKE_INSTALL_PREFIX)
    set(CMAKE_INSTALL_PREFIX "${BCL2FASTQ_PREFIX}")
endif (NOT CMAKE_INSTALL_PREFIX)

if (CMAKE_INSTALL_PREFIX)
    string(SUBSTRING "${CMAKE_INSTALL_PREFIX}" 0 1  slash)
    string(COMPARE NOTEQUAL "/" "${slash}" isslash)
    if    (isslash)
        install(CODE "
            string(LENGTH \"\$ENV{DESTDIR}\" ddlen)
            if (ddlen)
                message (FATAL_ERROR \"Installation prefix must begin with '/' if DESTDIR is set.\")
            endif (ddlen)
        ")
    endif (isslash)
endif (CMAKE_INSTALL_PREFIX)

if    (NOT BCL2FASTQ_EXEC_PREFIX)
    set(BCL2FASTQ_EXEC_PREFIX "." CACHE PATH "Install bin directory" FORCE)
endif (NOT BCL2FASTQ_EXEC_PREFIX)

if    (NOT BCL2FASTQ_BINDIR)
    set(BCL2FASTQ_BINDIR "${BCL2FASTQ_EXEC_PREFIX}/bin" 
        CACHE PATH "Install bin directory" FORCE)
endif (NOT BCL2FASTQ_BINDIR)

if    (NOT BCL2FASTQ_LIBDIR)
    set(BCL2FASTQ_LIBDIR
    "${BCL2FASTQ_EXEC_PREFIX}/lib/${BCL2FASTQ_NAME_SHORT}-${BCL2FASTQ_VERSION}"
        CACHE PATH "Install lib directory" FORCE)
endif (NOT BCL2FASTQ_LIBDIR)

if    (NOT BCL2FASTQ_LIBEXECDIR)
    set(BCL2FASTQ_LIBEXECDIR "${BCL2FASTQ_EXEC_PREFIX}/libexec/${BCL2FASTQ_NAME_SHORT}-${BCL2FASTQ_VERSION}"
        CACHE PATH "Install libexec directory" FORCE)
endif (NOT BCL2FASTQ_LIBEXECDIR)

if    (NOT BCL2FASTQ_INCLUDEDIR)
    set(BCL2FASTQ_INCLUDEDIR "include/${BCL2FASTQ_NAME_SHORT}-${BCL2FASTQ_VERSION}"
        CACHE PATH "Install include directory" FORCE)
endif (NOT BCL2FASTQ_INCLUDEDIR)

if    (NOT BCL2FASTQ_DATADIR)
    set(BCL2FASTQ_DATADIR "share"
        CACHE PATH "Install data directory" FORCE)
endif (NOT BCL2FASTQ_DATADIR)

if    (NOT BCL2FASTQ_ETCDIR)
    set(BCL2FASTQ_ETCDIR
    "${BCL2FASTQ_EXEC_PREFIX}/etc/${BCL2FASTQ_NAME_SHORT}-${BCL2FASTQ_VERSION}"
        CACHE PATH "Install sysconfdir directory" FORCE)
endif (NOT BCL2FASTQ_ETCDIR)

if    (NOT BCL2FASTQ_DOCDIR)
    set(BCL2FASTQ_DOCDIR "doc/${BCL2FASTQ_NAME_SHORT}-${BCL2FASTQ_VERSION}"
        CACHE PATH "Install doc directory" FORCE)
endif (NOT BCL2FASTQ_DOCDIR)

if    (NOT BCL2FASTQ_MANDIR)
    set(BCL2FASTQ_MANDIR "man" 
        CACHE PATH "Install man directory" FORCE)
endif (NOT BCL2FASTQ_MANDIR)

message (STATUS "install      prefix: ${BCL2FASTQ_PREFIX}")
message (STATUS "install exec prefix: ${BCL2FASTQ_EXEC_PREFIX}")
message (STATUS "install     bin dir: ${BCL2FASTQ_BINDIR}")
message (STATUS "install     lib dir: ${BCL2FASTQ_LIBDIR}")
message (STATUS "install libexec dir: ${BCL2FASTQ_LIBEXECDIR}")
message (STATUS "install include dir: ${BCL2FASTQ_INCLUDEDIR}")
message (STATUS "install    data dir: ${BCL2FASTQ_DATADIR}")
message (STATUS "install     doc dir: ${BCL2FASTQ_DOCDIR}")
message (STATUS "install     man dir: ${BCL2FASTQ_MANDIR}")

add_custom_target(BCL2FASTQ_CXX)
add_custom_target(BCL2FASTQ_OPT)

# Get around broken BOOST cmake macros if BOOST_ROOT set
set(CMAKEBOOSTROOT $ENV{BOOST_ROOT})
if( NOT ${CMAKEBOOSTROOT} STREQUAL "" )
   set(Boost_NO_SYSTEM_PATHS ON)
   message("BOOST_ROOT is set, Boost_NO_SYSTEM_PATHS is ${Boost_NO_SYSTEM_PATHS}")
endif( NOT ${CMAKEBOOSTROOT} STREQUAL "" )

set(BOOST_REDIST_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../redist")

# required boost libraries
set (BCL2FASTQ_BOOST_VERSION 1.54.0)
set (BCL2FASTQ_BOOST_COMPONENTS
    chrono
    date_time
    filesystem
    iostreams
    program_options
    regex
    serialization
    system
    timer
    thread
)
set (Boost_USE_MULTITHREAD ON)

# required libxml2 and libxslt libraries
set (BCL2FASTQ_LIBXML2_VERSION 2.7.8)
set (LIBXML2_REDIST_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../redist")
set (LIBXSLT_REDIST_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../redist")

set (LIBXML2_INSTAL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../redist")
set (LIBXSLT_INSTAL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../redist")

set (BCL2FASTQ_LIBXSLT_VERSION 1.1.26)
set (LIBXSLT_REDIST_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../redist")

# test
set (GOOGLE_TEST_MODULE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../modules/googletest/googletest")

# globals and macros
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")
set (BCL2FASTQ_GLOBALS_CMAKE "${CMAKE_SOURCE_DIR}/cmake/globals.cmake")
set (BCL2FASTQ_MACROS_CMAKE "${CMAKE_SOURCE_DIR}/cmake/macros.cmake")

# redist includes
include ("${CMAKE_SOURCE_DIR}/cmake/boost.cmake")
include ("${CMAKE_SOURCE_DIR}/cmake/bcl2fastq_redist_macros.cmake")

# subdirs
add_subdirectory (cmake/preInstall)
add_subdirectory (css)
add_subdirectory (cxx)
if (EXISTS ${GOOGLE_TEST_MODULE_DIR})
    add_subdirectory (test)
endif()
add_subdirectory (xsl)
add_subdirectory (cmake/postInstall)

