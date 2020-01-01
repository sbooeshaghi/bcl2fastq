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
## file cxxLibrary.cmake
##
## CMake configuration file for all the c++ libraries
##
## author Come Raczy
##
################################################################################


include_directories (${BCL2FASTQ_CXX_ALL_INCLUDES})
include_directories (${CMAKE_CURRENT_BINARY_DIR})
include_directories (${CMAKE_CURRENT_SOURCE_DIR})
include_directories (${BCL2FASTQ_CXX_CONFIG_H_DIR})

get_filename_component(BCL2FASTQ_CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
message (STATUS "Adding the cpp    library subdirectory: ${BCL2FASTQ_CURRENT_DIR_NAME}")

##
## Some generators (VS) require all targets to be unique across the project.
## Therefore, a unique prefix is needed to create the target names which are
## shared across libraries
##

string(REGEX REPLACE ${CMAKE_SOURCE_DIR}/cpp/ "" TMP1 ${CMAKE_CURRENT_SOURCE_DIR}/)
string(REGEX REPLACE "/" "_" BCL2FASTQ_UNIQUE_PREFIX ${TMP1})

##
## build the library
##
file(GLOB_RECURSE BCL2FASTQ_LIBRARY_SOURCES_WITH_CPPUNIT *.cpp *.c)

foreach (SOURCE_FILE ${BCL2FASTQ_LIBRARY_SOURCES_WITH_CPPUNIT})
    string(REGEX MATCH "cppunit" CPPUNIT_MATCH ${SOURCE_FILE} )
    if (NOT CPPUNIT_MATCH)
        set(BCL2FASTQ_LIBRARY_SOURCES ${BCL2FASTQ_LIBRARY_SOURCES} ${SOURCE_FILE})
    endif (NOT CPPUNIT_MATCH)
endforeach (SOURCE_FILE)

foreach (SOURCE_FILE ${BCL2FASTQ_LIBRARY_SOURCES})
    get_filename_component(SOURCE_NAME ${SOURCE_FILE} NAME_WE)
    if (${SOURCE_NAME}_COMPILE_FLAGS)
        set_source_files_properties(${SOURCE_FILE} PROPERTIES COMPILE_FLAGS ${${SOURCE_NAME}_COMPILE_FLAGS})
    endif (${SOURCE_NAME}_COMPILE_FLAGS)
endforeach (SOURCE_FILE)

#include_directories (${BCL2FASTQ_COMMON_INCLUDE} )
add_library         (bcl2fastq_${BCL2FASTQ_LIB_DIR} STATIC ${BCL2FASTQ_LIBRARY_SOURCES})
add_dependencies(bcl2fastq_${BCL2FASTQ_LIB_DIR} BCL2FASTQ_OPT)

##
## build the unit tests if any (this should be mandatory really)
##

if (HAVE_CPPUNIT AND BCL2FASTQ_UNIT_TESTS)
    find_path(${CMAKE_CURRENT_SOURCE_DIR}_CPPUNIT_DIR cppunit PATHS ${CMAKE_CURRENT_SOURCE_DIR} NO_DEFAULT_PATH)
    if (${CMAKE_CURRENT_SOURCE_DIR}_CPPUNIT_DIR)
        message (STATUS "Adding the cppunit subdirectory for ${BCL2FASTQ_LIB_DIR}")
        add_subdirectory (cppunit)
    endif (${CMAKE_CURRENT_SOURCE_DIR}_CPPUNIT_DIR)
endif(HAVE_CPPUNIT AND BCL2FASTQ_UNIT_TESTS)


