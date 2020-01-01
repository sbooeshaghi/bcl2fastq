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
## file cppunit.cmake
##
## Configuration file for the cppunit subfolders
##
## author Come Raczy
##
################################################################################

##
## the location of the cppunit shared libraries will be needed for
## LD_LIBRARY_PATH
##

if (HAVE_CPPUNIT)
    get_filename_component(CPPUNIT_LOCATION ${HAVE_CPPUNIT} PATH)
    include_directories(${CPPUNIT_INCLUDE_DIR})
else (HAVE_CPPUNIT)
    message(FATAL_ERROR "cppunit not found")
endif (HAVE_CPPUNIT)

##
## find all the source files
##

file (GLOB BCL2FASTQ_TEST_SOURCE_LIST "*.cpp")

##
## create the targets to build the tests
##

set(BCL2FASTQ_CPPUNIT_TEST_NAME cppunitTest)
add_executable(${BCL2FASTQ_UNIQUE_PREFIX}${BCL2FASTQ_CPPUNIT_TEST_NAME} ${BCL2FASTQ_TEST_SOURCE_LIST})
set_target_properties(${BCL2FASTQ_UNIQUE_PREFIX}${BCL2FASTQ_CPPUNIT_TEST_NAME} PROPERTIES OUTPUT_NAME ${BCL2FASTQ_CPPUNIT_TEST_NAME})
add_test(cppunit_${BCL2FASTQ_UNIQUE_PREFIX}${BCL2FASTQ_CPPUNIT_TEST_NAME} "${CMAKE_CURRENT_BINARY_DIR}/${BCL2FASTQ_CPPUNIT_TEST_NAME}")

include_directories   (${BCL2FASTQ_CXX_ALL_INCLUDES} ${BCL2FASTQ_OPT_INC} "${CMAKE_SOURCE_DIR}/cxx/unittest")

set(BCL2FASTQ_LINK_LIBRARIES "-lpthread")
if    (HAVE_ZLIB)
    set(BCL2FASTQ_LINK_LIBRARIES "${BCL2FASTQ_LINK_LIBRARIES} -lz")
endif (HAVE_ZLIB)
#if    (HAVE_BZLIB)
#    set(BCL2FASTQ_LINK_LIBRARIES "${BCL2FASTQ_LINK_LIBRARIES} -lbz2")
#endif (HAVE_BZLIB)
if    (NOT BCL2FASTQ_FORCE_STATIC_LINK)
    set(BCL2FASTQ_LINK_LIBRARIES "${BCL2FASTQ_LINK_LIBRARIES} -ldl")
endif (NOT BCL2FASTQ_FORCE_STATIC_LINK)

target_link_libraries (${BCL2FASTQ_UNIQUE_PREFIX}${BCL2FASTQ_CPPUNIT_TEST_NAME}
                       bcl2fastq_${BCL2FASTQ_LIB_DIR} ${BCL2FASTQ_AVAILABLE_LIBRARIES} 
                       bcl2fastq_cppunit ${Boost_LIBRARIES}
                       ${BCL2FASTQ_LINK_LIBRARIES} ${BCL2FASTQ_DEP_LIB} ${CPPUNIT_LIBRARY})

##
## Run some sanity check on the source file
##
foreach(BCL2FASTQ_CPPUNIT_SOURCE_FILE ${BCL2FASTQ_TEST_SOURCE_LIST})
    get_filename_component(FILE_NAME ${BCL2FASTQ_CPPUNIT_SOURCE_FILE} NAME)
    set(BCL2FASTQ_CPPUNIT_BINARY_FILE "${CMAKE_CURRENT_BINARY_DIR}/${FILE_NAME}")
    add_custom_command(TARGET ${BCL2FASTQ_UNIQUE_PREFIX}${BCL2FASTQ_CPPUNIT_TEST_NAME} 
                       PRE_BUILD
                       COMMAND ${CMAKE_SOURCE_DIR}/cxx/unittest/check-source.sh ARGS ${BCL2FASTQ_CPPUNIT_BINARY_FILE}.checked ${BCL2FASTQ_CPPUNIT_SOURCE_FILE}
                       COMMENT "Sanity check on ${BCL2FASTQ_CPPUNIT_SOURCE_FILE}")
endforeach(BCL2FASTQ_CPPUNIT_SOURCE_FILE)

add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/RegistryNames.txt 
                   COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/RegistryNames.txt ${CMAKE_CURRENT_BINARY_DIR}
		   DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/RegistryNames.txt
                   COMMENT "Copying RegistryNames.txt for ${BCL2FASTQ_UNIQUE_PREFIX}${BCL2FASTQ_CPPUNIT_TEST_NAME}")

##
## create the targets to run the tests
##
add_custom_command(OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/${BCL2FASTQ_CPPUNIT_TEST_NAME}.passed 
		   COMMAND export LD_LIBRARY_PATH=$$LD_LIBRARY_PATH:${CPPUNIT_LOCATION} && ./${BCL2FASTQ_CPPUNIT_TEST_NAME}
	           COMMAND touch ${BCL2FASTQ_CPPUNIT_TEST_NAME}.passed  
		   DEPENDS ${BCL2FASTQ_UNIQUE_PREFIX}${BCL2FASTQ_CPPUNIT_TEST_NAME} ${CMAKE_CURRENT_BINARY_DIR}/RegistryNames.txt
		   COMMENT "Running unit tests ${BCL2FASTQ_UNIQUE_PREFIX}${BCL2FASTQ_CPPUNIT_TEST_NAME}")
add_custom_target(${BCL2FASTQ_UNIQUE_PREFIX}passed ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${BCL2FASTQ_CPPUNIT_TEST_NAME}.passed)

##
## Copy the data directory from the source tree if available
##

find_path(${CMAKE_CURRENT_SOURCE_DIR}_DATA_DIR data PATHS ${CMAKE_CURRENT_SOURCE_DIR} NO_DEFAULT_PATH)
if (${CMAKE_CURRENT_SOURCE_DIR}_DATA_DIR)
message (STATUS "Adding the data subdirectory for the cppunits under ${BCL2FASTQ_LIB_DIR}")
    add_subdirectory (data)
    add_dependencies(${BCL2FASTQ_UNIQUE_PREFIX}passed ${BCL2FASTQ_UNIQUE_PREFIX}data)
endif (${CMAKE_CURRENT_SOURCE_DIR}_DATA_DIR)
