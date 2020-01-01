################################
# Unit Tests
################################
#

include (${BCL2FASTQ_GLOBALS_CMAKE})

set (GTEST_MODULE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../modules/googletest/googletest")
set (BCL2FASTQ_INCLUDE "${CMAKE_CURRENT_SOURCE_DIR}/../cxx/include")

##
## Collect source files
## 
file(GLOB_RECURSE GTEST_LIBRARY_SOURCES *.cpp *.c)
foreach (SOURCE_FILE ${GTEST_LIBRARY_SOURCES})
    set(BCL2FASTQ_GTEST_SOURCES ${BCL2FASTQ_GTEST_SOURCES} ${SOURCE_FILE})
endforeach (SOURCE_FILE)
message(STATUS "Sources ${BCL2FASTQ_GTEST_SOURCES}")

file(COPY "testData" DESTINATION "${CMAKE_BINARY_DIR}")

##
## The include directories
##

include_directories(${GTEST_MODULE_DIR}/include ${GTEST_MODULE_DIR})
include_directories(${BCL2FASTQ_INCLUDE}) 

# Add test cpp file
add_executable( unitTest ${BCL2FASTQ_GTEST_SOURCES} ${GTEST_MODULE_DIR}/src/gtest-all.cc ${GTEST_MODULE_DIR}/src/gtest_main.cc)

add_dependencies(unitTest 
    bcl2fastq_common 
    bcl2fastq_config 
    bcl2fastq_io 
    bcl2fastq_data 
    bcl2fastq_layout 
    bcl2fastq_stats 
    bcl2fastq_conversion)

# Link test executable against gtest & gtest_main
target_link_libraries(unitTest 
    bcl2fastq_conversion
    bcl2fastq_stats 
    bcl2fastq_data 
    bcl2fastq_layout 
    bcl2fastq_data 
    bcl2fastq_io
    bcl2fastq_config
    bcl2fastq_common
    rt)

find_package( Boost 1.54.0 REQUIRED 
    thread
    chrono
    date_time
    filesystem
    program_options
    regex
    serialization
    system
    timer
    iostreams)
if (Boost_FOUND)
    target_link_libraries( unitTest ${Boost_LIBRARIES} )
endif(Boost_FOUND)

find_package( ZLIB REQUIRED )
if ( ZLIB_FOUND )
    include_directories( ${ZLIB_INCLUDE_DIRS} )
    target_link_libraries( unitTest ${ZLIB_LIBRARIES} )
endif( ZLIB_FOUND )

add_test( unitTest unitTest )
