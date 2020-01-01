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
## file cxxConfigure.cmake
##
## CMake configuration file for c++ executables.
##
## author Come Raczy
##
################################################################################


include ("${BCL2FASTQ_MACROS_CMAKE}")


include(TestBigEndian)
TEST_BIG_ENDIAN(BCL2FASTQ_IS_BIG_ENDIAN)

include(CheckFunctionExists)

bcl2fastq_find_header_or_die(HAVE_INTTYPES_H inttypes.h)
bcl2fastq_find_header_or_die(HAVE_MALLOC_H malloc.h)
bcl2fastq_find_header_or_die(HAVE_MCHECK_H mcheck.h)
bcl2fastq_find_header_or_die(HAVE_MEMORY_H memory.h)
bcl2fastq_find_header_or_die(HAVE_SIGNAL_H signal.h)
bcl2fastq_find_header_or_die(HAVE_STDINT_H stdint.h)
bcl2fastq_find_header_or_die(HAVE_STDLIB_H stdlib.h)
bcl2fastq_find_header_or_die(HAVE_STRING_H string.h)
bcl2fastq_find_header_or_die(HAVE_STRINGS_H strings.h)
bcl2fastq_find_header_or_die(HAVE_TIME_H time.h)
bcl2fastq_find_header_or_die(HAVE_UNISTD_H unistd.h)
bcl2fastq_find_header_or_die(HAVE_SYS_STAT_H  sys/stat.h)
bcl2fastq_find_header_or_die(HAVE_SYS_TYPES_H sys/types.h)

# Math functions that might be missing in some flavors of c++
set (CMAKE_REQUIRED_LIBRARIES m)
check_function_exists(floorf HAVE_FLOORF)
check_function_exists(round  HAVE_ROUND)
check_function_exists(roundf HAVE_ROUNDF)
check_function_exists(powf HAVE_POWF)
check_function_exists(erf HAVE_ERF)
check_function_exists(erf HAVE_ERFF)
check_function_exists(erfc HAVE_ERFC)
check_function_exists(erfc HAVE_ERFCF)

# Systems calls
check_function_exists(sysconf HAVE_SYSCONF)
check_function_exists(clock HAVE_CLOCK)


# optional support for numa
bcl2fastq_find_library(NUMA numa.h numa)
if    (HAVE_NUMA)
    message(STATUS "NUMA supported")
    include_directories(BEFORE SYSTEM ${NUMA_INCLUDE_DIR})
    set(BCL2FASTQ_ADDITIONAL_LIB ${BCL2FASTQ_ADDITIONAL_LIB} "${NUMA_LIBRARY}")
else  (HAVE_NUMA)
    message(STATUS "No support for NUMA")
endif (HAVE_NUMA)

# optional support for gzip compression
bcl2fastq_find_library(ZLIB zlib.h z)
if    (HAVE_ZLIB)
    set  (BCL2FASTQ_ADDITIONAL_LIB ${BCL2FASTQ_ADDITIONAL_LIB} z)
    message(STATUS "gzip compression supported")
else  (HAVE_ZLIB)
    message(FATAL_ERROR "No support for gzip compression")
endif (HAVE_ZLIB)

# optional support for librt
bcl2fastq_find_library(LIBRT time.h rt)
if    (HAVE_LIBRT)
    set  (BCL2FASTQ_ADDITIONAL_LIB ${BCL2FASTQ_ADDITIONAL_LIB} rt)
    message(STATUS "librt supported")
else  (HAVE_LIBRT)
    message(FATAL_ERROR "No support for librt")
endif (HAVE_LIBRT)


bcl2fastq_find_boost(${BCL2FASTQ_BOOST_VERSION} "${BCL2FASTQ_BOOST_COMPONENTS}")

bcl2fastq_find_library(CPGPLOT cpgplot.h cpgplot)
bcl2fastq_find_library(PGPLOT cpgplot.h pgplot)
bcl2fastq_find_library(X11 X.h X11)

set(REINSTDIR ${CMAKE_BINARY_DIR}/bootstrap)

# CPPUNIT
bcl2fastq_find_library(CPPUNIT "cppunit/Test.h" cppunit${CPPUNIT_DEBUG})

# XML2 - bootstrap first (if necessary) so xslt can build against it 
# XSLT and EXSLT
if((NOT HAVE_LIBXML2) OR (NOT HAVE_LIBXSLT))
  find_package_version(LibXml2 ${BCL2FASTQ_LIBXML2_VERSION})
  find_package_version(LibXslt ${BCL2FASTQ_LIBXSLT_VERSION})
endif((NOT HAVE_LIBXML2) OR (NOT HAVE_LIBXSLT))

if((NOT HAVE_LIBXML2) OR (NOT HAVE_LIBXSLT))
  redist_package(LIBXML2 ${BCL2FASTQ_LIBXML2_VERSION} 
                 "--prefix=${REINSTDIR};--without-modules;--without-http;--without-ftp;--without-python;--without-threads;--without-schematron;--without-debug;--without-iconv")
  find_library_redist(LIBXML2 ${REINSTDIR} libxml/xpath.h xml2)
  redist_package(LIBXSLT ${BCL2FASTQ_LIBXSLT_VERSION} "--prefix=${REINSTDIR};--with-libxml-prefix=${REINSTDIR};--without-plugins;--without-crypto")
  find_library_redist(LIBEXSLT ${REINSTDIR} libexslt/exslt.h exslt)
  find_library_redist(LIBXSLT ${REINSTDIR} libxslt/xsltconfig.h xslt)
endif((NOT HAVE_LIBXML2) OR (NOT HAVE_LIBXSLT))

include_directories(BEFORE SYSTEM ${LIBXML2_INCLUDE_DIR})
include_directories(BEFORE SYSTEM ${LIBXSLT_INCLUDE_DIR})
include_directories(BEFORE SYSTEM ${LIBEXSLT_INCLUDE_DIR})
set(BCL2FASTQ_DEP_LIB ${BCL2FASTQ_DEP_LIB} "${LIBEXSLT_LIBRARIES}" "${LIBXSLT_LIBRARIES}" "${LIBXML2_LIBRARIES}")

#set (CMAKE_CXX_FLAGS "$ENV{CXX_FLAGS} $ENV{CXXFLAGS} -fopenmp -msse2 -Werror -Wall -Wextra -Wunused -Wno-long-long -Wsign-compare -Wpointer-arith" CACHE STRING "g++ flags" FORCE)
set (CMAKE_CXX_FLAGS "$ENV{CXX_FLAGS} $ENV{CXXFLAGS} -std=c++11 -fopenmp -msse2 -Wall -Wextra -Wunused -Wno-long-long -Wsign-compare -Wpointer-arith" CACHE STRING "g++ flags" FORCE)
#set (CMAKE_CXX_FLAGS_DEBUG "-O0 -g -pg -fprofile-arcs -ftest-coverage -D_GLIBCXX_DEBUG" CACHE STRING "g++ flags" FORCE)
set (CMAKE_CXX_FLAGS_DEBUG "-O0 -std=c++11 -g -pg -fprofile-arcs -ftest-coverage" CACHE STRING "g++ flags" FORCE)
#set (CMAKE_CXX_FLAGS_DEBUG "-std=c++11 -O0 -g -pg -fprofile-arcs -ftest-coverage" CACHE STRING "g++ flags" FORCE)
set (CMAKE_CXX_FLAGS_RELEASE "-O3 -std=c++11 -DNDEBUG" CACHE STRING "g++ flags" FORCE)
#set (CMAKE_CXX_FLAGS_RELEASE "-std=c++11 -O3 -DNDEBUG" CACHE STRING "g++ flags" FORCE)
set (CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -std=c++11 -g" CACHE STRING "g++ flags" FORCE)
#set (CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g -pg" CACHE STRING "g++ flags" FORCE)
#set (CMAKE_CXX_FLAGS_RELWITHDEBINFO "-std=c++11 -O3 -g" CACHE STRING "g++ flags" FORCE)
set (CMAKE_CXX_FLAGS_MINSIZEREL "-Os -std=c++11 -DNDEBUG" CACHE STRING "g++ flags" FORCE)
#set (CMAKE_CXX_FLAGS_MINSIZEREL "-std=c++11 -Os -DNDEBUG" CACHE STRING "g++ flags" FORCE)

# Force static linking
set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "")

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE version)
    string(STRIP ${version} version)

    string(REGEX REPLACE "^([0-9])\\.[0-9]\\.[0-9]" "\\1" major_version ${version})
    string(REGEX REPLACE "^[0-9]\\.([0-9])\\.[0-9]" "\\1" minor_version ${version})
    string(REGEX REPLACE "^[0-9]\\.[0-9]\\.([0-9])" "\\1" patch_version ${version})
    if    (major_version LESS 4 OR (major_version EQUAL 4 AND (minor_version LESS 1 OR (minor_version EQUAL 1 AND patch_version LESS 2) ) ) )
        message (FATAL_ERROR "Unsupported GNU C++ compiler: g++ version ${version}: "
                             "only g++ versions >= 4.1.2 are supported")
    endif (major_version LESS 4 OR (major_version EQUAL 4 AND (minor_version LESS 1 OR (minor_version EQUAL 1 AND patch_version LESS 2) ) ) )

    set("${CMAKE_CXX_COMPILER_ID}${major_version}" true)
    set("${CMAKE_CXX_COMPILER_ID}${major_version}${minor_version}" true)
    set("${CMAKE_CXX_COMPILER_ID}${major_version}${minor_version}${patch_version}" true)
    message (STATUS "using compiler: gcc version ${version}")

endif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")

##
## Suppress spurious warnings in less recent compilers
##
if    (NOT GNU42)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-parameter ")
endif (NOT GNU42)

if    (GNU412 OR GNU42 OR GNU43)
    ## Before 4.1.2, pedantic breaks on boost lambda expressions
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic ")
endif (GNU412 OR GNU42 OR GNU43)

if (CMAKE_SYSTEM_PROCESSOR MATCHES "^i[67]86$")
    ##
    ## Use scalar floating point instructions from the SSE instruction set.
    ## Note: Pentium3 SSE supports only single precision arithmetics
    ##
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse -mfpmath=sse")
endif (CMAKE_SYSTEM_PROCESSOR MATCHES "^i[67]86$")
if (CMAKE_SYSTEM_PROCESSOR MATCHES "^i[345]86$")
    ##
    ## Prevent using 80bits registers (more consistent rounding)
    ##
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffloat-store")
endif (CMAKE_SYSTEM_PROCESSOR MATCHES "^i[345]86$")

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/lib/common/config.h.in ${BCL2FASTQ_CXX_CONFIG_H_DIR}/config.h)
