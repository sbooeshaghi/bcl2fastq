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
## Configuration file for the cxx/lib subfolder.
##
## author Roman Petrovski
##
################################################################################


##
## List of BCL2FASTQ libraries
##
set (BCL2FASTQ_ALL_LIBRARIES
    common
    config
    io
    data
    layout
    stats
    conversion
)

##
## Build all the libraries for the project
## BCL2FASTQ_AVAILABLE_LIBRARIES is incrementally updated
##

set (BCL2FASTQ_AVAILABLE_LIBRARIES "")
set (BCL2FASTQ_ALL_LIBRARY_DIRS "")
foreach (BCL2FASTQ_LIB_DIR ${BCL2FASTQ_ALL_LIBRARIES})
    add_subdirectory(${BCL2FASTQ_LIB_DIR})
    set(BCL2FASTQ_AVAILABLE_LIBRARIES bcl2fastq_${BCL2FASTQ_LIB_DIR} ${BCL2FASTQ_AVAILABLE_LIBRARIES} )
endforeach (BCL2FASTQ_LIB_DIR)

set (BCL2FASTQ_AVAILABLE_LIBRARIES ${BCL2FASTQ_AVAILABLE_LIBRARIES} PARENT_SCOPE)


