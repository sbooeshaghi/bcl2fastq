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
## Configuration file for the cxx/bin subdirectory.
##
## author Come Raczy
##
################################################################################


include(${BCL2FASTQ_CXX_EXECUTABLE_CMAKE})

file (GLOB BCL2FASTQ_PROGRAM_SOURCE_LIST [a-zA-Z0-9]*.cpp)

##
## Generic rule for all the other programs
##
foreach(BCL2FASTQ_PROGRAM_SOURCE ${BCL2FASTQ_PROGRAM_SOURCE_LIST})
    get_filename_component(BCL2FASTQ_PROGRAM ${BCL2FASTQ_PROGRAM_SOURCE} NAME_WE)
    add_executable        (${BCL2FASTQ_PROGRAM} ${BCL2FASTQ_PROGRAM_SOURCE})
    target_link_libraries (${BCL2FASTQ_PROGRAM} ${BCL2FASTQ_AVAILABLE_LIBRARIES}
                           ${Boost_LIBRARIES}
                           ${BCL2FASTQ_DEP_LIB}
                           ${BCL2FASTQ_ADDITIONAL_LIB} )
    install(TARGETS ${BCL2FASTQ_PROGRAM} RUNTIME DESTINATION ${BCL2FASTQ_BINDIR})
endforeach(BCL2FASTQ_PROGRAM_SOURCE)


