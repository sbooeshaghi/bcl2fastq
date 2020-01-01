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
## file CMakeLists.txt
##
## Configuration file for the checkTargetPathsWritable subdirectory.
##
## author Roman Petrovski
##
################################################################################


message (STATUS "Verifying target directories access")

include ("${BCL2FASTQ_GLOBALS_CMAKE}")
install(
    CODE "set(BCL2FASTQ_DEST_DIRS \"\${BCL2FASTQ_DEST_ETCDIR}\" \"\${BCL2FASTQ_DEST_DATADIR}\" 
                                  \"\${BCL2FASTQ_DEST_BINDIR}\" \"\${BCL2FASTQ_DEST_LIBDIR}\" 
                                  \"\${BCL2FASTQ_DEST_LIBEXECDIR}\")"
    SCRIPT "${CMAKE_CURRENT_SOURCE_DIR}/checkTargetPathWritable.cmake"
)


