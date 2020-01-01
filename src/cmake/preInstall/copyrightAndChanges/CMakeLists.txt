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
## Configuration file for the COPYRIGHT and Changes file installation.
##
## author Roman Petrovski
##
################################################################################


# Installing top level components
install(FILES "${BCL2FASTQ_SOURCE_DIR}/Changes" "${CMAKE_SOURCE_DIR}/../COPYRIGHT" DESTINATION "${BCL2FASTQ_DATADIR}")


