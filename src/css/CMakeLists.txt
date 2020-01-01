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
## Configuration file for the css subfolder
##
## author Mauricio Varea
##
################################################################################
include ("${BCL2FASTQ_GLOBALS_CMAKE}")

install(CODE "
    include (\"${BCL2FASTQ_MACROS_CMAKE}\")
    configure_files_recursively (\"${CMAKE_CURRENT_SOURCE_DIR}\" \"${CMAKE_CURRENT_BINARY_DIR}\" \"*.css\")
    install_files_recursively (\"${CMAKE_CURRENT_BINARY_DIR}\" \"${BCL2FASTQ_ORIG_DATADIR}/css\" \"*.css\" \"\${BCL2FASTQ_LIBRARY_PERMISSIONS}\")
    configure_files_recursively (\"${CMAKE_CURRENT_SOURCE_DIR}\" \"${CMAKE_CURRENT_BINARY_DIR}\" \"*.css.xml\")
    install_files_recursively (\"${CMAKE_CURRENT_BINARY_DIR}\" \"${BCL2FASTQ_ORIG_DATADIR}/css\" \"*.css.xml\" \"\${BCL2FASTQ_LIBRARY_PERMISSIONS}\")
    ")

