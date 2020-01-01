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
## file checkTargetPathWritable.cmake
##
## Script checking writable state of installation path.
##
## author Roman Petrovski
##
################################################################################


foreach (BCL2FASTQ_DEST_DIR ${BCL2FASTQ_DEST_DIRS})
    message (STATUS "Testing access to ${BCL2FASTQ_DEST_DIR}...")
    execute_process(COMMAND "${CMAKE_COMMAND}" "-E" "make_directory" "${BCL2FASTQ_DEST_DIR}/test" RESULT_VARIABLE TMP_RESULT )
    if (TMP_RESULT)
        message (STATUS "ERROR: Directory is not writeable: ${BCL2FASTQ_DEST_DIR}")
        message (STATUS "If you don't have administrator access to the "
                         "target installation location, please use --prefix "
                         "command-line option when configuring BCL2FASTQ. "
                         "Please use configure --help for all installer "
                         "command-line options details.")
        message (FATAL_ERROR "ERROR: BCL2FASTQ installation cannot continue")
    else (TMP_RESULT)
        message (STATUS "Directory is writeable: ${BCL2FASTQ_DEST_DIR}")
    endif (TMP_RESULT)
endforeach (BCL2FASTQ_DEST_DIR ${BCL2FASTQ_DEST_DIRS})


