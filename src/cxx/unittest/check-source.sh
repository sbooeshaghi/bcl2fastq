#!/bin/bash
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
## file check-source.sh
##
## Basic sanity checks on the source code before running the cppuint tests
##
## author Come Raczy
##
################################################################################

target=$1
shift
good=yes
for file in $* ; do
    check=`grep -nH CPPUNIT_TEST_SUITE_NAMED_REGISTRATION $file | grep -v registryName`
    if [[ $check ]] ; then
        if [[ $good ]] ; then
            good=
            echo >&2
            echo use of unchecked registry names: >&2
        fi
        echo "    "$check >&2
        echo >&2
    fi
done
if [[ $good ]] ; then
    echo checked > $target
else 
    exit 1
fi


