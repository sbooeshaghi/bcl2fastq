#!/bin/bash
version=`grep "\<BCL2FASTQ_VERSION\>" build/cxx/common/config.h | cut -d " " -f 3 | tr -d "\042"`
