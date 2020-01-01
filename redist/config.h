/**
 ** BCL to FASTQ file converter
 ** Copyright (c) 2007-2017 Illumina, Inc.
 **
 ** This software is covered by the accompanying EULA
 ** and certain third party copyright/licenses, and any user of this
 ** source file is bound by the terms therein.
 **
 ** \file config.h
 **
 ** \brief Various system-specific definitions (configured by cmake)
 **
 **/

/* cpp/include/config.h.cmake. Manually edited */
/* cpp/include/config.h.in.  Generated from configure.ac by autoheader.  */

#ifndef HG_BCL2FASTQ_COMMON_CONFIG_H
#define HG_BCL2FASTQ_COMMON_CONFIG_H


/* Endianness of the architecture */
/* #undef BCL2FASTQ_IS_BIG_ENDIAN */

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 if you have the <malloc.h> header file. */
#define HAVE_MALLOC_H 1

/* Define to 1 if you have the <mcheck.h> header file. */
/* #undef HAVE_MCHECK_H */

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define to 1 if you have the <signal.h> header file. */
#define HAVE_SIGNAL_H 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <time.h> header file. */
#define HAVE_TIME_H 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the `floorf' function. */
#define HAVE_FLOORF 1

/* Define to 1 if you have the `round' function. */
#define HAVE_ROUND 1

/* Define to 1 if you have the `roundf' function. */
#define HAVE_ROUNDF 1

/* Define to 1 if you have the `powf' function. */
#define HAVE_POWF 1

/* Define to 1 if you have the `erf' function. */
#define HAVE_ERF 1

/* Define to 1 if you have the `erff' function. */
#define HAVE_ERFF 1

/* Define to 1 if you have the `erfc' function. */
#define HAVE_ERFC 1

/* Define to 1 if you have the `erfcf' function. */
#define HAVE_ERFCF 1

/* Define to 1 if you have the numa library. */
/* #undef HAVE_NUMA */

/* Define to 1 if you have the `zlib' library */
#define HAVE_ZLIB 1

/* Define to 1 if you have the `sysconf' library */
#define HAVE_SYSCONF 1

/* Define to 1 if you have the `clock' library */
#define HAVE_CLOCK 1

/* Define to 1 if you have the `bzip2' library */
/* #undef HAVE_BZIP2 */
/* #undef HAVE_BZLIB */

/* Define to 1 if you have the `fftw3f' library */
/* #undef HAVE_FFTW3F */

/* Define to 1 if you have the `cpgplot' library */
/* #undef HAVE_CPGPLOT */

/* Define to 1 if you have the `pgplot' library */
/* #undef HAVE_PGPLOT */

/* Define to 1 if you have the `X11' library */
#define HAVE_X11 1

/* Define to 1 if you have the `g2c' library */
/* #undef HAVE_G2C */

/* Define to 1 if you have the `boost_xxx_yyy' library
   (-lboost_xxx_yyy). */
#define HAVE_LIBBOOST_DATE_TIME 1
#define HAVE_LIBBOOST_FILESYSTEM 1
#define HAVE_LIBBOOST_IOSTREAMS 1
#define HAVE_LIBBOOST_PROGRAM_OPTIONS 1
/* #undef HAVE_LIBBOOST_PYTHON */
#define HAVE_LIBBOOST_REGEX 1
#define HAVE_LIBBOOST_SERIALIZATION 1
#define HAVE_LIBBOOST_SYSTEM 1

/* Define to 1 if you have the `cppunit' library (-lcppunit). */
#define HAVE_CPPUNIT 1

/* Name of package */
/* #undef PACKAGE */

/* Top level namespace */
/* #undef NAMESPACE */

/* Define to the address where bug reports for this package should be sent. */
/* #undef PACKAGE_BUGREPORT bcl2fastq_bug@illumina.com */

/* Short name */
#define BCL2FASTQ_NAME_SHORT "bcl2fastq"

/* Long name */
#define BCL2FASTQ_NAME_LONG "BCL to FASTQ file converter"

/* Version number */
#define BCL2FASTQ_VERSION "2.18.1.13W"

/* Version number: marketing */
#define BCL2FASTQ_VERSION_MAJOR "2"

/* Version number: year */
#define BCL2FASTQ_VERSION_MINOR "18"

/* Version number: month */
#define BCL2FASTQ_VERSION_PATCH "1"

/* Version number: day */
#define BCL2FASTQ_VERSION_BUILD "13"

/* Copyright one-liner */
#define BCL2FASTQ_COPYRIGHT "Copyright (c) 2007-2017 Illumina, Inc."

/* Location of installed share files */
#define BCL2FASTQ_DATADIR "share"

/* Location of installed bin files */
#define BCL2FASTQ_VERSION_FULL ""

/* Define to empty if `const' does not conform to ANSI C. */
#undef const


#endif /* HG_BCL2FASTQ_COMMON_CONFIG_H */


