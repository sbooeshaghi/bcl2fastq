/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file Utility.cpp
 *
 * \brief Implementation of various I/O related utilities.
 *
 * \author Roman Petorvski
 * \author Marek Balint
 */


#include <vector>

#include <boost/assign/list_of.hpp>

#include "io/Utility.hh"


namespace bcl2fastq {


namespace io {


namespace detail {


static const std::vector<const char *> iosBaseToStdioOpenModesTranslationTable = boost::assign::list_of<const char*>
                // inspired by /usr/include/c++/4.3.2/fstream
                // | ios_base Flag combination      |
                // |binary  in  out  trunc  app     |
    (0      )   // |                                |
    ("a"    )   // |                         +      |
    (0      )   // |                   +            |
    (0      )   // |                   +     +      |
    ("w"    )   // |             +                  |
    ("a"    )   // |             +           +      |
    ("w"    )   // |             +     +            |
    (0      )   // |             +     +     +      |
    ("r"    )   // |         +                      |
    ("a+"   )   // |         +               +      |
    (0      )   // |         +         +            |
    (0      )   // |         +         +     +      |
    ("r+"   )   // |         +   +                  |
    ("a+"   )   // |         +   +           +      |
    ("w+"   )   // |         +   +     +            |
    (0      )   // |         +   +     +     +      |
    (0      )   // |   +                            |
    ("ab"   )   // |   +                     +      |
    (0      )   // |   +               +            |
    (0      )   // |   +               +     +      |
    ("wb"   )   // |   +         +                  |
    ("ab"   )   // |   +         +           +      |
    ("wb"   )   // |   +         +     +            |
    (0      )   // |   +         +     +     +      |
    ("rb"   )   // |   +     +                      |
    ("a+b"  )   // |   +     +               +      |
    (0      )   // |   +     +         +            |
    (0      )   // |   +     +         +     +      |
    ("r+b"  )   // |   +     +   +                  |
    ("a+b"  )   // |   +     +   +           +      |
    ("w+b"  )   // |   +     +   +     +            |
    (0      )   // |   +     +   +     +     +      |
;


} // namespace detail


const char * iosFlagsToStdioMode(std::ios_base::openmode mode)
{
    const unsigned openModeIndex =
        !!(mode & std::ios_base::binary) << 4 |
        !!(mode & std::ios_base::in    ) << 3 |
        !!(mode & std::ios_base::out   ) << 2 |
        !!(mode & std::ios_base::trunc ) << 1 |
        !!(mode & std::ios_base::app   ) << 0
    ;

    return detail::iosBaseToStdioOpenModesTranslationTable.at(openModeIndex);
}


} // namespace io


} // namespace bcl2fastq


