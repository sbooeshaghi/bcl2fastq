/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file FileReader.hh
 *
 * \brief Declaration of FileReader file.
 *
 * \author Aaron Day
 */

#ifndef BCL2FASTQ_DATA_FILEREADER_HPP
#define BCL2FASTQ_DATA_FILEREADER_HPP

#include <boost/function.hpp>

namespace bcl2fastq
{
namespace data
{

template<typename HEADER>
bool BinaryFileReader::readHeader(HEADER& header)
{
    std::streamsize headerLength = this->read(reinterpret_cast<char *>(&header), sizeof(header));

    int errnum = errno;
    if (headerLength != sizeof(header))
    {
        boost::format errFormat(boost::format(
            "Unable to read header of '%s' file '%s' : bytes_read=%d bytes_expected=%d") %
            this->getFileTypeStr() % this->getPath().string() % headerLength % sizeof(header));

        if (ignoreErrors_)
        {
            BCL2FASTQ_LOG(common::LogLevel::WARNING)
                << errFormat
                << ":" << std::strerror(errnum) << " (" << errnum << ")" << std::endl;

            return false;
        }
        else
        {
            BOOST_THROW_EXCEPTION(
                common::InputDataError(errnum,
                                       errFormat.str()));
        }
    }

    return true;
}

} // namespace data
} // namespace bcl2fastq

#endif // BCL2FASTQ_DATA_FILEREADER_HPP

