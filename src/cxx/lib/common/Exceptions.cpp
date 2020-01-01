/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file Exceptions.cpp
 *
 * \brief Implementation of the common exception mechanism.
 *
 * \author Come Raczy
 * \author Marek Balint
 */


#include <cstring>

#include <boost/exception/diagnostic_information.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/format.hpp>

#include "common/Exceptions.hh"


namespace bcl2fastq {


namespace common {


Exception::Exception(int errorNumber, const std::string &message)
: boost::exception()
, errorNumber_(errorNumber)
, message_(message)
{
}

Exception::Exception(const Exception &that)
: boost::exception(that)
, errorNumber_(that.errorNumber_)
, message_(that.message_)
{
}

Exception::~Exception() throw()
{
}

// Restrict assignment operator.
#if 0
Exception & Exception::operator=(const Exception &)
{
}
#endif

int Exception::getErrorNumber() const
{
    return errorNumber_;
}

const std::string &Exception::getMessage() const
{
    return message_;
}

std::string Exception::getContext() const
{
    const std::string now = boost::posix_time::to_simple_string(
        boost::posix_time::second_clock::local_time()
    );

    return (boost::format("%s: %s (%d): %s")
        % now
        % std::string(std::strerror(errorNumber_))
        % errorNumber_
        % boost::diagnostic_information(*this)
    ).str();
}


RuntimeError::RuntimeError(int errorNumber, const std::string &message)
: Exception(errorNumber, message)
, std::runtime_error(message + ": unknown error " + boost::lexical_cast<std::string>(errorNumber))
{
}

RuntimeError::RuntimeError(const std::string &message)
: Exception(0, message)
, std::runtime_error(message)
{
}

RuntimeError::RuntimeError(const RuntimeError &that)
: Exception(that)
, std::runtime_error(that)
{
}

const char * RuntimeError::what() const throw()
{
    return this->getMessage().c_str();
}


LogicError::LogicError(int errorNumber, const std::string &message)
: Exception(errorNumber, message)
, std::logic_error(message)
{
}

LogicError::LogicError(const std::string &message)
: Exception(0, message)
, std::logic_error(message)
{
}

LogicError::LogicError(const LogicError &that)
: Exception(that)
, std::logic_error(that)
{
}

const char * LogicError::what() const throw()
{
    return this->getMessage().c_str();
}


IoError::IoError(int errorNumber, const std::string &message)
: Exception(errorNumber, message)
, std::ios_base::failure(message)
{
}

IoError::IoError(const std::string &message)
: Exception(0, message)
, std::ios_base::failure(message)
{
}

IoError::IoError(const IoError &that)
: Exception(that)
, std::ios_base::failure(that)
{
}

const char * IoError::what() const throw()
{
    return this->getMessage().c_str();
}


InputDataError::InputDataError(int errorNumber, const std::string &message)
: RuntimeError(errorNumber, message)
{
}

InputDataError::InputDataError(const std::string &message)
: RuntimeError(0, message)
{
}

InputDataError::InputDataError(const InputDataError &that)
: RuntimeError(that)
{
}


OutOfRange::OutOfRange(int errorNumber, const std::string &message)
: LogicError(errorNumber, message)
{
}

OutOfRange::OutOfRange(const std::string &message)
: LogicError(message)
{
}

OutOfRange::OutOfRange(const OutOfRange &that)
: LogicError(that)
{
}


} // namespace comomn


} // namespace bcl2fastq


