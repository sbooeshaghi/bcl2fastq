/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file Logger.cpp
 *
 * \brief Preprocessor-based logger mechanism.
 *
 * \author Roman Petrovski
 * \author Marek Balint
 */


#include "common/Logger.hh"


namespace bcl2fastq {


namespace common {

std::istream& operator>>(std::istream& is, LogLevel::value_type &l)
{
    std::string token;
    is >> token;

    if ((token == "none") || (token == "NONE") || (token == "0"))
    {
        l = LogLevel::NONE;
    }
    else if ((token == "fatal") || (token == "FATAL") || (token == "1"))
    {
        l = LogLevel::FATAL;
    }
    else if ((token == "error") || (token == "ERROR") || (token == "2"))
    {
        l = LogLevel::ERROR_TYPE;
    }
    else if ((token == "warning") || (token == "WARNING") || (token == "3"))
    {
        l = LogLevel::WARNING;
    }
    else if ((token == "info") || (token == "INFO") || (token == "4"))
    {
        l = LogLevel::INFO;
    }
    else if ((token == "debug") || (token == "DEBUG") || (token == "5"))
    {
        l = LogLevel::DEBUG;
    }
    else if ((token == "trace") || (token == "TRACE") || (token == "6"))
    {
        l = LogLevel::TRACE;
    }
    else
    {
        is.setstate(std::ios_base::failbit);
    }

    return is;
}

std::ostream& operator<<(std::ostream& os, LogLevel::value_type l)
{
    switch(l) {
        case LogLevel::NONE:
        {
            //os << "NONE";
            break;
        }
        case LogLevel::FATAL:
        {
            os << "FATAL";
            break;
        }
        case LogLevel::ERROR_TYPE:
        {
            os << "ERROR";
            break;
        }
        case LogLevel::WARNING:
        {
            os << "WARNING";
            break;
        }
        case LogLevel::INFO:
        {
            os << "INFO";
            break;
        }
        case LogLevel::DEBUG:
        {
            os << "DEBUG";
            break;
        }
        case LogLevel::TRACE:
        {
            os << "TRACE";
            break;
        }
    }

    return os;
}

LogLevel::value_type BCL2FASTQ_MIN_LOG_LEVEL = LogLevel::INFO;


namespace detail {


boost::recursive_mutex LogStream::clogMutex_;
unsigned int LogStream::numErrors_ = 0;
unsigned int LogStream::numWarnings_ = 0;

LogStream::LogStream()
: lock_(clogMutex_)
, ias_(std::clog)
, logLevel_(LogLevel::WARNING)
{
    ++numWarnings_;
}

LogStream::LogStream(LogLevel::value_type logLevel)
: lock_(clogMutex_)
, ias_(std::clog)
, logLevel_(logLevel)
{
    if (logLevel_ == ::bcl2fastq::common::LogLevel::WARNING)
    {
        ++numWarnings_;
    }
    else if (logLevel_ == ::bcl2fastq::common::LogLevel::ERROR_TYPE ||
             logLevel_ == ::bcl2fastq::common::LogLevel::FATAL)
    {
        ++numErrors_;
    }
}

LogStream::LogStream(const LogStream &that)
: lock_(clogMutex_)
, ias_(std::clog)
, logLevel_(that.logLevel_)
{
}

LogStream::operator bool () const
{
    return logLevel_ > BCL2FASTQ_MIN_LOG_LEVEL;
}

std::ostream & operator <<(std::ostream &os, const ThreadTimestamp &) {

    const boost::posix_time::ptime currentTime = boost::posix_time::second_clock::local_time();

    // IMPORTANT: this is the way to serialize date without causing any dynamic memory operations to occur
    const std::tm t = boost::posix_time::to_tm(currentTime.time_of_day());
    const boost::posix_time::ptime::date_type d = currentTime.date();
    return os
        << d.year() << '-'
        << std::setfill('0') << std::setw(2) << d.month().as_number() << '-'
        << std::setfill('0') << std::setw(2) << d.day() << ' '
        << std::setfill('0') << std::setw(2) << t.tm_hour << ':'
        << std::setfill('0') << std::setw(2) << t.tm_min << ':'
        << std::setfill('0') << std::setw(2) << t.tm_sec << ' '
        << '[' << boost::this_thread::get_id() << ']' << ' '
    ;
}


} // namespace detail


} // namespace comomn


} // namespace bcl2fastq


