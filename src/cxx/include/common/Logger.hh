/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file Logger.hh
 *
 * \brief Preprocessor-based logger mechanism.
 *
 * \author Roman Petrovski
 * \author Marek Balint
 */


#include <ios>
#include <string>

#include <boost/atomic.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/thread.hpp>


#ifndef BCL2FASTQ_COMMON_LOGGER_HH
#define BCL2FASTQ_COMMON_LOGGER_HH


/// \brief Logging stream.
/// \param level Log level to be used.
/// \todo Refactoring: Define separate macro for each logging level (e.g. @c BCL2FASTQ_LOG_ERROR, @c BCL2FASTQ_LOG_INFO, ...).
#define BCL2FASTQ_LOG(level)                                                                       \
    if (::bcl2fastq::common::detail::LogStream s = ::bcl2fastq::common::detail::LogStream(level))  \
    {                                                                                              \
    }                                                                                              \
    else                                                                                           \
        std::clog                                                                                  \
            << ::bcl2fastq::common::detail::ThreadTimestamp()                                      \
            /* << std::setfill(' ') << std::setw(7) */ << level << (level ? ": " : "")             \
            << (                                                                                   \
                   /*(level == ::bcl2fastq::common::LogLevel::TRACE)                                 \
                   ?                                                                               \
                       std::string(__FILE__) + ":"                                                 \
                       + boost::lexical_cast<std::string>(__LINE__) + ": "                         \
                   :*/                                                                               \
                       ""                                                                          \
               )                                                                                   \
        /* no semicolon (';') here */


namespace bcl2fastq {


namespace common {


/// \brief Logging level.
struct LogLevel
{
    enum value_type
    {
        NONE = 0,    ///< Suppress all logging.
        FATAL = 1,   ///< Fatal errors only.
        ERROR_TYPE = 2,   ///< All errors.
        WARNING = 3, ///< Warnings and all errors.
        INFO = 4,    ///< Informative messages and warnings/errors.
        DEBUG = 5,   ///< Debug messages, informative messages and warnings/errors.
        TRACE = 6    ///< Log everything.
    };
};

/// \brief Input operator for logging level.
/// \param is Input stream.
/// \param l Loging level to be input.
/// \return Input stream.
std::istream& operator>>(std::istream& is, LogLevel::value_type &l);

/// \brief Output operator for logging level.
/// \param os Output stream.
/// \param l Loging level to be output.
/// \return Output stream.
std::ostream& operator<<(std::ostream& os, LogLevel::value_type l);

/// \brief Minimum log level of the program.
extern LogLevel::value_type BCL2FASTQ_MIN_LOG_LEVEL;


namespace detail {

/// \brief Stream for logging.
class LogStream
{
public:

    /// \brief Default constructor.
    LogStream();

    /// \brief Constructor.
    /// \param logLevel Log level to be used.
    explicit LogStream(LogLevel::value_type logLevel);

    /// \brief Copy constructor.
    /// \param that The other instance.
    LogStream(const LogStream &that);

    /// \brief Conversion operator.
    operator bool () const;

    static unsigned int getNumWarnings() { return numWarnings_; }
    static unsigned int getNumErrors() { return numErrors_; }

private:

    /// \brief Mutex.
    static boost::recursive_mutex clogMutex_;
    static unsigned int numWarnings_;
    static unsigned int numErrors_;

    /// \brief Lock guard.
    boost::lock_guard<boost::recursive_mutex> lock_;

    /// \brief Log stream saved state.
    boost::io::ios_base_all_saver ias_;

    /// \brief Log level to be used.
    LogLevel::value_type logLevel_;
};


/// \brief Thread timestamp tag for output operator overload.
struct ThreadTimestamp
{
};

/// \brief Formats time stamp and thread id to simplify threaded logging.
/// \param os Output stream.
/// \return Output stream.
std::ostream & operator << (std::ostream &os, const ThreadTimestamp &);


} // namespace detail


} // namespace common


} // namespace bcl2fastq


#endif // BCL2FASTQ_COMMON_LOGGER_HH


