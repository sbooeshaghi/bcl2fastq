/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file Exceptions.hh
 *
 * \brief Declaration of the common exception mechanism.
 *
 * \author Come Raczy
 * \author Marek Balint
 */


#ifndef BCL2FASTQ_COMMON_EXCEPTIONS_HH
#define BCL2FASTQ_COMMON_EXCEPTIONS_HH


#include <string>
#include <stdexcept>
#include <ios>

#include <boost/exception/exception.hpp>
// definition of BOOST_EXCEPTION_THROW macro
// (not directly used by this file, included for user convenience)
#include <boost/throw_exception.hpp>


namespace bcl2fastq {


namespace common {


/// \brief Virtual base class to all the exception classes in Bcl2FastQ. 
/// \note Use BOOST_THROW_EXCEPTION to get the contect info
/// (file, function, line) at the throw site. 
class Exception : public boost::exception
{
public:

    /// \brief Constructor.
    /// \param errorNumber Error number.
    /// \param message Error mesasge.
    Exception(int errorNumber, const std::string &message);

    /// \brief Copy constructor.
    /// \param that Other exception object to copy from.
    Exception(const Exception &that);

    /// \brief Virtual destructor.
    /// \par Exception guarantee:
    /// no-throw
    virtual ~Exception() throw();

private:

    /// \brief Restricted assignment operator.
    Exception & operator=(const Exception &);

public:

    /// \brief Get error number.
    /// \return Error number.
    int getErrorNumber() const;

    /// \brief Gete error message.
    /// \return Error message.
    const std::string &getMessage() const;

    /// \brief Get error context.
    /// \return Error context.
    std::string getContext() const;

private:

    /// \brief Error number.
    const int errorNumber_;

    /// \brief Error message.
    const std::string message_;
};


/// \brief Exception for errors presumably detectable only when the program executes.
class RuntimeError : public Exception, public std::runtime_error
{
public:

    /// \brief Constructor.
    /// \param errorNumber Error number.
    /// \param message Error message.
    RuntimeError(int errorNumber, const std::string &message);

    /// \brief Constructor.
    /// \param message Error message.
    explicit RuntimeError(const std::string &message);

    /// \brief Copy constructor.
    /// \param that Other exception object to copy from.
    RuntimeError(const RuntimeError &that);

public:

    /// \brief Get error message.
    /// \par Exception guarantee:
    /// no-throw
    virtual const char * what() const throw();
};


/// \brief Exception for errors presumably detectable before the program executes, such as violations of logical preconditions or class invariants.
class LogicError : public Exception, public std::logic_error
{
public:

    /// \brief Constructor.
    /// \param errorNumber Error number.
    /// \param message Error message.
    LogicError(int errorNumber, const std::string &message);

    /// \brief Constructor.
    /// \param message Error message.
    explicit LogicError(const std::string &message);

    /// \brief Copy constructor.
    /// \param that Other exception object to copy from.
    LogicError(const LogicError &that);

public:

    /// \brief Get error message.
    /// \par Exception guarantee:
    /// no-throw
    virtual const char * what() const throw();
};


/// \brief Exception thrown when there are problems with the IO operations.
/// \todo Refactoring: In C++11, @c ios_base::failure is no longer derived
/// directly from @c std::exception, but rather from std::system_error,
/// which is it turn derived from std::runtime_error.
class IoError: public Exception, public std::ios_base::failure
{
public:

    /// \brief Constructor.
    /// \param errorNumber Error number.
    /// \param message Error message.
    IoError(int errorNumber, const std::string &message);

    /// \brief Constructor.
    /// \param message Error message.
    explicit IoError(const std::string &message);

    /// \brief Copy constructor.
    /// \param that Other exception object to copy from.
    IoError(const IoError &that);

public:

    /// \brief Get error message.
    /// \par Exception guarantee:
    /// no-throw
    virtual const char * what() const throw();
};


/// \brief Exception for input data errors.
class InputDataError : public RuntimeError
{
public:

    /// \brief Constructor.
    /// \param errorNumber Error number.
    /// \param message Error message.
    InputDataError(int errorNumber, const std::string &message);

    /// \brief Constructor.
    /// \param message Error message.
    explicit InputDataError(const std::string &message);

    /// \brief Copy constructor.
    /// \param that Other exception object to copy from.
    InputDataError(const InputDataError &that);
};


/// \brief Exception for out-of-range errors.
class OutOfRange : public LogicError
{
public:

    /// \brief Constructor.
    /// \param errorNumber Error number.
    /// \param message Error message.
    OutOfRange(int errorNumber, const std::string &message);

    /// \brief Constructor.
    /// \param message Error message.
    explicit OutOfRange(const std::string &message);

    /// \brief Copy constructor.
    /// \param that Other exception object to copy from.
    OutOfRange(const OutOfRange &that);
};




} // namespace common


} // namespace bcl2fastq


#endif // BCL2FASTQ_COMMON_EXCEPTIONS_HH


