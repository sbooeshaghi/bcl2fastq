/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file Options.cpp
 *
 * \brief Implementation of the program options processing.
 *
 * \author Come Raczy
 * \author Mauricio Varea
 * \author Marek Balint
 */


#include "common/Logger.hh"
#include "common/Options.hh"

#ifndef WIN32
#include <sys/ioctl.h>
#endif

namespace bcl2fastq {


namespace common {

static const unsigned short min_description_length = 50;

static unsigned short getStderrLineLength()
{
#ifndef WIN32
    winsize ret = {0,0,0,0};
    if (-1 == ioctl(STDERR_FILENO, TIOCGWINSZ, &ret))
    {
        return boost::program_options::options_description::m_default_line_length;
    }
    return ret.ws_col;
#else
    return boost::program_options::options_description::m_default_line_length;
#endif
}


Options::Options()
: optionsInitialized_(false)
, namedOptions_("Command-line options",getStderrLineLength(),getStderrLineLength()-min_description_length)
, unnamedOptions_()
, positionalOptions_()
{
}

Options::~Options()
{
}

void Options::initOptions()
{
    if (optionsInitialized_)
    {
        return;
    }

    namedOptions_.add_options()(
        "help,h",
        "produce help message and exit"
    );
    namedOptions_.add_options()(
        "version,v",
        "print program version information"
    );
    namedOptions_.add_options()(
        "min-log-level,l",
        boost::program_options::value<LogLevel::value_type>(&BCL2FASTQ_MIN_LOG_LEVEL)
            ->default_value(BCL2FASTQ_MIN_LOG_LEVEL),
        "minimum log level\n"
        "recognized values: NONE, FATAL, ERROR, WARNING, INFO, DEBUG, TRACE"
    );

    this->initNamedOptions(namedOptions_);
    this->initUnnamedOptions(unnamedOptions_);
    this->initPositionalOptions(positionalOptions_);

    optionsInitialized_ = true;
}

std::string Options::usage() const
{
    std::ostringstream os;

    os << this->usagePrefix() << std::endl;
    os << namedOptions_ << std::endl;
    os << this->usageSuffix() << std::endl;

    return os.str();
}

std::string Options::usageSuffix() const
{
    return "";
}

Options::Action::value_type Options::parse(int argc, const char *argv[])
{
    boost::program_options::variables_map vm;

    try
    {
        this->initOptions();
        boost::program_options::options_description allOptions("Allowed options");
        allOptions.add(namedOptions_);
        allOptions.add(unnamedOptions_);

        boost::program_options::store(
            boost::program_options::command_line_parser(argc, argv)
                .options(allOptions)
                .positional(positionalOptions_)
                .run(),
            vm
        );

        boost::program_options::notify(vm);

        if (vm.count("help"))
        {
            return Action::HELP;
        }
        else if (vm.count("version"))
        {
            return Action::VERSION;
        }

        std::stringstream log;
        log << "Command-line invocation: ";
        std::copy(&argv[0], &argv[argc], std::ostream_iterator<std::string>(log, " "));
        BCL2FASTQ_LOG(LogLevel::NONE) << log.str() << std::endl;
        BCL2FASTQ_LOG(LogLevel::INFO) << "Minimum log level: " << BCL2FASTQ_MIN_LOG_LEVEL << std::endl;
    }
    catch(const boost::program_options::multiple_values &e)
    {
        std::clog << this->usage();
        std::clog
            << "Failed to parse the options: " << e.what() << ": "
            << e.get_option_name() << std::endl
        ;
        return Action::ABORT;
    }
    catch (const boost::program_options::multiple_occurrences &e)
    {
        std::clog << this->usage() << std::endl;
        std::clog
            << "Failed to parse the options: " << e.what() << ": "
            << e.get_option_name() << std::endl
        ;
        return Action::ABORT;
    }
    catch (const boost::program_options::required_option &e)
    {
        std::clog << this->usage() << std::endl;
        std::clog
            << "Failed to parse the options: " << e.what() << ": "
            << e.get_option_name() << std::endl
        ;
        return Action::ABORT;
    }
    catch (const std::exception &e)
    {
        std::clog << this->usage() << std::endl;
        std::clog
            << "Failed to parse the options: "
            << e.what() << std::endl
        ;
        return Action::ABORT;
    }

    // If we were using a standard ini file, we could call this:
    // store(parse_config_file("example.cfg", desc), vm);
    this->storeSampleSheetSettings(vm);

    this->postProcess(vm);

    return Action::RUN;
}


std::ostream& operator<<(std::ostream& os, Options::Action::value_type a)
{
    switch(a) {
        case Options::Action::RUN:
        {
            os << "RUN";
            break;
        }
        case Options::Action::HELP:
        {
            os << "HELP";
            break;
        }
        case Options::Action::VERSION:
        {
            os << "VERSION";
            break;
        }
        case Options::Action::ABORT:
        {
            os << "ABORT";
            break;
        }
    }

    return os;
}


} // namespace comomn


} // namespace bcl2fastq


