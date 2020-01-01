/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file Program.hpp
 *
 * \brief Implementation of the skeleton of all c++ programs.
 *
 * \author Come Raczy
 * \author Mauricio Varea
 * \author Marek Balint
 */


#ifndef BCL2FASTQ_COMMON_PROGRAM_HPP
#define BCL2FASTQ_COMMON_PROGRAM_HPP


#include <iostream>
#include <sstream>
#include <algorithm>
#include <string>
#include <cstdlib>

#include <boost/cstdlib.hpp>
#include <boost/program_options.hpp>

#include "common/ProgramInfo.hh"
#include "common/Exceptions.hh"
#include "common/Logger.hh"


namespace bcl2fastq {


namespace common {


template<typename O>
void run(void(*callback)(const ProgramInfo &, O &), int argc, const char *argv[])
{
#ifndef WIN32
    setenv("LC_ALL", "C", 1);
#endif

    try {
        ProgramInfo info(argv[0], boost::program_options::options_description::m_default_line_length);
        std::clog << info << std::endl;

        O options;
        const typename O::Action::value_type action = options.parse(argc, argv);

        switch(action) {
            case O::Action::RUN:
            {
                callback(info, options);

                BCL2FASTQ_LOG(LogLevel::NONE) << "Processing completed with " << ::bcl2fastq::common::detail::LogStream::getNumErrors()
                          << " errors and " << ::bcl2fastq::common::detail::LogStream::getNumWarnings()
                          << " warnings." << std::endl;
                break;
            }
            case O::Action::HELP:
            {
                std::clog << options.usage();
                std::exit(boost::exit_success);
            }
            case O::Action::VERSION:
            {
                std::exit(boost::exit_success);
            }
            case O::Action::ABORT:
            {
                std::exit(boost::exit_failure);
            }
        }
    }
    catch (const Exception &e)
    {
        BCL2FASTQ_LOG(LogLevel::ERROR_TYPE)
            << "bcl2fastq::common::Exception: "
            << e.getContext()
            << std::endl
        ;
        std::exit(boost::exit_failure);
    }
    catch (const boost::exception &e)
    {
        BCL2FASTQ_LOG(LogLevel::ERROR_TYPE)
            << "boost::exception: "
            << boost::diagnostic_information(e)
            << std::endl
        ;
        std::exit(boost::exit_failure);
    }
    catch (const std::runtime_error &e)
    {
        BCL2FASTQ_LOG(LogLevel::ERROR_TYPE)
            << "std::runtime_error: "
            << e.what()
            << std::endl
        ;
        std::exit(boost::exit_failure);
    }
    catch (const std::logic_error &e)
    {
        BCL2FASTQ_LOG(LogLevel::ERROR_TYPE)
            << "std::logic_error: "
            << e.what()
            << std::endl
        ;
        std::exit(boost::exit_failure);
    }
    catch(...) {
        BCL2FASTQ_LOG(LogLevel::FATAL)
            << "Unknown exception"
            << std::endl
        ;
        std::exit(boost::exit_failure);
    }
}


} // namespace common


} // namespace bcl2fastq


#endif // BCL2FASTQ_COMMON_PROGRAM_HH


