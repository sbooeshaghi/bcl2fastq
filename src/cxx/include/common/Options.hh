/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file Options.hh
 *
 * \brief Declaration of the program options processing.
 *
 * \author Come Raczy
 * \author Mauricio Varea
 * \author Marek Balint
 */


#ifndef BCL2FASTQ_COMMON_OPTIONS_HH
#define BCL2FASTQ_COMMON_OPTIONS_HH


#include <boost/utility.hpp>
#include <boost/program_options.hpp>


namespace bcl2fastq {


namespace common {


/// \brief Encapsulation of the processing of the command line options.
/// \todo Features: Add config file and environment options.
class Options : private boost::noncopyable
{

public:

    /// \brief Program action.
    struct Action
    {
            enum value_type
            {
                RUN,     ///< Run.
                HELP,    ///< Display help and exit.
                VERSION, ///< Display version and exit.
                ABORT    ///< Abort.
            };
    };

public:

    /// \brief Default constructor.
    Options();

    /// \brief Virtual destructor.
    virtual ~Options() = 0;

private:

    /// \brief Initialize program options.
    void initOptions();

protected:

    /// \brief Template method: initialize named options.
    virtual void initNamedOptions(boost::program_options::options_description &options) = 0;

    /// \brief Template method: initialize unnamed options.
    virtual void initUnnamedOptions(boost::program_options::options_description &options) = 0;

    /// \brief Tmeplate method: initialize positional options.
    virtual void initPositionalOptions(boost::program_options::positional_options_description &options) = 0;

public:

    /// \brief Get program usage help.
    /// \return String containing program usage help.
    std::string usage() const;

protected:

    /// \brief Get program usage prefix.
    /// \return String containing program usage prefix.
    virtual std::string usagePrefix() const = 0;

    /// \brief Get program usage suffix.
    /// \return String containing program usage suffix.
    virtual std::string usageSuffix() const;

public:

    /// \brief Parse command line arguments.
    /// \param argc Number of command line arguments.
    /// \param argv Command line arguments.
    /// \return Action to be taked based on arguments.
    Options::Action::value_type parse(int argc, const char *argv[]);

protected:

    /// \brief Store the options from the sample sheet.
    /// \param vm Options map where the settings are stored.
    virtual void storeSampleSheetSettings(boost::program_options::variables_map &vm) = 0;

    /// \brief Post-process parsed command line arguments.
    /// \param vm Command line arguments to be post-process.
    virtual void postProcess(boost::program_options::variables_map &vm) = 0;

private:

    /// \brief Flag determining whether option description has been initialized.
    bool optionsInitialized_;

    /// \brief Named options.
    boost::program_options::options_description namedOptions_;

    /// \brief Unnamed options.
    boost::program_options::options_description unnamedOptions_;

    /// \brief Positional options.
    boost::program_options::positional_options_description positionalOptions_;
};


/// \brief Output operator for program Action.
/// \param os Output stream.
/// \param a Program action to be output.
/// \return Output stream.
std::ostream& operator<<(std::ostream& os, Options::Action::value_type a);


} // namespace common


} // namespace bcl2fastq


#endif // BCL2FASTQ_COMMON_OPTIONS_HH


