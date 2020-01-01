#ifndef HELPER_BCL2FASTQTESTHELPER_HH

#include <vector>
#include <string>

namespace bcl2fastq {
namespace helper {

    struct Bcl2FastqOptionsOrator;
    typedef std::shared_ptr<Bcl2FastqOptionsOrator> Bcl2FastqOptionsOratorPtr;

    struct Bcl2FastqOptionsOrator 
    {
        config::Bcl2FastqOptions options_;
        common::LogLevel::value_type oldLevel_ = common::BCL2FASTQ_MIN_LOG_LEVEL;

        virtual ~Bcl2FastqOptionsOrator() {
            common::BCL2FASTQ_MIN_LOG_LEVEL = oldLevel_;
        }

        Bcl2FastqOptionsOrator(
            const std::vector<std::string> &arguments,
            common::LogLevel::value_type logLevel = common::LogLevel::WARNING)
        : options_()
        {
            common::BCL2FASTQ_MIN_LOG_LEVEL = logLevel;

            const char **argv = new const char* [arguments.size()+2];
            int argc = 0;
            argv [argc++] = "bcl2fastq";
            for (std::string value : arguments) {
                argv [argc++] = value.c_str();
            }

            argv [argc++] = NULL;  // end of arguments sentinel is NULL

            options_.parse(argc-1, argv);
        }

        static Bcl2FastqOptionsOratorPtr createOptionsWithRunFolder(const std::string &runFolder) {
             std::vector<std::string> arguments( {"--runfolder-dir", runFolder} );
             return std::make_shared<helper::Bcl2FastqOptionsOrator>( arguments );
        }

    }; // end struct Bcl2FastqOptionsOrator

} // end helper namespace
} // end bcl2fastq namespace

#endif // HELPER_BCL2FASTQTESTHELPER_HH
