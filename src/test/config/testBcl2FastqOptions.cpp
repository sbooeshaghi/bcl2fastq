// testBcl2FastqOptions

#include <vector>
#include <string>
#include <boost/format.hpp>
#include <iostream>
#include <ios>

#include "common/Debug.hh"
#include "common/SystemCompatibility.hh"
#include "common/Options.hh"
#include "config/Bcl2FastqOptions.hh"
#include "io/SyncFile.hh"
#include "common/Logger.hh"
#include "data/InteropFile.hh"
#include "gtest/gtest.h"

using namespace bcl2fastq;

namespace {

    typedef  bool (config::Bcl2FastqOptions::*BooleanFn0)() const;
    typedef  const boost::filesystem::path& (config::Bcl2FastqOptions::*PathFn0)() const;

    struct TestOptions {

        bool operator()(
            const std::vector<std::string> &arguments,
            config::Bcl2FastqOptions &options,
            BooleanFn0 memberFn) {

            const char **argv = new const char* [arguments.size()+2];
            int argc = 0;
            argv [argc++] = "bcl2fastq";
            for (std::string value : arguments) {
                argv [argc++] = value.c_str();
            }

            argv [argc++] = NULL;  // end of arguments sentinel is NULL

            options.parse(argc-1, argv);
            return (options.*memberFn)();
        }

        std::string operator()(
            const std::vector<std::string> &arguments,
            config::Bcl2FastqOptions &options,
            PathFn0 memberFn) {

            const char **argv = new const char* [arguments.size()+2];
            int argc = 0;
            argv [argc++] = "bcl2fastq";
            for (std::string value : arguments) {
                argv [argc++] = value.c_str();
            }

            argv [argc++] = NULL;  // end of arguments sentinel is NULL

            options.parse(argc-1, argv);
            return (options.*memberFn)().string();
        }

    };

} // end unnamed namespace


//=========================================
//
// TESTING A
//
//=========================================

//TEST(Bcl2fastqOptions, IgnoreMissing) {
//
//    TestOptions tester;
//    config::Bcl2FastqOptions options;
//
//    common::LogLevel::value_type oldLevel = common::BCL2FASTQ_MIN_LOG_LEVEL;
//    common::BCL2FASTQ_MIN_LOG_LEVEL = common::LogLevel::NONE;
//
//    EXPECT_FALSE(tester({                             }, options, &config::Bcl2FastqOptions::getIgnoreMissingDataFiles));
//    EXPECT_TRUE (tester({"--ignore-missing-data-files"}, options, &config::Bcl2FastqOptions::getIgnoreMissingDataFiles));
//
//    // Deprecated options
//    EXPECT_FALSE(tester({                             }, options, &config::Bcl2FastqOptions::getIgnoreMissingDataFiles));
//    EXPECT_TRUE (tester({"--ignore-missing-bcls"      }, options, &config::Bcl2FastqOptions::getIgnoreMissingDataFiles));
//    EXPECT_TRUE (tester({"--ignore-missing-filter"    }, options, &config::Bcl2FastqOptions::getIgnoreMissingDataFiles));
//    EXPECT_TRUE (tester({"--ignore-missing-positions" }, options, &config::Bcl2FastqOptions::getIgnoreMissingDataFiles));
//    EXPECT_TRUE (tester({"--ignore-missing-controls"  }, options, &config::Bcl2FastqOptions::getIgnoreMissingDataFiles));
//    EXPECT_FALSE(tester({                             }, options, &config::Bcl2FastqOptions::getIgnoreMissingDataFiles));
//
//    common::BCL2FASTQ_MIN_LOG_LEVEL = oldLevel;
//}

TEST(Bcl2fastqOptions, SettingRunfolder) {

    TestOptions tester;
    config::Bcl2FastqOptions options;

    common::LogLevel::value_type oldLevel = common::BCL2FASTQ_MIN_LOG_LEVEL;
    common::BCL2FASTQ_MIN_LOG_LEVEL = common::LogLevel::NONE;

//    EXPECT_EQ ("\"$CWD\"",
//        tester({ },
//        options, &config::Bcl2FastqOptions::getRunfolderDir));
    EXPECT_EQ ("/tmp/runfolder",
        tester({"--runfolder-dir", "/tmp/runfolder" },
        options, &config::Bcl2FastqOptions::getRunfolderDir));

    common::BCL2FASTQ_MIN_LOG_LEVEL = oldLevel;
}

TEST(Bcl2fastqOptions, SettingInputDir) {

    TestOptions tester;
    config::Bcl2FastqOptions options;

    common::LogLevel::value_type oldLevel = common::BCL2FASTQ_MIN_LOG_LEVEL;
    common::BCL2FASTQ_MIN_LOG_LEVEL = common::LogLevel::NONE;

    EXPECT_EQ ("/tmp/runfolder/Data/Intensities/BaseCalls/",
        tester({"--runfolder-dir", "/tmp/runfolder"},
        options, &config::Bcl2FastqOptions::getInputDir));
    EXPECT_EQ ("/tmp/runfolder/BaseCalls",
        tester({"--input-dir", "/tmp/runfolder/BaseCalls" },
        options, &config::Bcl2FastqOptions::getInputDir));

    common::BCL2FASTQ_MIN_LOG_LEVEL = oldLevel;
}

TEST(Bcl2fastqOptions, SettingIntensitiesDir) {

    TestOptions tester;
    config::Bcl2FastqOptions options;

    common::LogLevel::value_type oldLevel = common::BCL2FASTQ_MIN_LOG_LEVEL;
    common::BCL2FASTQ_MIN_LOG_LEVEL = common::LogLevel::NONE;

    EXPECT_EQ ("/tmp/runfolder/Data/Intensities/",
        tester({"--runfolder-dir", "/tmp/runfolder"},
        options, &config::Bcl2FastqOptions::getIntensitiesDir));
    EXPECT_EQ ("/tmp/runfolder/Data/Intensities",
        tester({"--intensities-dir", "/tmp/runfolder/Data/Intensities" },
        options, &config::Bcl2FastqOptions::getIntensitiesDir));

    common::BCL2FASTQ_MIN_LOG_LEVEL = oldLevel;
}

TEST(Bcl2fastqOptions, SettingOutputDir) {

    TestOptions tester;
    config::Bcl2FastqOptions options;

    common::LogLevel::value_type oldLevel = common::BCL2FASTQ_MIN_LOG_LEVEL;
    common::BCL2FASTQ_MIN_LOG_LEVEL = common::LogLevel::NONE;

    EXPECT_EQ ("/tmp/runfolder/Data/Intensities/BaseCalls/",
        tester({"--runfolder-dir", "/tmp/runfolder"},
        options, &config::Bcl2FastqOptions::getOutputDir));
    EXPECT_EQ ("/tmp/output",
        tester({"--output-dir", "/tmp/output" },
        options, &config::Bcl2FastqOptions::getOutputDir));

    common::BCL2FASTQ_MIN_LOG_LEVEL = oldLevel;
}

TEST(Bcl2fastqOptions, SettingStatsDir) {

    TestOptions tester;
    config::Bcl2FastqOptions options;

    common::LogLevel::value_type oldLevel = common::BCL2FASTQ_MIN_LOG_LEVEL;
    common::BCL2FASTQ_MIN_LOG_LEVEL = common::LogLevel::NONE;

//    EXPECT_EQ ("$CWD/Stats",
//        tester({ },
//        options, &config::Bcl2FastqOptions::getStatsDir));
    EXPECT_EQ ("/tmp/output/Stats/",
        tester({"--output-dir", "/tmp/output" },
        options, &config::Bcl2FastqOptions::getStatsDir));
    EXPECT_EQ ("/tmp/output/stats2",
        tester({"--stats-dir", "/tmp/output/stats2" },
        options, &config::Bcl2FastqOptions::getStatsDir));

    common::BCL2FASTQ_MIN_LOG_LEVEL = oldLevel;
}

TEST(Bcl2fastqOptions, SettingReportsDir) {

    TestOptions tester;
    config::Bcl2FastqOptions options;

    common::LogLevel::value_type oldLevel = common::BCL2FASTQ_MIN_LOG_LEVEL;
    common::BCL2FASTQ_MIN_LOG_LEVEL = common::LogLevel::NONE;

    EXPECT_EQ ("/tmp/output/Reports/",
        tester({"--output-dir", "/tmp/output" },
        options, &config::Bcl2FastqOptions::getReportsDir));
    EXPECT_EQ ("/tmp/output/report2",
        tester({"--reports-dir", "/tmp/output/report2" },
        options, &config::Bcl2FastqOptions::getReportsDir));

    common::BCL2FASTQ_MIN_LOG_LEVEL = oldLevel;
}

TEST(Bcl2fastqOptions, SettingInterOpDir) {

    TestOptions tester;
    config::Bcl2FastqOptions options;

    common::LogLevel::value_type oldLevel = common::BCL2FASTQ_MIN_LOG_LEVEL;
    common::BCL2FASTQ_MIN_LOG_LEVEL = common::LogLevel::NONE;

    EXPECT_EQ ("/tmp/runfolder/InterOp/",
        tester({"--runfolder-dir", "/tmp/runfolder"},
        options, &config::Bcl2FastqOptions::getInteropDir));
    EXPECT_EQ ("/tmp/interop",
        tester({"--interop-dir", "/tmp/interop" },
        options, &config::Bcl2FastqOptions::getInteropDir));

    common::BCL2FASTQ_MIN_LOG_LEVEL = oldLevel;
}

TEST(Bcl2fastqOptions, SampleSheetPath) {

    TestOptions tester;
    config::Bcl2FastqOptions options;

    common::LogLevel::value_type oldLevel = common::BCL2FASTQ_MIN_LOG_LEVEL;
    common::BCL2FASTQ_MIN_LOG_LEVEL = common::LogLevel::NONE;

    EXPECT_EQ ("/tmp/MySampleSheet.csv",
        tester({"--sample-sheet", "/tmp/MySampleSheet.csv" },
        options, &config::Bcl2FastqOptions::getSampleSheetPath));

    common::BCL2FASTQ_MIN_LOG_LEVEL = oldLevel;
}

