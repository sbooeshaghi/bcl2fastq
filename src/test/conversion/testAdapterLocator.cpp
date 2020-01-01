#include "data/BclBuffer.hh"
#include "conversion/AdapterLocator.hh"

#include "gtest/gtest.h"

using namespace bcl2fastq;

namespace {

void transposeToBclContainer(const std::string& read, data::BclBuffer &bclBuffer) {

    bclBuffer.bcls_.clear();
    bclBuffer.bcls_.resize(read.size());

    for (size_t i = 0; i < read.size(); ++i) {
        switch (read[i]) {

            case 'A':
                bclBuffer.bcls_[i].bcls_.push_back('@');
                break;
            case 'C':
                bclBuffer.bcls_[i].bcls_.push_back('A');
                break;
            case 'G':
                bclBuffer.bcls_[i].bcls_.push_back('B');
                break;
            case 'T':
                bclBuffer.bcls_[i].bcls_.push_back('C');
                break;
            case 'N':
                bclBuffer.bcls_[i].bcls_.push_back('\0');
                break;
            default:
                BOOST_ASSERT(false);
                break;
        }
    }
}

class AdapterFinder {

public:

    bool indels_;

    AdapterFinder(bool indel = false) 
    : indels_(indel)
    {}

    size_t operator()(
        const std::string &read,
        const std::string &trimAdapter,
        float              trimAdapterStringency,
        uint32_t           minimumTrimmedReadLength = -1) {

        data::BclBuffer bclBuffer;
        transposeToBclContainer(read, bclBuffer);

        unsigned int i = 0;
        conversion::FastqConstIterator<common::NumBasesPerByte::ONE> begin(bclBuffer, i);
        conversion::FastqConstIterator<common::NumBasesPerByte::ONE> end(bclBuffer, i);
        end += bclBuffer.bcls_.size();

        size_t distance = 0;
        if (indels_) {
   
            conversion::AdapterLocatorWithIndels<common::NumBasesPerByte::ONE> trimAdapterLocator(
                trimAdapter, trimAdapterStringency, minimumTrimmedReadLength);
            conversion::FastqConstIterator<common::NumBasesPerByte::ONE> pos = trimAdapterLocator.identifyAdapter(begin, end);
            distance = std::distance(begin, pos);
        }
        else {

            conversion::AdapterLocatorSlidingWindow<common::NumBasesPerByte::ONE> trimAdapterLocator(
                trimAdapter, trimAdapterStringency);
            conversion::FastqConstIterator<common::NumBasesPerByte::ONE> pos = trimAdapterLocator.identifyAdapter(begin, end);
            distance = std::distance(begin, pos);
        }

        return distance;
    }

}; // end class AdapterFinder

} // end unnamed namespace

TEST(SlidingWindowAdapterLocator, FindAdapter) { 

    // NOTE: The N used in these tests could be any sequence but it's easier to see the test like this. 
    AdapterFinder finder(false);

    const std::string a1 = "CTGTCTCTTATACACATCT";

    // Find trim adapter //0....+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8....+....9....+....0
    EXPECT_EQ( 77, finder("NNNNNNNNNNNNNNNNNNNNNNNNNGGCCGGGAGTTGGGCGAGTACGGGCTGCAGGCATACACTGAAGTGAAAACTGCTGTCTCTTATACACATCTTCCG", a1, 0.8));
                         //                                                                             -------------------
    EXPECT_EQ(  7, finder("NNNNNNNNNNNNNNNNNNNNNNNNNTATCCTCTTTCTTCTCCTTTCCCTTTTCTTCCCCCTTNNNNNNNNCNNNNNNNNNNNNNNCTGTCTCTTATACAC", a1, 0.8));
                         //       -------------------                                                           
    EXPECT_EQ( 64, finder("NNNNNNNNNNNNNNNNNNNNNNNNNGATCCTCTTTCTTCTCCTTTCCCTTTTCTTCCCCCTTNNNNNNNNCNNNNNNNNNNNNNNCTGTCTCTTATACAC", a1, 0.8));
                         //                                                                -------------------
    EXPECT_EQ( 85, finder("NNNNNNNNNNNNNNNNNNNNNNNNNGATCCTCTTTCTTCTCCTTTCCCTTTTCTTCCCCCTTNNNNNNNNCAAANNNNNNNNNNNCTGTCTCTTATACAC", a1, 0.8));
                         //                                                                                     ---------------...
}

TEST(IndelAdapterLocator, FindAdapter) {

    // NOTE: The N used in these tests could be any sequence but it's easier to see the test like this. 
    AdapterFinder finder(true);

    // Find trim adapter //0....+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8....+....9....+....0
    const std::string a1 = "CTGTCTCTTATACACATCT";
    EXPECT_EQ( 77, finder("NNNNNNNNNNNNNNNNNNNNNNNNNGGCCGGGAGTTGGGCGAGTACGGGCTGCAGGCATACACTGAAGTGAAAACTGCTGTCTCTTATTACACATCTCCG", a1, 0.8, 30));
                         //                                                                             -------------------
    EXPECT_EQ( 85, finder("NNNNNNNNNNNNNNNNNNNNNNNNNTATCCTCTTTCTTCTCCTTTCCCTTTTCTTCCCCCTTNNNNNNNNCNNNNNNNNNNNNNNCTGTCTCTTATACAC", a1, 0.8, 30));
                         //                                                                                     ---------------...
    EXPECT_EQ(100, finder("NNNNNNNNNNNNNNNNNNNNNNNNNNACGTAGGTATTCTTTCTTTCTTCTTTTTCCTTTTCTTCCTTNCTCCTNNNNNNNNNNNNNNNNNNNNNNNNNNN", a1, 0.8, 30));
                         //                                                                                                    ...
    EXPECT_EQ( 54, finder("NNNNNNNNNNNNNNNNNNNNNNNNNCGTCCACACCCTCCCCGGCCAGCCATGCACCGTTCTCTTATACACATCTCCGAGCCCACGAGACAGGGAGAAATT", a1, 0.8, 30));
                         //                                                      -------------------
    EXPECT_EQ( 67, finder("NNNNNNNNNNNNNNNNNNNNNNNNNAAGCACAGGGCAATCATATTGTGCTTTTCACTGTCTCCAGGACATGTCTCTTATACACATCTCCGAGCCCACGAG", a1, 0.8, 30));
                         //                                                                   -------------------
    EXPECT_EQ( 38, finder("NNNNNNNNNNNNNNNNNNNNNNNGCCTCATTGCCACAGCTGTTCTCTTATACACATCTCCGAGCCCACGAGACTAAGGCGAATCTCGTATGCCGTCTTCT", a1, 0.8, 30));
                         //                                      -------------------
    EXPECT_EQ( 91, finder("NNNNNNNNNNNNNNNNNNNNNNNAAAACCACTTTGTCCTCATCTCATTCCCCAAAACCGGCAGGCCCTCTCAGTTCTTTTTCCTGTCTTTACTTCTCTTA", a1, 0.8, 30));
                         //                                                                                           ---------...
    EXPECT_EQ( 99, finder("NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNC", a1, 0.8, 30));
                         //                                                                                                   -...

    const std::string a2 = "CTGTCTCTTATACACATCTACTGCTAGTCGCGTTACGGTT";
    // Long trim adapter //0....+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8....+....9....+....0
    EXPECT_EQ( 99, finder("NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNC", a2, 0.8, 30));
    EXPECT_EQ( 45, finder("NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCTGTCTCTTATACACATCTACTGCTAGTCGCGTTACGGTTNNNNNNNNNNNNNNN", a2, 0.8, 30));

    // Trim adapter at beginning
    const std::string a3 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
    EXPECT_EQ(  0, finder("ATCGGAAGACACACGTCTGAACTCCAGTCACATTCCTCGATCTCGTATGCCGTTTTCTGCTTGAAAAAACCAACAACACCTGGGAACAAGCGAATCACCN", a3, 0.9, 35));
                      //...---------d--------------------
}

TEST(IndelAdapterLocator, InvalidAdapter) {

    AdapterFinder finder;
    // Interesting...? SampleSheet loader would fail with this trim adaptor
    EXPECT_EQ(100, finder("ATCGGAAGACACACGTCTGAACTCCAGTCACATTCCTCGATCTCGTATGCCGTTTTCTGCTTGAAAAAACCAACAACACCTGGGAACAAGCGAATCACCT", "SDFGHJK", 0.9, 30));
}

