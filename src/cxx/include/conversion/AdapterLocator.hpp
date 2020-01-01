/**
* BCL to FASTQ file converter
* Copyright (c) 2007-2017 Illumina, Inc.
*
* This software is covered by the accompanying EULA
* and certain third party copyright/licenses, and any user of this
* source file is bound by the terms therein.
*
* \file AdapterLocator.hpp
*
* \brief Implementation of adapter locator.
*
* \author Aaron Day
*/

#ifndef BCL2FASTQ_CONVERSION_ADAPTER_LOCATOR_HPP
#define BCL2FASTQ_CONVERSION_ADAPTER_LOCATOR_HPP

#include "conversion/AdapterLocator.hh"

#include "conversion/BclBaseConversion.hh"
#include "common/Debug.hh"

namespace bcl2fastq
{
namespace conversion
{

template<common::NumBasesPerByte numBasesPerByte>
AdapterLocator<numBasesPerByte>::AdapterLocator(const std::string& adapter,
                                                float              adapterStringency)
: adapter_(adapter)
, adapterStringency_(adapterStringency)
{
}

template<common::NumBasesPerByte numBasesPerByte>
AdapterLocatorSlidingWindow<numBasesPerByte>::AdapterLocatorSlidingWindow(const std::string& adapter,
                                                                          float              adapterStringency)
: AdapterLocator<numBasesPerByte>(adapter, adapterStringency)
{
}

template<common::NumBasesPerByte numBasesPerByte>
FastqConstIterator<numBasesPerByte> AdapterLocatorSlidingWindow<numBasesPerByte>::identifyAdapter(const FastqConstIterator<numBasesPerByte>& basesBegin,
                                                                                                  const FastqConstIterator<numBasesPerByte>& basesEnd) const
{
    const std::string::size_type adapterLength = this->adapter_.size();

    for (FastqConstIterator<numBasesPerByte> basesIter = basesBegin; basesIter != basesEnd; ++basesIter)
    {
        std::size_t mismatchCount = 0;
        std::size_t matchCount = 0;
        const size_t checkCycles = std::min(static_cast<typename FastqConstIterator<numBasesPerByte>::difference_type>(adapterLength), basesEnd - basesIter);
        for (std::size_t adapterOffset = 0; adapterOffset < checkCycles; ++adapterOffset)
        {
            const char base = convertBcl2FastqBase(basesIter[adapterOffset]);
            if (base == 'N')
            {
                // Ns are neither matches nor mismatches
                continue;
            }
            if (base == this->adapter_.at(adapterOffset))
            {
                ++matchCount;
            }
            else
            {
                ++mismatchCount;
            }

            if ((mismatchCount > 1) && (mismatchCount > matchCount)) // Shortcut
            {
                // If at any point there are 2 or more mismatches and the number of mismatches is
                // greater than the number of matches, stop scanning and start a new scan at
                // the next initial position in the read sequence
                break;
            }
        }

        float sequenceIdentity = matchCount / (float)(matchCount + mismatchCount);
        if (matchCount > 0 && sequenceIdentity >= this->adapterStringency_)
        {
            return basesIter;
        }
    }

    return basesEnd;
}

template<common::NumBasesPerByte numBasesPerByte>
AdapterLocatorWithIndels<numBasesPerByte>::AdapterLocatorWithIndels(const std::string& adapter,
                                                                    float              adapterStringency,
                                                                    uint32_t           minimumTrimmedReadLength)
: AdapterLocator<numBasesPerByte>(adapter, adapterStringency)
, minimumTrimmedReadLength_(minimumTrimmedReadLength)
, minLengthMatch_(0)
, mismatchesToAccountFor_(0)
, maxAllowedMismatches_(adapter.size())
, adapterCallsSubsets_(adapter.size())
, peqs_(adapter.size())
, pv0s_(adapter.size())
{
    BCL2FASTQ_ASSERT_MSG(adapter.size() <= 64,
                         "Adapter size must be 64 or less. Adapter: " << adapter);

    std::vector<unsigned char> adapterCalls;
    binarizeAdapterCalls(adapter,
                         adapterCalls);

    size_t adapterLength = adapter.size();

    minLengthMatch_ = std::min(adapterLength, static_cast<size_t>(std::ceil(adapterStringency * 10)));
    mismatchesToAccountFor_ = std::ceil(adapterLength * (1.0 - adapterStringency));

    // calculate all values between 1 and adapterLength (need shorter than
    // minLengthMatch if mismatches found that reduce the search range to less then minLengthMatch)
    for (size_t i = 0; i < adapterLength; ++i)
    {
        //maximum number of allowed mismatches for a full length match
        maxAllowedMismatches_[i] = std::floor((i+1) * (1.0 - this->adapterStringency_));  //rounding down

        adapterCallsSubsets_[i] = std::vector<unsigned char>(adapterCalls.begin(), adapterCalls.begin()+(i+1));

        //calculate Peqs and Pv0s
        peqs_[i] = getPeq(adapterCallsSubsets_[i]);

        pv0s_[i] = (2L << i) - 1;
    }
}

template<common::NumBasesPerByte numBasesPerByte>
FastqConstIterator<numBasesPerByte> AdapterLocatorWithIndels<numBasesPerByte>::identifyAdapter(const FastqConstIterator<numBasesPerByte>& basesBegin,
                                                                                               const FastqConstIterator<numBasesPerByte>& basesEnd) const
{
    if (basesBegin == basesEnd)
    {
        return basesEnd;
    }

    uint32_t minScorePartialMatch = this->adapter_.size();

    FastqConstIterator<numBasesPerByte> pos(findFullMatch(basesBegin,
                                            basesEnd,
                                            minScorePartialMatch));

    if (pos == basesEnd)
    {
        // try to find partial match at end since we did not find full match
        return findPartialMatch(basesBegin,
                                basesEnd,
                                minScorePartialMatch);
    }

    return pos;
}

template<common::NumBasesPerByte numBasesPerByte>
FastqConstIterator<numBasesPerByte> AdapterLocatorWithIndels<numBasesPerByte>::findPartialMatch(const FastqConstIterator<numBasesPerByte>& basesBegin,
                                                                                                const FastqConstIterator<numBasesPerByte>& basesEnd,
                                                                                                uint32_t minScorePartialMatch) const
{
    FastqConstIterator<numBasesPerByte> trimPosition = basesEnd;

    size_t adapterLength = this->adapter_.size();
    uint64_t minScore = adapterLength;

    // number of bases to consider for partial match
    uint32_t tmpLength = (minScorePartialMatch > mismatchesToAccountFor_ ?
        adapterLength - minScorePartialMatch + mismatchesToAccountFor_ - 1 :
        adapterLength - 1);

    if (tmpLength > std::distance(basesBegin, basesEnd))
    {
        return basesEnd;
    }

    while (trimPosition == basesEnd)
    {
        // only exact matches possible
        if (minLengthMatch_ >= tmpLength)
        {
            return findLongestExactSuffixMatch(basesBegin,
                                               basesEnd,
                                               adapterCallsSubsets_[tmpLength - 1],
                                               tmpLength);
        }

        // mismatches still allowed
        MyersValues myersValues(tmpLength, pv0s_[tmpLength - 1]);
        uint32_t maxAllowedMismatch = maxAllowedMismatches_[tmpLength - 1];

        BOOST_REVERSE_FOREACH(char base, std::make_pair(basesBegin, basesEnd))
        {
            myersValues.updateValues(base, peqs_[tmpLength - 1]);

            if (myersValues.score_ <= minScore)
            {
                minScore = myersValues.score_;

               // partial adapter match
               if (maxAllowedMismatch >= myersValues.score_)
               {
                   trimPosition = basesEnd - myersValues.count_;
               }
            }

            if (myersValues.count_ > myersValues.compareLength_)
            {
                break; // end foreach
            }
        }

        if (trimPosition != basesEnd || minScore == tmpLength) // match found or all bases are mismatching
        {
            break; //end while loop
        }

        tmpLength -= (minScore > maxAllowedMismatch ?
                      (minScore - maxAllowedMismatch) :
                      1);
    }

    return trimPosition;
}

template<common::NumBasesPerByte numBasesPerByte>
FastqConstIterator<numBasesPerByte> AdapterLocatorWithIndels<numBasesPerByte>::findFullMatch(const FastqConstIterator<numBasesPerByte>& basesBegin,
                                                                                             const FastqConstIterator<numBasesPerByte>& basesEnd,
                                                                                             uint32_t&                                  minScorePartialMatch) const
{
    size_t adapterLength = this->adapter_.size();
    MyersValues myersValues(adapterLength, pv0s_[adapterLength - 1]);

    uint32_t minScore = adapterLength;
    uint32_t cycleCount = std::distance(basesBegin, basesEnd);
    int trimPosition = -1;
    uint32_t lengthWithMismatch = cycleCount + maxAllowedMismatches_[adapterLength - 1];

    //search for full matching adapters
    for (FastqConstIterator<numBasesPerByte> pos = basesEnd - 1; pos >= basesBegin; --pos)
    {
        myersValues.updateValues(*pos, peqs_[adapterLength - 1]);

        //full adapter match
        if (myersValues.score_ <= maxAllowedMismatches_[adapterLength - 1] &&
            (myersValues.score_ <= minScore || myersValues.count_ > (lengthWithMismatch - trimPosition)))
        {
            trimPosition = cycleCount - myersValues.count_;

            // NOTE: There used to be code here to break out of this loop 
            // early, assuming that it could stop if there was a good enough
            // match that will trim the read shorter than the required minimum
            // read length. This was flawed, because indels and N's mean there
            // could be more matches that would mask more.
        }

        //match with less mismatches than previous best match
        if (myersValues.score_ <= minScore)
        {
            minScore = myersValues.score_;
            //maximum number of matches within adapter length from end (used to find partial match)
            if (myersValues.count_ < adapterLength)
            {
                minScorePartialMatch = minScore;
            }

        }
    }

    return (trimPosition == -1 ? basesEnd : basesBegin + trimPosition);
}

template<common::NumBasesPerByte numBasesPerByte>
AdapterLocatorWithIndels<numBasesPerByte>::MyersValues::MyersValues(uint32_t length,
                                                                    uint64_t pv)
: mv_(0)
, count_(0)
, score_(length)
, compareLength_(length-1)
, pv_(pv)
{
}

template<common::NumBasesPerByte numBasesPerByte>
void AdapterLocatorWithIndels<numBasesPerByte>::MyersValues::updateValues(unsigned char b, const std::vector<uint64_t>& peq)
{
    // all zero for N and otherwise take the two right-most bits for basecalls (others are for quality calls)
    uint64_t eq_ = peq[b];
    uint64_t d0 =  (((eq_ & pv_) + pv_) ^ pv_) | eq_ | mv_;

    uint64_t ph = mv_ | ~(d0 | pv_);
    uint64_t mh = pv_ & d0;

    score_ += (ph >> compareLength_ & 1) - (mh >> compareLength_ & 1);

    ph <<= 1;
    pv_ = (mh << 1) | ~(d0 | ph);
    mv_ = ph & d0;

    //zero-based position from 3'-end or step count
    ++count_;
}

template<common::NumBasesPerByte numBasesPerByte>
std::vector<uint64_t> AdapterLocatorWithIndels<numBasesPerByte>::getPeq(const std::vector<unsigned char>& adapterCalls) const
{
    std::vector<uint64_t> peq(256, 0);
    for (int peqBase = 1; peqBase <= 255; ++peqBase)
    {
        uint64_t& value = peq[peqBase];
        BOOST_FOREACH(unsigned char base, std::make_pair(adapterCalls.begin(), adapterCalls.end()))
        {
            value <<= 1;

            if ((peqBase & 3) == base)
            {
                value |= 1L;
            }
        }
    }

    return peq;
}

template<common::NumBasesPerByte numBasesPerByte>
FastqConstIterator<numBasesPerByte> AdapterLocatorWithIndels<numBasesPerByte>::findLongestExactSuffixMatch(const FastqConstIterator<numBasesPerByte>& basesBegin,
                                                                                                           const FastqConstIterator<numBasesPerByte>& basesEnd,
                                                                                                           const std::vector<unsigned char>& adapterSubstring,
                                                                                                           int tagLength) const
{
    int cycleCount = std::distance(basesBegin, basesEnd);

    for (int i = tagLength - 1; i >= 0; --i)
    {
        int match = 0;
        int tmp = cycleCount - 1 - i;
        for (int j = 0; j <= i; ++j)
        {
            if (( *(basesBegin + tmp + j) & 3) == adapterSubstring[j])
            {
                ++match;
            }
            else
            {
                break;
            }
        }

        if (match > i)
        {
            return basesEnd - match;
        }
    }

    return basesEnd;
}

template<common::NumBasesPerByte numBasesPerByte>
void AdapterLocatorWithIndels<numBasesPerByte>::binarizeAdapterCalls(const std::string& adapter,
                                                                     std::vector<unsigned char>& adapterCalls) const
{
    adapterCalls.resize(adapter.size());
    for (size_t i = 0; i < adapter.size(); ++i)
    {
        switch (adapter[i])
        {
            case 'A':
            case 'a':
                adapterCalls[i] = 0;
                break;
            case 'C':
            case 'c':
                adapterCalls[i] = 1;
                break;
            case 'G':
            case 'g':
                adapterCalls[i] = 2;
                break;
            case 'T':
            case 't':
                adapterCalls[i] = 3;
                break;
            default:
                break;
        }
    }
}

} // namepsace conversion
} // namespace bcl2fastq

#endif // BCL2FASTQ_CONVERSION_ADAPTER_LOCATOR_HPP

