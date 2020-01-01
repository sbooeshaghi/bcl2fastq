/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file AdapterLocator.hh
 *
 * \brief Declaration of Adapter Locator.
 *
 * \author Aaron Day
 */


#ifndef BCL2FASTQ_CONVERSION_ADAPTER_LOCATOR_HH
#define BCL2FASTQ_CONVERSION_ADAPTER_LOCATOR_HH

#include "conversion/FastqIterator.hh"
#include "common/Types.hh"

#include <boost/noncopyable.hpp>

namespace bcl2fastq
{
namespace conversion
{

template<common::NumBasesPerByte numBasesPerByte>
class AdapterLocator : private boost::noncopyable
{
public:
    AdapterLocator(const std::string& adapter,
                   float              adapterStringency);

    virtual ~AdapterLocator() { }

    virtual FastqConstIterator<numBasesPerByte> identifyAdapter(const FastqConstIterator<numBasesPerByte>& basesBegin,
                                                                const FastqConstIterator<numBasesPerByte>& basesEnd) const = 0;

protected:
    const std::string adapter_;
    float             adapterStringency_;
};

template<common::NumBasesPerByte numBasesPerByte>
class AdapterLocatorSlidingWindow : public AdapterLocator<numBasesPerByte>
{
public:
    AdapterLocatorSlidingWindow(const std::string& adapter,
                                float              adapterStringency);

    virtual ~AdapterLocatorSlidingWindow() { }

    virtual FastqConstIterator<numBasesPerByte> identifyAdapter(const FastqConstIterator<numBasesPerByte>& basesBegin,
                                                                const FastqConstIterator<numBasesPerByte>& basesEnd) const;
};

template<common::NumBasesPerByte numBasesPerByte>
class AdapterLocatorWithIndels : public AdapterLocator<numBasesPerByte>
{
public:
    AdapterLocatorWithIndels(const std::string& adapter,
                             float              adapterStringency,
                             uint32_t           minimumTrimmedReadLength);

    virtual ~AdapterLocatorWithIndels() { }

    virtual FastqConstIterator<numBasesPerByte> identifyAdapter(const FastqConstIterator<numBasesPerByte>& basesBegin,
                                                                const FastqConstIterator<numBasesPerByte>& basesEnd) const;

private:
    struct MyersValues : private boost::noncopyable
    {
        MyersValues(uint32_t length,
                    uint64_t pv);

        void updateValues(unsigned char b, const std::vector<uint64_t>& peq);

        uint64_t mv_;
        uint64_t count_;
        uint64_t score_;
        uint64_t compareLength_;
        uint64_t pv_;
    };

    FastqConstIterator<numBasesPerByte> findPartialMatch(const FastqConstIterator<numBasesPerByte>& basesBegin,
                                                         const FastqConstIterator<numBasesPerByte>& basesEnd,
                                                         uint32_t                                   minScorePartialMatch) const;

    FastqConstIterator<numBasesPerByte> findFullMatch(const FastqConstIterator<numBasesPerByte>& basesBegin,
                                                      const FastqConstIterator<numBasesPerByte>& basesEnd,
                                                      uint32_t&                                  minScorePartialMatch) const;

    std::vector<uint64_t> getPeq(const std::vector<unsigned char>& adapterCalls) const;

    FastqConstIterator<numBasesPerByte> findLongestExactSuffixMatch(const FastqConstIterator<numBasesPerByte>& basesBegin,
                                                                    const FastqConstIterator<numBasesPerByte>& basesEnd,
                                                                    const std::vector<unsigned char>& adapterSubstring,
                                                                    int tagLength) const;

    void binarizeAdapterCalls(const std::string&          adapter,
                              std::vector<unsigned char>& adapterCalls) const;

    uint32_t                 minimumTrimmedReadLength_;
    uint32_t                 minLengthMatch_;
    uint32_t                 mismatchesToAccountFor_;
    std::vector<uint32_t>    maxAllowedMismatches_;
    std::vector< std::vector<unsigned char> > adapterCallsSubsets_;

    std::vector< std::vector<uint64_t> > peqs_;
    std::vector<uint64_t>                pv0s_;
};

} // namespace conversion
} // namespace bcl2fastq

#include "conversion/AdapterLocator.hpp"

#endif // BCL2FASTQ_CONVERSION_ADAPTER_LOCATOR_HH
