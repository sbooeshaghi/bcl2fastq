/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file BclBaseConversion.hh
 *
 * \brief Declaration of Bcl base conversion functions.
 *
 * \author Aaron Day
 */

#ifndef BCL2FASTQ_CONVERSION_BCLBASECONVERSION_HH
#define BCL2FASTQ_CONVERSION_BCLBASECONVERSION_HH

namespace bcl2fastq
{
namespace conversion
{

static const char bases[256] =
    {'N', 'N', 'N', 'N', 'A', 'C', 'G', 'T',
     'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T',
     'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T',
     'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T',
     'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T',
     'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T',
     'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T',
     'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T',
     'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T',
     'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T',
     'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T',
     'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T',
     'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T',
     'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T',
     'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T',
     'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T',
     'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T',
     'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T',
     'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T',
     'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T',
     'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T',
     'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T',
     'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T',
     'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T',
     'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T',
     'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T',
     'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T',
     'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T',
     'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T',
     'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T',
     'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T',
     'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T'};

static const char noQscoreBases[4] = { 'A', 'C', 'G', 'T' };

static const char reverseComplementBases[256] =
    {'N', 'N', 'N', 'N', 'T', 'G', 'C', 'A',
     'T', 'G', 'C', 'A', 'T', 'G', 'C', 'A',
     'T', 'G', 'C', 'A', 'T', 'G', 'C', 'A',
     'T', 'G', 'C', 'A', 'T', 'G', 'C', 'A',
     'T', 'G', 'C', 'A', 'T', 'G', 'C', 'A',
     'T', 'G', 'C', 'A', 'T', 'G', 'C', 'A',
     'T', 'G', 'C', 'A', 'T', 'G', 'C', 'A',
     'T', 'G', 'C', 'A', 'T', 'G', 'C', 'A',
     'T', 'G', 'C', 'A', 'T', 'G', 'C', 'A',
     'T', 'G', 'C', 'A', 'T', 'G', 'C', 'A',
     'T', 'G', 'C', 'A', 'T', 'G', 'C', 'A',
     'T', 'G', 'C', 'A', 'T', 'G', 'C', 'A',
     'T', 'G', 'C', 'A', 'T', 'G', 'C', 'A',
     'T', 'G', 'C', 'A', 'T', 'G', 'C', 'A',
     'T', 'G', 'C', 'A', 'T', 'G', 'C', 'A',
     'T', 'G', 'C', 'A', 'T', 'G', 'C', 'A',
     'T', 'G', 'C', 'A', 'T', 'G', 'C', 'A',
     'T', 'G', 'C', 'A', 'T', 'G', 'C', 'A',
     'T', 'G', 'C', 'A', 'T', 'G', 'C', 'A',
     'T', 'G', 'C', 'A', 'T', 'G', 'C', 'A',
     'T', 'G', 'C', 'A', 'T', 'G', 'C', 'A',
     'T', 'G', 'C', 'A', 'T', 'G', 'C', 'A',
     'T', 'G', 'C', 'A', 'T', 'G', 'C', 'A',
     'T', 'G', 'C', 'A', 'T', 'G', 'C', 'A',
     'T', 'G', 'C', 'A', 'T', 'G', 'C', 'A',
     'T', 'G', 'C', 'A', 'T', 'G', 'C', 'A',
     'T', 'G', 'C', 'A', 'T', 'G', 'C', 'A',
     'T', 'G', 'C', 'A', 'T', 'G', 'C', 'A',
     'T', 'G', 'C', 'A', 'T', 'G', 'C', 'A',
     'T', 'G', 'C', 'A', 'T', 'G', 'C', 'A',
     'T', 'G', 'C', 'A', 'T', 'G', 'C', 'A',
     'T', 'G', 'C', 'A', 'T', 'G', 'C', 'A'};

static const char noQscoreReverseComplementBases[4] = { 'T', 'G', 'C', 'A' };

/// \brief Convert BCL byte to FASTQ base byte.
/// \param bcl BCL byte to be converted.
/// \return FASTQ base byte.
inline char convertBcl2FastqBase(unsigned char bcl)
{
    return bases[bcl];
}

/// \brief Convert BCL byte to FASTQ base byte.
/// \param bcl BCL byte to be converted.
/// \return FASTQ base byte.
inline char convertBcl2FastqBaseNoQscore(unsigned char bcl)
{
    return noQscoreBases[bcl];
}

/// \brief Convert BCL byte to FASTQ base complement byte.
/// \param bcl BCL byte to be converted.
/// \return FASTQ base complement byte.
inline char convertBcl2FastqBaseComplement(unsigned char bcl)
{
    return reverseComplementBases[bcl];
}

/// \brief Convert BCL byte to FASTQ base complement byte.
/// \param bcl BCL byte to be converted.
/// \return FASTQ base complement byte.
inline char convertBcl2FastqBaseComplementNoQscore(unsigned char bcl)
{
    return noQscoreReverseComplementBases[bcl];
}

/// \brief Convert BCL byte to FASTQ quality byte.
/// \param bcl BCL byte to be converted.
/// \return FASTQ quality byte.
inline char convertBcl2FastqQuality(unsigned char bcl)
{
    const unsigned char quality = bcl >> 2;
    return quality + 33 + (quality==0)*2;
}

/// \brief Convert BCL byte to FASTQ quality byte.
/// \param bcl BCL byte to be converted.
/// \param remappedQscores Map to actual qscore numbers for given qscore bin
/// \return FASTQ quality byte.
inline char convertBcl2FastqQuality2(unsigned char bcl, const std::vector<uint32_t>& remappedQscores)
{
    const unsigned char quality = remappedQscores[bcl >> 2];
    return quality + 33 + (quality==0)*2;
}

/// \brief Convert BCL byte to FASTQ quality byte.
/// \param bcl BCL byte to be converted.
/// \return FASTQ quality byte.
inline char convertBcl2FastqRawQuality(unsigned char bcl)
{
    return bcl >> 2;
}

/// \brief Convert BCL byte to FASTQ quality byte.
/// \param bcl BCL byte to be converted.
/// \param remappedQscores Map to actual qscore numbers for given qscore bin
/// \return FASTQ quality byte.
inline char convertBcl2FastqRawQuality2(unsigned char bcl, const std::vector<uint32_t>& remappedQscores)
{
    return remappedQscores[bcl >> 2];
}

} // namespace bcl2fastq
} // namespace conversion

#endif // BCL2FASTQ_CONVERSION_BCLBASECONVERSION_HH

