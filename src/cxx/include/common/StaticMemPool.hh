/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file StaticMemPool.hh
 *
 * \brief Declaration of simple static memory pool.
 *
 * \author Marek Balint
 */


#ifndef BCL2FASTQ_COMMON_STATICMEMPOOL_HH
#define BCL2FASTQ_COMMON_STATICMEMPOOL_HH


#include <cstddef>

#include <boost/utility.hpp>


namespace bcl2fastq {


namespace common {


/// \brief Simple static memory pool.
template<std::size_t PageSize, std::size_t PageCount>
class StaticMemPool
{
public:

    /// \brief Size type.
    typedef std::size_t size_type;

public:

    /// \brief Default constructor.
    StaticMemPool();

    /// \brief Copy constructor.
    /// \param that Other instance to copy from.
    StaticMemPool(const StaticMemPool<PageSize, PageCount> &that);

    /// \brief Assignment operator.
    /// \param rhs Right-hand-side parameter.
    StaticMemPool & operator=(const StaticMemPool<PageSize, PageCount> &rhs);

    /// \brief Destructor.
    ~StaticMemPool();

public:

    /// \brief Allocate memory from pool.
    /// \param size Size of memory to be allocated.
    /// \return Pointer to allocated memory or @c NULL in case of error.
    void * allocate(size_type size);

    /// \brief Free previously allocated memory.
    /// \param ptr Pointer to memory previously allocated by this pool instanc.
    void deallocate(void *ptr);

private:

    /// \brief Memory pool.
    char pool_[PageCount][PageSize];

    /// \brief Memory pool metadata.
    /// \note Currently only boolean flag determining whether the given page
    /// is free (@c true) or not (@c false).
    bool metadata_[PageCount];
};


} // namespace common


} // namespace bcl2fastq


#include "common/StaticMemPool.hpp"


#endif // BCL2FASTQ_COMMON_STATICMEMPOOL_HH


