/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file StaticMemPool.hpp
 *
 * \brief Implementation of simple static memory pool.
 *
 * \author Marek Balint
 */


#ifndef BCL2FASTQ_COMMON_STATICMEMPOOL_HPP
#define BCL2FASTQ_COMMON_STATICMEMPOOL_HPP


#include <algorithm>

#include <cstring>

#include "common/Debug.hh"
#include "common/Logger.hh"
#include "common/StaticMemPool.hh"


namespace bcl2fastq {


namespace common {


template<std::size_t PageSize, std::size_t PageCount>
StaticMemPool<PageSize, PageCount>::StaticMemPool()
: pool_()
, metadata_()
{
    std::fill(&metadata_[0], &metadata_[PageCount], true);
}

template<std::size_t PageSize, std::size_t PageCount>
StaticMemPool<PageSize, PageCount>::StaticMemPool(const StaticMemPool<PageSize, PageCount> &that)
: pool_()
, metadata_()
{
    std::copy(&that.pool_[0][0], &that.pool[PageCount-1][PageSize], &pool_[0][0]);
    std::copy(&that.metadata_[0], &that.metadata_[PageCount], &metadata_[0]);
}

template<std::size_t PageSize, std::size_t PageCount>
StaticMemPool<PageSize, PageCount> & StaticMemPool<PageSize, PageCount>::operator=(const StaticMemPool<PageSize, PageCount> &rhs)
{
    if (this != &rhs)
    {
        std::copy(&rhs.pool_[0][0], &rhs.pool[PageCount-1][PageSize], &pool_[0][0]);
        std::copy(&rhs.metadata_[0], &rhs.metadata_[PageCount], &metadata_[0]);
    }
    return *this;
}

template<std::size_t PageSize, std::size_t PageCount>
StaticMemPool<PageSize, PageCount>::~StaticMemPool()
{
}

template<std::size_t PageSize, std::size_t PageCount>
void * StaticMemPool<PageSize, PageCount>::allocate(typename StaticMemPool<PageSize, PageCount>::size_type size)
{
    if (size > PageSize)
    {
        return NULL;
    }

    bool *freePage = std::find(&metadata_[0], &metadata_[PageCount], true);
    if (freePage == &metadata_[PageCount])
    {
        return NULL;
    }

    size_type freePageIdx = freePage - &metadata_[0];
    void *ret = static_cast<void *>(&pool_[freePageIdx][0]);
    metadata_[freePageIdx] = false;

    return ret;
}

template<std::size_t PageSize, std::size_t PageCount>
void StaticMemPool<PageSize, PageCount>::deallocate(void *ptr)
{
    char *char_ptr = static_cast<char *>(ptr);
    BCL2FASTQ_ASSERT_MSG(char_ptr >= &pool_[0][0], "Invalid pointer: " << ptr);
    BCL2FASTQ_ASSERT_MSG(char_ptr < &pool_[PageCount][0], "Invalid pointer: " << ptr);
    BCL2FASTQ_ASSERT_MSG((char_ptr - &pool_[0][0]) % PageSize == 0, "Invalid pointer: " << ptr);

    size_type pageIdx = (char_ptr - &pool_[0][0]) / PageSize;
    if (metadata_[pageIdx] == true)
    {
        BCL2FASTQ_LOG(LogLevel::DEBUG) << "Memory already free: " << ptr << std::endl;
        return;
    }
    metadata_[pageIdx] = true;
}


} // namespace common


} // namespace bcl2fastq


#endif // BCL2FASTQ_COMMON_STATICMEMPOOL_HPP




