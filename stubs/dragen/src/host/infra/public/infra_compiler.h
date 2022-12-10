// Copyright 2018 Edico Genome Corporation. All rights reserved.
// Copyright 2018-2020 Illumina, Inc. All rights reserved.

//
// This file contains confidential and proprietary information of the Edico
// Genome Corporation and is protected under the U.S. and international
// copyright and other intellectual property laws.
//
// $Id$
// $Author$
// $Change$
// $DateTime$
//

#pragma once

#include <stdint.h>

#define POWEROF2(a) (((a) & ((a)-1)) == 0)

//----------------------------------------------------------------------
// GCC
//----------------------------------------------------------------------
#if defined(__GNUC__)
#define __LIKELY(x) __builtin_expect((x), 1)
#define __UNLIKELY(x) __builtin_expect((x), 0)
#define __PREFETCH(addr, rw, loc) __builtin_prefetch(addr, rw, loc)

#define __CTZ(x) __builtin_ctz(x)
#define __CTZL(x) __builtin_ctzl(x)
#define __CTZLL(x) __builtin_ctzll(x)

#define __CLZ(x) __builtin_clz(x)
#define __CLZL(x) __builtin_clzl(x)
#define __CLZLL(x) __builtin_clzll(x)

#define __FFS(x) __builtin_ffs(x)
#define __FFSL(x) __builtin_ffsl(x)
#define __FFSLL(x) __builtin_ffsll(x)

#define __PARITY(x) __builtin_parity(x)
#define __PARITYL(x) __builtin_parityl(x)
#define __PARITYLL(x) __builtin_parityll(x)

#define __POPCOUNT(x) __builtin_popcount(x)
#define __POPCOUNTL(x) __builtin_popcountl(x)
#define __POPCOUNTLL(x) __builtin_popcountll(x)

#define __OFFSETOF(type, field) __builtin_offsetof(type, field)
#define __SIZEOF(type, field) sizeof(((type *)0)->field)

#define __BSWAP32(x) __builtin_bswap32(x)
#define __BSWAP64(x) __builtin_bswap64(x)

//----------------------------------------------------------------------
// x86
#if defined(_TARGET_X86_)
#define mb() asm volatile("mfence" ::: "memory")
#define rmb() asm volatile("lfence" ::: "memory")
#define read_barrier_depends()                                                 \
  do {                                                                         \
  } while (0)
#define wmb() asm volatile("sfence" ::: "memory")
#define compiler_barrier()                                                     \
  do {                                                                         \
    asm volatile("" : : : "memory");                                           \
  } while (0)
// write/read/generic barrier between vcpus/lcores
#define smp_wmb() compiler_barrier()
#define smp_rmb() compiler_barrier()
#define smp_mb()                                                               \
  do {                                                                         \
    asm volatile("lock addl $0, -128(%%rsp); " ::: "memory");                  \
  } while (0)

#define CACHE_LINE_SIZE 64
#define CACHE_LINE_SIZE_LOG2 6

//----------------------------------------------------------------------
// powerpc
#elif defined(_TARGET_PPC_)
#define mb() asm volatile("sync" ::: "memory")
#define rmb() asm volatile("sync" ::: "memory")
#define read_barrier_depends()                                                 \
  do {                                                                         \
  } while (0)
#define wmb() asm volatile("sync" ::: "memory")
// write/read/generic barrier between vcpus/lcores
#define smp_wmb() wmb()
#define smp_rmb() rmb()
#define smp_mb() mb()

#define CACHE_LINE_SIZE 128
#define CACHE_LINE_SIZE_LOG2 7

#elif defined(_TARGET_WIN_)
#warning porting needed for Windows
#define mb()
#define rmb()
#define read_barrier_depends()
#define wmb()

//----------------------------------------------------------------------
// ARM
#elif defined(_TARGET_ARM_)
#define isb() asm volatile("isb" ::: "memory")
#define dmb(opt) asm volatile("dmb " #opt ::: "memory")
#define read_barrier_depends()                                                 \
  do {                                                                         \
  } while (0)
#define dsb(opt) asm volatile("dsb " #opt ::: "memory")

#define mb() dsb(sy)
#define rmb() dsb(ld)
#define wmb() dsb(st)
#define dma_rmb() dmb(oshld)
#define dma_wmb() dmb(oshst)
// write/read/generic barrier between vcpus/lcores
#define smp_wmb() dmb(ishst)
#define smp_rmb() dmb(ishld)
#define smp_mb() dmb(ish)

#define CACHE_LINE_SIZE 64
#define CACHE_LINE_SIZE_LOG2 6

#endif

#if POWEROF2(CACHE_LINE_SIZE)
// ok
#else
#error CACHE_LINE_SIZE must be a power of 2
#endif

// macros for pointer and size alignement
#define __ALIGNED(a) __attribute__((__aligned__(a)))
#define __CACHE_ALIGNED __ALIGNED(CACHE_LINE_SIZE)
#define __PTR_ADD(ptr, x) ((void *)((uintptr_t)(ptr) + (x)))

// floor + ceil for sizes
#define __ALIGN_CACHE_SIZE_FLOOR(val)                                          \
  (typeof(val))((val) & (~((typeof(val))(CACHE_LINE_SIZE - 1))))
#define __ALIGN_CACHE_SIZE_CEIL(val)                                           \
  __ALIGN_CACHE_SIZE_FLOOR((val) + ((typeof(val))(CACHE_LINE_SIZE - 1)))
// floor + ceil for pointers
#define __ALIGNPTR_CACHE_SIZE_FLOOR(ptr)                                       \
  (typeof(ptr))(__ALIGN_CACHE_SIZE_FLOOR((uintptr_t)ptr))
#define __ALIGNPTR_CACHE_SIZE_CEIL(ptr)                                        \
  __PTRALIGN_CACHE_SIZE_FLOOR((typeof(ptr))__PTR_ADD(ptr, CACHE_LINE_SIZE - 1))

static inline int atomic16_cmpset(uint16_t volatile *ptr, const uint16_t exp,
                                  const uint16_t val) {
  return __sync_bool_compare_and_swap(ptr, exp, val);
}
static inline uint32_t atomic16_xchg(uint16_t volatile *ptr,
                                     const uint16_t val) {
  return __atomic_exchange_n(ptr, val, __ATOMIC_SEQ_CST);
}
static inline int32_t atomic16_read(int16_t volatile *ptr) { return *ptr; }
static inline void atomic16_write(int16_t volatile *ptr, const int16_t val) {
  *ptr = val;
}
static inline void atomic16_add(int16_t volatile *ptr, const int16_t val) {
  __sync_fetch_and_add(ptr, val);
}
static inline void atomic16_sub(int16_t volatile *ptr, const int16_t val) {
  __sync_fetch_and_add(ptr, val);
}
static inline void atomic16_inc(int16_t volatile *ptr) { atomic16_add(ptr, 1); }
static inline void atomic16_dec(int16_t volatile *ptr) { atomic16_sub(ptr, 1); }
static inline int32_t atomic16_addret(int16_t volatile *ptr,
                                      const int16_t val) {
  return __sync_add_and_fetch(ptr, val);
}
static inline int32_t atomic16_subret(int16_t volatile *ptr,
                                      const int16_t val) {
  return __sync_sub_and_fetch(ptr, val);
}
static inline int atomic16_incntst(int16_t volatile *ptr) {
  return __sync_add_and_fetch(ptr, 1) == 0;
}
static inline int atomic16_decntst(int16_t volatile *ptr) {
  return __sync_sub_and_fetch(ptr, 1) == 0;
}
static inline int atomic16_tstnset(uint16_t volatile *ptr) {
  return atomic16_cmpset(ptr, 0, 1);
}

static inline int atomic32_cmpset(uint32_t volatile *ptr, const uint32_t exp,
                                  const uint32_t val) {
  return __sync_bool_compare_and_swap(ptr, exp, val);
}
static inline uint32_t atomic32_xchg(uint32_t volatile *ptr,
                                     const uint32_t val) {
  return __atomic_exchange_n(ptr, val, __ATOMIC_SEQ_CST);
}
static inline int32_t atomic32_read(int32_t volatile *ptr) { return *ptr; }
static inline void atomic32_write(int32_t volatile *ptr, const int32_t val) {
  *ptr = val;
}
static inline void atomic32_add(int32_t volatile *ptr, const int32_t val) {
  __sync_fetch_and_add(ptr, val);
}
static inline void atomic32_sub(int32_t volatile *ptr, const int32_t val) {
  __sync_fetch_and_add(ptr, val);
}
static inline void atomic32_inc(int32_t volatile *ptr) { atomic32_add(ptr, 1); }
static inline void atomic32_dec(int32_t volatile *ptr) { atomic32_sub(ptr, 1); }
static inline int32_t atomic32_addret(int32_t volatile *ptr,
                                      const int32_t val) {
  return __sync_add_and_fetch(ptr, val);
}
static inline int32_t atomic32_subret(int32_t volatile *ptr,
                                      const int32_t val) {
  return __sync_sub_and_fetch(ptr, val);
}
static inline int atomic32_incntst(int32_t volatile *ptr) {
  return __sync_add_and_fetch(ptr, 1) == 0;
}
static inline int atomic32_decntst(int32_t volatile *ptr) {
  return __sync_sub_and_fetch(ptr, 1) == 0;
}
static inline int atomic32_tstnset(uint32_t volatile *ptr) {
  return atomic32_cmpset(ptr, 0, 1);
}

static inline int atomic64_cmpset(uint64_t volatile *ptr, const uint64_t exp,
                                  const uint32_t val) {
  return __sync_bool_compare_and_swap(ptr, exp, val);
}
static inline uint64_t atomic64_xchg(uint64_t volatile *ptr,
                                     const uint64_t val) {
  return __atomic_exchange_n(ptr, val, __ATOMIC_SEQ_CST);
}
static inline int64_t atomic64_read(int64_t volatile *ptr) { return *ptr; }
static inline void atomic64_write(int64_t volatile *ptr, const int64_t val) {
  *ptr = val;
}
static inline void atomic64_add(int64_t volatile *ptr, const int64_t val) {
  __sync_fetch_and_add(ptr, val);
}
static inline void atomic64_sub(int64_t volatile *ptr, const int64_t val) {
  __sync_fetch_and_add(ptr, val);
}
static inline void atomic64_inc(int64_t volatile *ptr) { atomic64_add(ptr, 1); }
static inline void atomic64_dec(int64_t volatile *ptr) { atomic64_sub(ptr, 1); }
static inline int64_t atomic64_addret(int64_t volatile *ptr,
                                      const int64_t val) {
  return __sync_add_and_fetch(ptr, val);
}
static inline int64_t atomic64_subret(int64_t volatile *ptr,
                                      const int64_t val) {
  return __sync_sub_and_fetch(ptr, val);
}
static inline int atomic64_incntst(int64_t volatile *ptr) {
  return __sync_add_and_fetch(ptr, 1) == 0;
}
static inline int atomic64_decntst(int64_t volatile *ptr) {
  return __sync_sub_and_fetch(ptr, 1) == 0;
}
static inline int atomic64_tstnset(uint64_t volatile *ptr) {
  return atomic64_cmpset(ptr, 0, 1);
}

//----------------------------------------------------------------------
// MS
//----------------------------------------------------------------------
#elif defined(MSVC)
#error MS toolchain not supported

//----------------------------------------------------------------------
// ARM
//----------------------------------------------------------------------
#elif defined(__ARM__)
#error ARM toolchain not supported

#else
#error Toolchain not supported
#endif
