/**
 ** DRAGEN Open Source Software
 ** Copyright (c) 2019-2020 Illumina, Inc.
 ** All rights reserved.
 **
 ** This software is provided under the terms and conditions of the
 ** GNU GENERAL PUBLIC LICENSE Version 3
 **
 ** You should have received a copy of the GNU GENERAL PUBLIC LICENSE Version 3
 ** along with this program. If not, see
 ** <https://github.com/illumina/licenses/>.
 **
 **/

#ifndef COMMON_SYSTEM_COMPATIBILITY_HPP
#define COMMON_SYSTEM_COMPATIBILITY_HPP

#include <boost/filesystem.hpp>

#if defined(_WIN32) && !defined(DRAGEN_OS_CYGWIN)

#define THREAD_LOCAL thread_local
#define FILE_THAT_ALWAYS_EXISTS "nul"
#define posix_fadvise(FILE, X, Y, Z) false
#define TSTRING(s) L##s
// 2048 is maximum you can have on win32. And it is not a default
#define SET_MAX_FILES _setmaxstdio(2048)
#define DRAGEN_OS_TRACE_STAT(prefix)

#define __builtin_popcountll _popcnt64

// vectorizing 128 bit unsigned integer
typedef struct __declspec(intrin_type) __declspec(align(16)) vint128_t {
  uint64_t   m128i_u64[2];
  vint128_t& operator|=(const vint128_t& that)
  {
    m128i_u64[0] |= that.m128i_u64[0];
    m128i_u64[1] |= that.m128i_u64[1];
    return *this;
  }

  vint128_t operator~() const
  {
    vint128_t ret = {~m128i_u64[0], ~m128i_u64[1]};
    return ret;
  }

  vint128_t operator&(const vint128_t& that) const
  {
    vint128_t ret = {m128i_u64[0] & that.m128i_u64[0], m128i_u64[1] & that.m128i_u64[1]};
    return ret;
  }

} vint128_t;

#else

#define THREAD_LOCAL __thread
#define FILE_THAT_ALWAYS_EXISTS "/dev/null"
#define TSTRING(s) s
// nothing to set here. users do it with ulimit
#define SET_MAX_FILES

// vectorizing 128 bit unsigned integer
typedef unsigned long long vint128_t __attribute__((__vector_size__(16)));

#endif

namespace dragenos {
namespace common {
using PathCharType   = boost::filesystem::path::value_type;
using PathStringType = boost::filesystem::path::string_type;

std::string pathStringToStdString(const PathStringType& pathString);

int deleteFile(const PathCharType* filename);

int getTerminalWindowSize(unsigned short int& ws_row, unsigned short int& ws_col);

void truncateFile(const PathCharType* filePath, uint64_t dataSize);

/// Maximum number of files that a process can have opened at the same time
unsigned int getMaxOpenFiles();

/// File size in bytes as returned by stat
uint64_t getFileSize(const PathCharType* filePath);

/// Determine the processor time
int64_t clock();

/// Check if the architecture is little endian
bool isLittleEndian();

/// limit virtual memory size available to process. (equivalent of ulimit -v)
bool ulimitV(const uint64_t availableMemory);
/// retrieves the current ulimit -v
bool ulimitV(uint64_t* pLimit);

/**
 * \brief Sets a hook that monitors memory allocations.
 *
 * \param hook  Hook to set. If hook returns false, allocation will attempt to terminateWithCoreDump and fail.
 *              Calls to hook happen under the boost::unique_lock<boost::mutex>
 */
//void hookMalloc(bool (*hook)(size_t size, const void *caller));

/**
 * \brief Removes block set by blockMalloc. Returns the number of allocations made since last blockMalloc
 * call.
 *
 * \param hook  Used to ensure that the hook removed is the one that was previously set.
 */
//unsigned unhookMalloc(bool (*hook)(size_t size, const void *caller));

/**
 * \brief Generate a core dump with a meaningful backtrace
 */
void terminateWithCoreDump();

/**
 * \brief Disables memory management optimizations that are detrimental to the access pattern used in
 * high-performance parts of the product
 */
void configureMemoryManagement(const bool disableMultipleArenas, const bool disableFastbins);

boost::filesystem::path getModuleFileName();

/**
 * \brief calls linux-specific fallocate with FALLOC_FL_KEEP_SIZE to pre-allocate file on disk in a way that
 * does not mess up reopening for append. Notice that posix_fallocate does not do the job.
 */
int linuxFallocate(int fd, std::size_t offset, std::size_t len);

/**
 * \brief calls linux-specific ftruncate. Only compiles if fallocate is available
 */
int linuxFtruncate(int fd, std::size_t len);

}  // namespace common
}  // namespace dragenos

#endif  // #ifndef COMMON_SYSTEM_COMPATIBILITY_HPP
