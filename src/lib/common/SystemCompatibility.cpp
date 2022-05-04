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

#include <stdio.h>

#include <iostream>
#include <new>

#include <boost/format.hpp>
#include <boost/thread.hpp>

#include "common/Debug.hpp"
#include "common/Exceptions.hpp"
#include "common/SystemCompatibility.hpp"

// TODO: add a proper configuration system as needed
#define HAVE_SIGNAL_H
#define HAVE_MALLOC_H
#define HAVE_TIME_H
#define HAVE_SYS_STAT_H
#define HAVE_SYSCONF
#define HAVE_CLOCK
#define HAVE_STAT

#ifdef HAVE_SIGNAL_H
#include <signal.h>
#else  // #ifdef HAVE_SIGNAL_H
#error Only POSIX systems are supported. The header <signal.h> is required.
#endif  // #ifdef HAVE_SIGNAL_H

/*
 * TODO: separate the implementation into different files, depending on the system
 */
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#else  // #ifdef HAVE_MALLOC_H
#error Only POSIX systems are supported. The header <malloc.h> is required.
#endif  // #ifdef HAVE_MALLOC_H

#ifdef HAVE_TIME_H
#include <time.h>
#else  // #ifdef HAVE_TIME_H
#error The header <time.h> is required.
#endif  // #ifdef HAVE_TIME_H

#ifdef HAVE_SYS_STAT_H
#include <sys/stat.h>
#else  // #ifdef HAVE_SYS_STAT_H
#error The header <sys/stat.h> is required.
#endif  // #ifdef HAVE_SYS_STAT_H

#ifdef HAVE_SYS_IOCTL_H
#include <sys/ioctl.h>
#endif  // #ifdef HAVE_SYS_IOCTL_H

// either linux/falloc.h or fcntl provide fallocate
#ifdef HAVE_LINUX_FALLOC_H
#include <linux/falloc.h>
#endif  // #ifdef HAVE_LINUX_FALLOC_H

// either linux/falloc.h or fcntl provide fallocate
#ifdef HAVE_FCNTL_H
#include <fcntl.h>
#endif  // #ifdef HAVE_FCNTL_H

#include <cassert>

namespace dragenos {
namespace common {

std::string pathStringToStdString(const PathStringType& pathString)
{
  return boost::filesystem::path(pathString).string();
}

int getTerminalWindowSize(unsigned short int& ws_row, unsigned short int& ws_col)
{
#ifdef HAVE_SYS_IOCTL_H
  winsize ws = {0, 0, 0, 0};
  if (!ioctl(STDERR_FILENO, TIOCGWINSZ, &ws)) {
    ws_row = ws.ws_row;
    ws_col = ws.ws_col;
    return 0;
  }
#else   // just return whatever suitable constant dimensions
  static const int VT100_ROWS = 24;
  static const int VT100_COLS = 80;
  ws_row                      = VT100_ROWS;
  ws_col                      = VT100_COLS;
#endif  // HAVE_SYS_IOCTL_H

  return -1;
}

}  // namespace common
}  // namespace dragenos

#ifdef _WIN32

#include <windows.h>

namespace dragenos {
namespace common {

int deleteFile(const PathCharType* filename)
{
  return _wunlink(filename);
}

void truncateFile(const PathCharType* filePath, uint64_t dataSize)
{
  HANDLE file =
      ::CreateFileW(filePath, GENERIC_WRITE, 0, nullptr, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, nullptr);
  if (file != INVALID_HANDLE_VALUE) {
    LARGE_INTEGER newLen;
    newLen.QuadPart = dataSize;

    if (::SetFilePointerEx(file, newLen, nullptr, FILE_BEGIN)) {
      ::SetEndOfFile(file);
    }

    ::CloseHandle(file);
  }
}

unsigned int getMaxOpenFiles()
{
  return _getmaxstdio();
}

int64_t clock()
{
  return ::clock();
}

bool isLittleEndian()
{
  const uint64_t             v = 0x0706050403020100;
  const unsigned char* const p = reinterpret_cast<const unsigned char*>(&v);
  for (unsigned i = 0; i < sizeof(v); ++i) {
    if (p[i] != i) {
      return false;
    }
  }
  return true;
}

void terminateWithCoreDump()
{
  terminate();
}

uint64_t getFileSize(const PathCharType* filePath)
{
  struct __stat64 buffer;
  if (0 != _wstat64(filePath, &buffer)) {
    BOOST_THROW_EXCEPTION(
        common::IoException(errno, std::string("Failed to stat file ") + pathStringToStdString(filePath)));
  }

  return buffer.st_size;
}

boost::filesystem::path getModuleFileName()
{
  const DWORD bufsize = 10240;
  char        buffer[bufsize];
  ::GetModuleFileNameA(nullptr, buffer, bufsize);

  return boost::filesystem::path(buffer);
}

void configureMemoryManagement(const bool disableMultipleArenas, const bool disableFastbins) {}

/**$
 * \brief Exception thrown when there is insufficient resources to perform an operation. For example$
 *        if the adjusting the soft ulimit fails due to a set hard limit$
 */
class Win32Exception : public std::exception, public ExceptionData {
public:
  Win32Exception(DWORD dwLastError, const std::string& message)
    : ExceptionData(
          errno,
          message + ": " + (boost::format("0x%x") % dwLastError).str() + ": " + getLastErrorText(dwLastError))
  {
  }

  static std::string getLastErrorText(DWORD nErrorCode)
  {
    char* msg;
    // Ask Windows to prepare a standard message for a GetLastError() code:
    if (!FormatMessage(
            FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS,
            NULL,
            nErrorCode,
            MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
            (LPSTR)&msg,
            0,
            NULL)) {
      return (boost::format("FormatMessage failed for error code: 0x%x") % nErrorCode).str();
    } else {
      return msg;
    }
  }
};

bool ulimitV(uint64_t availableMemory)
{
  HANDLE jobObject = CreateJobObject(NULL, NULL);
  if (INVALID_HANDLE_VALUE == jobObject) {
    BOOST_THROW_EXCEPTION(Win32Exception(GetLastError(), (boost::format("Failed to CreateJobObject")).str()));
  }

  JOBOBJECT_EXTENDED_LIMIT_INFORMATION extendedLimits = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  extendedLimits.BasicLimitInformation.LimitFlags     = JOB_OBJECT_LIMIT_PROCESS_MEMORY;
  extendedLimits.ProcessMemoryLimit                   = availableMemory;
  BOOL res                                            = SetInformationJobObject(
      jobObject, JobObjectExtendedLimitInformation, &extendedLimits, sizeof(extendedLimits));

  if (!res) {
    CloseHandle(jobObject);
    BOOST_THROW_EXCEPTION(
        Win32Exception(GetLastError(), (boost::format("Failed to SetInformationJobObject")).str()));
  }
  res = AssignProcessToJobObject(jobObject, GetCurrentProcess());
  if (!res) {
    CloseHandle(jobObject);
    BOOST_THROW_EXCEPTION(
        Win32Exception(GetLastError(), (boost::format("AssignProcessToJobObject failed")).str()));
  }

  CloseHandle(jobObject);

  return true;
}

bool ulimitV(uint64_t* pLimit)
{
  *pLimit = 0;
  return true;
}

void hookMalloc(bool (*hook)(size_t size, const void* caller)) {}

unsigned unhookMalloc(bool (*hook)(size_t size, const void* caller))
{
  return 0;
}

int shmget(int key, size_t size, int shmflg)
{
  return 0;
}

void* shmat(int shmid, const void* shmaddr, int shmflg)
{
  return nullptr;
}

int shmctl(int shmid, int cmd, struct shmid_ds* buf)
{
  return 0;
}

int shmdt(const void* shmaddr)
{
  return 0;
}

}  // namespace common
}  // namespace dragenos

#else

#include <sys/resource.h>

namespace dragenos {
namespace common {

int deleteFile(const PathCharType* filename)
{
  return unlink(filename);
}

void truncateFile(const char* filePath, uint64_t dataSize)
{
  if (-1 == truncate(filePath, dataSize)) {
    BOOST_THROW_EXCEPTION(IoException(
        errno, (boost::format("Failed to truncate file '*s' to *i") % filePath % dataSize).str()));
  }
}

unsigned int getMaxOpenFiles()
{
#ifdef HAVE_SYSCONF
  assert(0 < sysconf(_SC_OPEN_MAX));
  return sysconf(_SC_OPEN_MAX);
#else
#error 'sysconf' is required
#endif
}

int64_t clock()
{
#ifdef HAVE_CLOCK
  return ::clock();
#else
#error 'clock' is required (from <time.h>)
#endif
}

bool isLittleEndian()
{
  const uint64_t v = 0x0706050403020100;
  const unsigned char* const p = reinterpret_cast<const unsigned char*>(&v);
  for (unsigned i = 0; i < sizeof(v); ++i) {
    if (p[i] != i) {
      return false;
    }
  }
  return true;
}

void terminateWithCoreDump()
{
  raise(SIGSEGV);
}

uint64_t getFileSize(const PathCharType* filePath)
{
#ifdef HAVE_STAT
  struct stat s;
  if (0 != stat(filePath, &s)) {
    BOOST_THROW_EXCEPTION(common::IoException(errno, std::string("Failed to stat file ") + filePath));
  }
  return s.st_size;
#else
#error 'stat' is required
#endif
}

boost::filesystem::path getModuleFileName()
{
  char szBuffer[10240];
  int readBytes = readlink("/proc/self/exe", szBuffer, sizeof(szBuffer));
  DRAGEN_OS_ASSERT_MSG(-1 != readBytes, "TODO: handle the readlink error: " << errno);
  // readlink does not zero-terminate the string.
  szBuffer[readBytes] = 0;
  return boost::filesystem::path(szBuffer);
}

}  // namespace common
}  // namespace dragenos

#ifndef DRAGEN_OS_CYGWIN

namespace dragenos {
namespace common {
void configureMemoryManagement(const bool disableMultipleArenas, const bool disableFastbins)
{
  if (disableMultipleArenas) {
    // By default linux creates  ((NUMBER_OF_CPU_CORES) * (sizeof(int64_t) == 4 ? 2 : 8)) arenas to allow for
    // (I guess) faster concurrent memory management. As iSAAC pre-allocates all the memory and spends
    // majority of time doing something other than memory allocation/deallocation, having more than one arena
    // is purely a waste of virtual memory. Virtual memory is important for iSAAC as it normally tries to use
    // all of the virtual memory that is allowed to the process. So, this trick should save in the range of 2
    // gigabytes.
    // TODO: Interestingly enough it actually saves more and the value seems to depend on the number of
    // threads. Would be nice to figure out why...
    mallopt(M_ARENA_MAX, 1);
  }

  if (disableFastbins) {
    // in Linux, fastbins are enabled by default in order to improve performance of
    // applications that allocate and free small chunks of memory on multiple threads often.
    // This causes unnecessary memory overhead and fragmentation which easily amounts to
    // a loss of 10-20 gigabytes of RAM in 400-plex bam generation.
    // As dragenos, where it matters, allocates memory up-front and controls concurrency,
    // fastbins are disabled.
    mallopt(M_MXFAST, 0);
  }
}

bool ulimitV(uint64_t availableMemory)
{
  const rlimit rl = {availableMemory, availableMemory};
  if (setrlimit(RLIMIT_AS, &rl)) {
    BOOST_THROW_EXCEPTION(dragenos::common::ResourceException(
        errno,
        (boost::format("Failed to set the memory consumption limit to: %d bytes") % availableMemory).str()));
  }
  return true;
}

bool ulimitV(uint64_t* pLimit)
{
  rlimit rl = {0, 0};
  if (-1 == getrlimit(RLIMIT_AS, &rl)) {
    return false;
  }
  *pLimit = RLIM_INFINITY == rl.rlim_cur ? -1 : rl.rlim_cur;
  return true;
}

#if 0

static boost::mutex block_malloc_hook_mutex_;

static void* (*old_malloc_hook_)(size_t, const void*) = 0;
static bool (*user_hook_)(size_t size, const void *caller) = 0;

unsigned mallocCount_(0);
static void * malloc_hook(size_t size, const void *caller)
{
    boost::unique_lock<boost::mutex> lock(block_malloc_hook_mutex_);
    ++mallocCount_;

    __malloc_hook = old_malloc_hook_;

    assert(user_hook_);
    if (!user_hook_(size, caller))
    {
        terminateWithCoreDump();
        // in case termination did not terminate...
        std::cerr << "ERROR: blocked allocation of " << size << " bytes. Returning 0 \n";
        __malloc_hook = malloc_hook;
        return 0;
    }
    else
    {
        void *result  = malloc(size);
        __malloc_hook = malloc_hook;
        return result;
    }
}

void hookMalloc(bool (*hook)(size_t size, const void *caller))
{
    boost::unique_lock<boost::mutex> lock(block_malloc_hook_mutex_);

    assert(__malloc_hook != malloc_hook);
    assert(!user_hook_);
    old_malloc_hook_ = __malloc_hook;
    user_hook_ = hook;
    mallocCount_ = 0;
    __malloc_hook = malloc_hook;

    //std::cerr << "malloc blocked\n";
}

unsigned unhookMalloc(bool (*hook)(size_t size, const void *caller))
{
    boost::unique_lock<boost::mutex> lock(block_malloc_hook_mutex_);

    assert(__malloc_hook == malloc_hook);
    assert(user_hook_ == hook);

    __malloc_hook = old_malloc_hook_;
    user_hook_ = 0;

//    std::cerr << "malloc unblocked\n";
    return mallocCount_;
}

#endif  // #if 0

}  // namespace common
}  // namespace dragenos

#else  // DRAGEN_OS_CYGWIN

#include <w32api/windows.h>

namespace dragenos {
namespace common {

void configureMemoryManagement(const bool disableMultipleArenas, const bool disableFastbins)
{
  if (disableMultipleArenas) {
    // By default linux creates  ((NUMBER_OF_CPU_CORES) * (sizeof(int64_t) == 4 ? 2 : 8)) arenas to allow for
    // (I guess) faster concurrent memory management. As iSAAC pre-allocates all the memory and spends
    // majority of time doing something other than memory allocation/deallocation, having more than one arene
    // is purely a waste of virtual memory. Virtual memory is important for iSAAC as it normally tries to use
    // all of the virtual memory that is allowed to the process. So, this trick should save in the range of 2
    // gigabytes.
    // TODO: Interestingly enough it actually saves more and the value seems to depend on the number of
    // threads. Would be nice to figure out why...
    //mallopt(M_ARENA_MAX, 1);
    DRAGEN_OS_THREAD_CERR << "WARNING: M_ARENA_MAX is not available in CYGWIN." << std::endl;
  }

  if (disableFastbins) {
    // in Linux, fastbins are enabled by default in order to improve performance of
    // applications that allocate and free small chunks of memory on multiple threads often.
    // This causes unnecessary memory overhead and fragmentation which easily amounts to
    // a loss of 10-20 gigabytes of RAM in 400-plex bam generation.
    // As dragenos, where it matters, allocates memory up-front and controls concurrency,
    // fastbins are disabled.
    mallopt(M_MXFAST, 0);
  }
}

/**$
 * \brief Exception thrown when there is insufficient resources to perform an operation. For example$
 *        if the adjusting the soft ulimit fails due to a set hard limit$
 */
class Win32Exception : public std::exception, public ExceptionData {
public:
  Win32Exception(DWORD dwLastError, const std::string& message)
    : ExceptionData(
          errno,
          message + ": " + (boost::format("0x%x") % dwLastError).str() + ": " + getLastErrorText(dwLastError))
  {
  }

  static std::string getLastErrorText(DWORD nErrorCode)
  {
    char* msg;
    // Ask Windows to prepare a standard message for a GetLastError() code:
    if (!FormatMessage(
            FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS,
            NULL,
            nErrorCode,
            MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
            (LPSTR)&msg,
            0,
            NULL)) {
      return (boost::format("FormatMessage failed for error code: 0x%x") % nErrorCode).str();
    } else {
      return msg;
    }
  }
};

bool ulimitV(uint64_t availableMemory)
{
  HANDLE hJobObject = CreateJobObject(NULL, NULL);
  if (INVALID_HANDLE_VALUE == hJobObject) {
    BOOST_THROW_EXCEPTION(Win32Exception(GetLastError(), (boost::format("Failed to CreateJobObject")).str()));
  }

  JOBOBJECT_EXTENDED_LIMIT_INFORMATION extendedLimits = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  extendedLimits.BasicLimitInformation.LimitFlags     = JOB_OBJECT_LIMIT_PROCESS_MEMORY;
  extendedLimits.ProcessMemoryLimit                   = availableMemory;
  BOOL res                                            = SetInformationJobObject(
      hJobObject, JobObjectExtendedLimitInformation, &extendedLimits, sizeof(extendedLimits));

  if (!res) {
    CloseHandle(hJobObject);
    BOOST_THROW_EXCEPTION(
        Win32Exception(GetLastError(), (boost::format("Failed to SetInformationJobObject")).str()));
  }
  res = AssignProcessToJobObject(hJobObject, GetCurrentProcess());
  if (!res) {
    CloseHandle(hJobObject);
    BOOST_THROW_EXCEPTION(
        Win32Exception(GetLastError(), (boost::format("AssignProcessToJobObject failed")).str()));
  }

  CloseHandle(hJobObject);

  return true;
}

bool ulimitV(uint64_t* pLimit)
{
  *pLimit = 0;
  return true;
}

void hookMalloc(bool (*hook)(size_t size, const void* caller))
{
  // memory control is not supported under cygwin
}

unsigned unhookMalloc(bool (*hook)(size_t size, const void* caller))
{
  return 0;
}

}  // namespace common
}  // namespace dragenos

#endif  // DRAGEN_OS_CYGWIN

#endif  // #ifdef _WIN32
