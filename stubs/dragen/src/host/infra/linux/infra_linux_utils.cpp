// Copyright 2015-2018 Edico Genome Corporation. All rights reserved.
// Copyright 2018-2020 Illumina, Inc. All rights reserved.
//
// This file contains confidential and proprietary information of the Edico Genome
// Corporation and is protected under the U.S. and international copyright and other
// intellectual property laws.
//
// $Id$
// $Author$
// $Change$
// $DateTime$
//

#include <fcntl.h>
#include <pwd.h>
#include <signal.h>
#include <sys/sysinfo.h>
#include <sys/types.h>
#include <sys/utsname.h>
#include <unistd.h>

#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <cstring>
#include <fstream>
#include <iostream>
#include <mutex>
#include <sstream>
#include <thread>

//#include "infra_assert.hpp"
//#include "infra_filesystem_utils.hpp"
//#include "infra_linux_utils.hpp"
//#include "infra_string_utils.hpp"

namespace infra {

static std::string g_kernelVersionStr;
static int         g_kernelVersion;

// static constexpr int MKVER(const int major, const int minor, const int patch)
// {
//   return major * 10000 + minor * 100 + patch;
// }

// static constexpr int MKVER(const int major, const int minor)
// {
//   return major * 10000 + minor * 100;
// }

//------------------------------------------------------------------------alain
int GetDmiValue(const std::string& label, std::string& value)
{
  value.clear();

  FILE* dmiOutput = popen("sudo /usr/sbin/dmidecode -t 2", "r");
  if (dmiOutput == NULL) {
    perror("dmidecode popen");
    pclose(dmiOutput);
    return -1;
  }

  // get a line and process
  char*  line     = NULL;
  size_t linesize = 0;
  while (getline(&line, &linesize, dmiOutput) > 0) {
    if (linesize > 0) {
      std::string s(line);
      std::transform(s.begin(), s.end(), s.begin(), ::tolower);

      // search for label
      if (std::strstr(s.c_str(), label.c_str()) != NULL) {
        // label found
        if (!value.empty()) {
          // duplicate value
          return -2;
        }
        // search for value after ':'
        const char* ptr = std::strchr(line, ':');
        if (ptr != NULL) {
          // skip ':'
          ++ptr;
          // skip heading spaces
          while (*ptr && isspace(*ptr)) {
            ptr++;
          }
          // assign value
          value.assign(ptr);
          boost::algorithm::trim(value);
        }
      }
    }

    free(line);
    line     = NULL;
    linesize = 0;
  }

  pclose(dmiOutput);
  return value.empty() ? -3 : 0;
}

//------------------------------------------------------------------------alain
std::string getExecutablePath()
{
  char buf[PATH_MAX] = {'\0'};

  ssize_t len = ::readlink("/proc/self/exe", buf, sizeof(buf));
  if (len != -1) {
    buf[len]  = '\0';
    char* ptr = std::strrchr(buf, '/');
    if (ptr != NULL) {
      ptr[0] = '\0';
    }
  } else {
    buf[0] = '\0';
  }
  return std::string(buf);
}

//------------------------------------------------------------------------alain
std::string getExecutableName()
{
  char buf[PATH_MAX] = {'\0'};

  ssize_t len = ::readlink("/proc/self/exe", buf, sizeof(buf));
  if (len != -1) {
    buf[len] = '\0';
  }
  char* ptr = std::strrchr(buf, '/');
  if (ptr) {
    return std::string(ptr + 1);
  }
  return std::string(buf);
}

//------------------------------------------------------------------------alain
std::string getExecutableDirectory()
{
  boost::filesystem::path exedir(getExecutablePath());
  return exedir.string();
}

//------------------------------------------------------------------------alain
uint32_t GetNumSystemThreads()
{
  return std::thread::hardware_concurrency();
}

//------------------------------------------------------------------------alain
std::string GetUserId()
{
  // determine username
  struct passwd* pw;
  uid_t          uid;

  uid = geteuid();
  pw  = getpwuid(uid);
  if (pw) {
    return std::string(pw->pw_name);
  }

  std::stringstream ss;
  ss << "none(" << std::dec << uid << ")";
  return std::string(ss.str());
}

//------------------------------------------------------------------------alain
// to get the real userid when running sudo
std::string GetRealUserId()
{
  auto ptr = std::getenv("SUDO_USER");
  if (ptr) {
    if (*ptr != '\0') {
      return std::string(ptr);
    } else {
      ptr = std::getenv("LOGNAME");
      if (ptr and *ptr != '\0') {
        return std::string(ptr);
      }
    }
  }
  // fallback on default
  return GetUserId();
}

//------------------------------------------------------------------------alain
uint64_t GetSystemMemorySize()
{
  struct sysinfo si;
  sysinfo(&si);
  return si.totalram;
}

//------------------------------------------------------------------------alain
uint32_t GetCacheLineSize()
{
  const char*   path = "/sys/devices/system/cpu/cpu0/cache/index0/coherency_line_size";
  std::ifstream ifs(path, std::ios::in);
  if (ifs.is_open()) {
    uint32_t size{};
    ifs >> std::dec >> size;
    if (ifs.good()) {
      return size;
    }
    // unable to read
    return 0u;
  }

  // unable to open file
  return 0u;
}

//------------------------------------------------------------------------alain
std::string GetKernelName()
{
#if 0
  return std::string("Linux");
#else
  struct utsname unameData;
  int            result = uname(&unameData);

  if (result != 0) {
    return std::string("?");
  }

  return std::string(unameData.sysname);
#endif
}

//------------------------------------------------------------------------alain
std::string GetKernelArch()
{
  struct utsname unameData;
  int            result = uname(&unameData);

  if (result != 0) {
    return std::string("?");
  }

  return std::string(unameData.machine);
}

////------------------------------------------------------------------------alain
//std::string GetKernelVersionFull()
//{
//  static std::once_flag g_fullVersionOnceFlag;
//
//  auto rdfunc = []() {
//    struct utsname unameData;
//    int            result = uname(&unameData);
//
//    if (result != 0) {
//      return std::string("?");
//    }
//
//    g_kernelVersionStr = std::string(unameData.release);
//
//    int major{};
//    int minor{};
//    int patch{};
//    infra::StringUtils::DecodeSwVersion(g_kernelVersionStr, major, minor, patch);
//    g_kernelVersion = MKVER(major, minor, patch);
//
//    return g_kernelVersionStr;
//  };
//
//  std::call_once(g_fullVersionOnceFlag, rdfunc);
//  return g_kernelVersionStr;
//}
//
////------------------------------------------------------------------------alain
//// major+minor
//std::string GetKernelVersion()
//{
//  auto       vfull = GetKernelVersionFull();
//  const auto pos1  = vfull.find('.');
//  if (pos1 == std::string::npos) {
//    return vfull;
//  }
//
//  const auto pos2 = vfull.find('.', pos1 + 1);
//  if (pos2 == std::string::npos) {
//    return vfull;
//  }
//  return vfull.substr(0, pos2);
//}
//
////------------------------------------------------------------------------alain
//bool KernelVersionAtLeast(const int major, const int minor, const int patch)
//{
//  (void)GetKernelVersionFull();
//  return g_kernelVersion >= MKVER(major, minor, patch);
//}

//------------------------------------------------------------------------alain
int GetKernelRandomizeAddressStatus()
{
  std::ifstream ifs("/proc/sys/kernel/randomize_va_space", std::ios::in);
  if (!ifs.is_open()) {
    return -1;
  }
  int val{-1};
  if (!(ifs >> val)) {
    return -1;
  }
  return val;
}
//
////------------------------------------------------------------------------alain
//static bool        g_inDocker{};
//static bool        g_inVM{};
//static std::string g_VMenv{};
//static void        DetectSlice()
//{
//  static std::once_flag g_detectSliceOnceFlag{};
//
//  auto detectfunc = []() {
//    std::ifstream ifs("/proc/self/cgroup", std::ios::in);
//    if (!ifs.is_open()) {
//      // cgroup not supported in kernel?
//      g_inDocker = false;
//    } else {
//      for (std::string line; std::getline(ifs, line);) {
//        auto fields = infra::StringUtils::Split(line, ':');
//        if (fields.size() == 3) {
//          if (infra::StringUtils::StartsWith(fields[2], "/docker")) {
//            g_inDocker = true;
//          } else {
//            g_inDocker = false;
//          }
//          break;
//        }
//      }
//    }
//    ifs.close();
//
//    // running in VM shows when CPU flags in CPU info list 'hypervisor'
//    ifs.open("/proc/cpuinfo", std::ios::in);
//    if (!ifs.is_open()) {
//      g_inVM = false;
//    } else {
//      for (std::string line; std::getline(ifs, line);) {
//        if (infra::StringUtils::StartsWith(line, "flags")) {
//          if (strstr(line.c_str(), "hypervisor") != nullptr) {
//            g_inVM = true;
//          } else {
//            g_inVM = false;
//          }
//          break;
//        }
//      }
//    }
//  };
//
//  std::call_once(g_detectSliceOnceFlag, detectfunc);
//}
//
//bool RunningInDocker()
//{
//  DetectSlice();
//  return g_inDocker;
//}
//
//bool RunningInVM()
//{
//  DetectSlice();
//  return g_inVM;
//}
//
//std::string GetVM()
//{
//  static std::once_flag g_detectVMenv{};
//
//  auto detectvm = []() {
//    FILE* fd = popen("/usr/bin/systemd-detect-virt", "r");
//    if (!fd) {
//      g_VMenv = "?";
//      return;
//    }
//
//    char line[64];
//    if (!fgets(line, sizeof(line), fd)) {
//      (void)pclose(fd);
//      g_VMenv = "?";
//      return;
//    }
//    (void)pclose(fd);
//    for (auto i = 0u; i < sizeof(line) and line[i]; ++i) {
//      if (std::isprint(line[i])) {
//        g_VMenv.push_back(line[i]);
//      }
//    }
//  };
//
//  std::call_once(g_detectVMenv, detectvm);
//  return g_VMenv;
//}
//
////------------------------------------------------------------------------alain
//bool VfioDriverPresent()
//{
//  return infra::FileUtils::FileExists("/dev/vfio/vfio");
//}

//------------------------------------------------------------------------alain
bool VfioDriverNoIommuMode()
{
  std::ifstream ifs("/sys/module/vfio/parameters/enable_unsafe_noiommu_mode", std::ios::in);

  if (ifs.is_open()) {
    char yesno;
    ifs >> yesno;
    if (yesno == 'Y' or yesno == 'y') {
      return true;
    }
    if (yesno == 'N' or yesno == 'n') {
      return false;
    }
  }

  //HWAL_FAIL("Invalid value for enable_unsafe_noiommu_mode");
  return false;
}

//------------------------------------------------------------------------alain
std::string VfioDriverVersion()
{
  std::string   verstr;
  std::ifstream ifs("/sys/bus/pci/drivers/vfio-pci/module/version", std::ios::in);
  if (ifs.is_open()) {
    ifs >> verstr;
  }
  return verstr;
}

//------------------------------------------------------------------------alain
uint32_t GetKernelPageSize()
{
  return sysconf(_SC_PAGESIZE);
}

//------------------------------------------------------------------------alain
uint32_t RoundupToKernelPageSize(const uint32_t size)
{
  const auto mask = GetKernelPageSize() - 1;
  return (size + mask) & ~mask;
}

//------------------------------------------------------------------------alain
size_t RoundupToKernelPageSize(const size_t size)
{
  const auto mask = static_cast<size_t>(GetKernelPageSize() - 1);
  return (size + mask) & ~mask;
}
//
////------------------------------------------------------------------------alain
//uint32_t GetHugePageSize()
//{
//  std::ifstream ifs("/proc/meminfo", std::ios::in);
//  if (!ifs) {
//    return 0u;
//  }
//
//  for (std::string line; std::getline(ifs, line);) {
//    if (StringUtils::StartsWithIcmp(line, "hugepagesize:")) {
//      const auto         ptrsize = &line[13];
//      std::istringstream iss(ptrsize);
//
//      // parse count
//      uint32_t cnt;
//      iss >> std::dec >> cnt;
//      if (!iss) {
//        return 0u;
//      }
//      // parse unit
//      std::string unit;
//      iss >> unit;
//      if (!iss) {
//        return 0u;
//      }
//      if (StringUtils::Iequals(unit, "b")) {
//        return cnt * 1;
//      } else if (StringUtils::Iequals(unit, "kb")) {
//        return cnt * 1024;
//      } else if (StringUtils::Iequals(unit, "mb")) {
//        return cnt * 1024 * 1024;
//      } else {
//        // unknown unit
//        INFRA_FAIL("Huge page size not recognized in '", ptrsize, "'");
//      }
//    }
//  }
//  ifs.close();
//
//  return 0u;
//}
//
////------------------------------------------------------------------------alain
//bool GetHugePageStatus(uint32_t& totalcnt, uint32_t& freecnt)
//{
//  return GetHugePageStatusForNode(-1, totalcnt, freecnt);
//}
//
////------------------------------------------------------------------------alain
//bool GetHugePageStatusForNode(const int node, uint32_t& totalcnt, uint32_t& freecnt)
//{
//  constexpr uint32_t HP_TOTAL    = (1u << 0);
//  constexpr uint32_t HP_FREE     = (1u << 1);
//  constexpr uint32_t HP_FOUNDALL = HP_TOTAL | HP_FREE;
//
//  std::string path;
//  std::string prefix;
//  if (node < 0) {
//    // all nodes
//    path   = "/proc/meminfo";
//    prefix = "";
//
//  } else {
//    // node specific
//    path   = std::string("/sys/devices/system/node/node") + std::to_string(node) + "/meminfo";
//    prefix = "node " + std::to_string(node) + " ";
//  }
//
//  std::ifstream ifs(path, std::ios::in);
//  if (!ifs.is_open()) {
//    return false;
//  }
//
//  uint32_t maskfound = 0u;
//  for (std::string line; std::getline(ifs, line);) {
//    if (!prefix.empty()) {
//      if (!StringUtils::StartsWithIcmp(line, prefix)) {
//        // doesn't start with prefix, skip
//        continue;
//      }
//      // remove prefix from line
//      line.erase(0, prefix.size());
//    }
//
//    // check total count
//    if (StringUtils::StartsWithIcmp(line, "hugepages_total:")) {
//      std::istringstream iss(&line[16]);
//      iss >> std::dec >> totalcnt;
//      maskfound |= HP_TOTAL;
//
//      // check free count
//    } else if (StringUtils::StartsWithIcmp(line, "hugepages_free:")) {
//      std::istringstream iss(&line[15]);
//      iss >> std::dec >> freecnt;
//      maskfound |= HP_FREE;
//    }
//
//    if (maskfound == HP_FOUNDALL) {
//      return true;
//    }
//  }
//
//  return false;
//}

//------------------------------------------------------------------------alain
bool GetHugePageRootInfo(std::string& dir, std::string& opts)
{
  std::ifstream ifs("/proc/mounts", std::ios::in);
  if (!ifs.is_open()) {
    return false;
  }
  for (std::string line; std::getline(ifs, line) and !line.empty();) {
    std::istringstream iss(line);

    // parse first field: filesystem type
    std::string dev;
    if (!(iss >> dev)) {
      continue;
    }
    // parse second field: mount point
    std::string mntpt;
    if (!(iss >> mntpt)) {
      continue;
    }

    // save root directory for caller
    dir = mntpt;
    // third field: filesystem type
    std::string fstype;
    if (!(iss >> fstype)) {
      opts.clear();
    }
    if (fstype.compare("hugetlbfs") != 0) {
      // skip non hugetlbfs entries
      continue;
    }

    // supply fourth field as well: mount options
    iss >> opts;
    return true;
  }

  // hugetlbfs not found, hence not mounted
  return false;
}

//------------------------------------------------------------------------alain
static std::string GetOneLine(const std::string& cmd)
{
  FILE* fd = popen(cmd.c_str(), "r");

  char line[1024];
  if (!fgets(line, sizeof(line), fd)) {
    (void)pclose(fd);
    return std::string("");
  }

  (void)pclose(fd);
  return std::string(line, strlen(line) - 1);
}

std::string GetIpAddress()
{
  return GetOneLine("/usr/bin/hostname -i");
}

std::string GetAllIpAddresses()
{
  return GetOneLine("/usr/bin/hostname -I");
}

std::string GetHostName()
{
  return GetOneLine("/usr/bin/hostname -f");
}

std::string GetDomainName()
{
  return GetOneLine("/usr/bin/hostname -d");
}

//------------------------------------------------------------------------alain
int GetSignalFromName(const char* name)
{
  // strsignal() could be used to be more generic, but I'd rather have an
  // error at compile-time than run-time
  // so use a hardcoded table, and just limit to a list of common signals
  // that exists in all Linux variants

  // clang-format off
  struct {
    const char* name;
    const int   num;
  } signals_table[] = {
      {"SIGHUP",   SIGHUP}
    , {"SIGINT",   SIGINT}
    , {"SIGQUIT",  SIGQUIT}
    , {"SIGILL",   SIGILL}
    , {"SIGTRAP",  SIGTRAP}
    , {"SIGABRT",  SIGABRT}
    , {"SIGIOT",   SIGIOT}
    , {"SIGBUS",   SIGBUS}
    , {"SIGFPE",   SIGFPE}
    , {"SIGKILL",  SIGKILL}
    , {"SIGUSR1",  SIGUSR1}
    , {"SIGSEGV",  SIGSEGV}
    , {"SIGUSR2",  SIGUSR2}
    , {"SIGPIPE",  SIGPIPE}
    , {"SIGALRM",  SIGALRM}
    , {"SIGTERM",  SIGTERM}
    , {"SIGCONT",  SIGCONT}
    , {"SIGSTOP",  SIGSTOP}
    , {"SIGTSTP",  SIGTSTP}
    , {"SIGTTIN",  SIGTTIN}
    , {"SIGTTOU",  SIGTTOU}
    , {"SIGURG",   SIGURG}
    , {"SIGIO",    SIGIO}
    , {nullptr,    -1}
  };
  // clang-format on

  for (const auto* sptr = &signals_table[0]; sptr->name; ++sptr) {
    if (strcmp(sptr->name, name) == 0) {
      return sptr->num;
    }
  }
  return -1;
}

int Kill(const int pid, const int signal)
{
  const auto rc = kill(pid, signal);
  return rc;
}

//------------------------------------------------------------------------alain
static bool ServiceCommand(const char* name, const char* action)
{
  const std::string cmd(std::string("service ") + name + " " + action);

  FILE* cmdpipe = popen(cmd.c_str(), "r");
  if (cmdpipe == NULL) {
    printf("  popen failed: %s\n", strerror(errno));
    return false;
  }

  {
    char line[256];
    while (fgets(line, sizeof(line), cmdpipe)) {
      //printf("%s\n", line);
    }
  }

  int rc = pclose(cmdpipe);
  if (rc == -1) {
    printf("  pclose failed: %s\n", strerror(errno));
    return false;
  }

  return true;
}

bool ServiceStart(const char* name)
{
  return ServiceCommand(name, "start");
}

bool ServiceStop(const char* name)
{
  return ServiceCommand(name, "stop");
}

bool ServiceRestart(const char* name)
{
  return ServiceCommand(name, "restart");
}
//
////------------------------------------------------------------------------alain
////------------------------------------------------------------------------alain
//AddressTranslator::AddressTranslator() : m_fd{-1}, m_pagesize{GetKernelPageSize()}
//{
//  static const char* PAGEMAP_FILE = "/proc/self/pagemap";
//  m_fd                            = open(PAGEMAP_FILE, O_RDONLY);
//}
//
//AddressTranslator::~AddressTranslator()
//{
//  if (m_fd >= 0) {
//    close(m_fd);
//  }
//}
//
//uint64_t AddressTranslator::PhysAddrFor(const void* ptr)
//{
//  memdesc_t mdesc;
//
//  if (m_fd < 0 or !m_pagesize) {
//    return BAD_PHYS_ADDR;
//  }
//
//  // compute page frame number
//  const uint64_t vaddr = reinterpret_cast<uint64_t>(ptr);
//  const uint64_t pfn   = vaddr / m_pagesize;
//  // read pfn'th memory descriptor, each being 64-bit
//  const off_t offset = static_cast<off_t>(pfn) * sizeof(memdesc_t);
//  if (lseek(m_fd, offset, SEEK_SET) == static_cast<off_t>(-1)) {
//    return BAD_PHYS_ADDR;
//  }
//  const int rc = read(m_fd, &mdesc, sizeof(mdesc));
//  if (rc < 0) {
//    return BAD_PHYS_ADDR;
//  } else if (rc != sizeof(memdesc_t)) {
//    return BAD_PHYS_ADDR;
//  }
//  if ((mdesc & 0x7fffffffffffffull) == 0) {
//    return BAD_PHYS_ADDR;
//  }
//
//  // compute physical address
//  return ((mdesc & 0x7fffffffffffffull) * m_pagesize) + (vaddr % m_pagesize);
//}

}  // namespace infra
