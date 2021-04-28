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

#pragma once

#include <cstdint>
#include <string>

namespace infra {

int         GetDmiValue(const std::string& label, std::string& value);
std::string getExecutablePath();
std::string getExecutableName();
std::string getExecutableDirectory();

uint32_t    GetNumSystemThreads();
std::string GetUserId();
std::string GetRealUserId();
uint64_t    GetSystemMemorySize();
uint32_t    GetCacheLineSize();

std::string GetKernelName();
std::string GetKernelArch();
std::string GetKernelVersionFull();
std::string GetKernelVersion();
int         GetKernelRandomizeAddressStatus();
bool        KernelVersionAtLeast(const int major, const int minor = 0, const int patch = 0);
bool        RunningInDocker();
bool        RunningInVM();
std::string GetVM();
bool        VfioDriverPresent();
bool        VfioDriverNoIommuMode();
std::string VfioDriverVersion();

uint32_t GetKernelPageSize();
uint32_t RoundupToKernelPageSize(const uint32_t size);
size_t   RoundupToKernelPageSize(const size_t size);
uint32_t GetHugePageSize();
bool     GetHugePageStatus(uint32_t& totalcnt, uint32_t& freecnt);
bool     GetHugePageStatusForNode(const int node, uint32_t& totalcnt, uint32_t& freecnt);
bool     GetHugePageRootInfo(std::string& dir, std::string& opts);

std::string GetAllIpAddresses();
std::string GetIpAddress();
std::string GetHostName();
std::string GetDomainName();

int GetSignalFromName(const char* name);
int Kill(const int pid, const int signal);

bool ServiceStart(const char* name);
bool ServiceStop(const char* name);
bool ServiceRestart(const char* name);

//// convert from virtual address to physical address
//class AddressTranslator {
//  using memdesc_t = uint64_t;
//
//public:
//  static constexpr uint64_t BAD_PHYS_ADDR = 0ull;
//
//  AddressTranslator();
//  ~AddressTranslator();
//  uint64_t PhysAddrFor(const void* vaddr);
//
//private:
//  int      m_fd;
//  uint32_t m_pagesize;
//};

}  // namespace infra
