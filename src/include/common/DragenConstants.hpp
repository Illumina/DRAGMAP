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

#ifndef COMMON_DRAGEN_CONSTANTS_HPP
#define COMMON_DRAGEN_CONSTANTS_HPP

#include <unistd.h>  // for _SC_PAGE_SIZE

#include <cstdint>

namespace DragenConstants {

//
// General
//
const uint64_t DRAGEN_1KB   = 0x400;
const uint64_t DRAGEN_4KB   = 0x1000;
const uint64_t DRAGEN_8KB   = 0x2000;
const uint64_t DRAGEN_16KB  = 0x4000;
const uint64_t DRAGEN_32KB  = 0x8000;
const uint64_t DRAGEN_64KB  = 0x10000;
const uint64_t DRAGEN_128KB = 0x20000;
const uint64_t DRAGEN_256KB = 0x40000;
const uint64_t DRAGEN_512KB = 0x80000;
const uint64_t DRAGEN_1MB   = 1 << 20;
const uint64_t DRAGEN_1GB   = 1 << 30;

// We have a hard-coded limit to the number of read groups we can handle,
// due to a corresponding hardware limit
const uint16_t DRAGEN_MAX_READ_GROUPS = 1024;

// The line for whether a dataset is an exome or a genome
const double DRAGEN_EXOME_CUTOFF                = 50.0;
const double DRAGEN_GERMLINE_EXOME_CUTOFF       = 50.0;
const double DRAGEN_SOMATIC_EXOME_CUTOFF        = 40.0;
const double DRAGEN_SOMATIC_NORMAL_EXOME_CUTOFF = 80.0;

typedef enum vcType { NONE, GERMLINE, SOMATIC, SOMATIC_NORMAL } vcType_t;
}  // namespace DragenConstants

#define DRAGEN_PAGE_SIZE sysconf(_SC_PAGE_SIZE)

#endif  // #ifndef COMMON_DRAGEN_CONSTANTS_HPP
