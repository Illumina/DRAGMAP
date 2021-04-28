// Copyright 2013 Edico Genome Corporation. All rights reserved.
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

#ifndef HOST_VERSION_HASHTABLE_H
#define HOST_VERSION_HASHTABLE_H

#include <stdint.h>

//--------------------------------------------------------------------------------------------------
// Return the host software version:
//
//   <HW Major Version>.<HW Minor Version>.<HW Build Version>.
//   <SW Major Version>.<SW Minor Version>.<Branch ID>.<Git Commit ID>
//
// The perforce change list is the highest changelist associated with the file revisions
// in the current workspace and is not included in customer releases.
//

#ifdef __cplusplus
extern "C" {
#endif

//--------------------------------------------------------------------------------------------------
//
// Hardware Revision Description:
//
// 32-bit revision code will be MM.mmm.bbb/BBBB.mcs, where:
//
// - MM  = Major revision.  This will change at Edico's discretion for major releases, feature addition, etc.
// - mmm = Minor revision.  This will change when the new HW will not be backwards compatible with the
// previous.
// - bbb = Build revision.  This will increment for each new bitstream that is generated.
// - BBBB = Build type.  This indicates what is present in the bitstream.  0000 = map/align, 0001 = diag, 0003
// = hmm, 0004= combo
//
// If the major/minor revision is updated, then all sub-revisions will reset to 0.
//
//--------------------------------------------------------------------------------------------------

#define HW_MAJOR_VERSION (1)
#define HW_MINOR_VERSION (2)
#define HW_BUILD_VERSION (54)
// this is for the rpm build, please update version if HW major, minor or build version changes
// #define HW_VERSION           01.002.054

#define HW_BASE_MAJOR_VERSION (1)
#define HW_BASE_MINOR_VERSION (2)
#define HW_BASE_BUILD_VERSION (50)
// this is for the rpm build, please update version if HW base major, minor or build version changes
// #define HW_BASE_VERSION             01.002.050

#define HW_DNA_MAPPER_MAJOR_VERSION (1)
#define HW_DNA_MAPPER_MINOR_VERSION (2)
#define HW_DNA_MAPPER_BUILD_VERSION (54)
// this is for the rpm build, please update version if HW DNA mapper major, minor or build version changes
// #define HW_DNA_MAPPER_VERSION       01.002.054

#define HW_RNA_MAPPER_MAJOR_VERSION (1)
#define HW_RNA_MAPPER_MINOR_VERSION (2)
#define HW_RNA_MAPPER_BUILD_VERSION (54)
// this is for the rpm build, please update version if HW RNA mapper major, minor or build version changes
// #define HW_RNA_MAPPER_VERSION       01.002.054

#define HW_HMM_MAJOR_VERSION (1)
#define HW_HMM_MINOR_VERSION (2)
#define HW_HMM_BUILD_VERSION (50)
// this is for the rpm build, please update version if HW HMM major, minor or build version changes
// #define HW_HMM_VERSION              01.002.050

//--------------------------------------------------------------------------------------------------
// Return the minimum hardware version supported with the host software.
//
//

// Return the HW version
const char* getHardwareVersion();

// Return the HW base version
const char* getHardwareBaseVersion();

// Return the DNA mapper version
const char* getHardwareDnaMapperVersion();

// Return the RNA mapper version
const char* getHardwareRnaMapperVersion();

// Return the HMM version
const char* getHardwareHmmVersion();

// Return the host software version
const char* getHostVersion();

#ifdef __cplusplus
}  // extern "C"
#endif

#endif
