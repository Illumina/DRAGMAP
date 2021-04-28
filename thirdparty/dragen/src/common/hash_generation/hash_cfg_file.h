// Copyright 2013-2015 Edico Genome Corporation. All rights reserved.
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
//

#ifndef _hash_cfg_file_h_
#define _hash_cfg_file_h_

#include "gen_hash_table.h"

#ifdef __cplusplus
extern "C" {
#endif

// Return any error message reported during hash table processing
const char* getHashCfgError();

// Return message with instructions for rebuilding the hash table
const char* getBuildHashTableCmdline(const char* binFile, const char* liftoverFile);

// Write binary hash table configuration file
void writeHashCfgBin(const char* filename, const hashTableConfig_t* config);
// Write text hash table configuration file
void writeHashCfgTxt(const char* filename, const hashTableConfig_t* config);

// Read binary hash table configuration file
int readHashCfgBin(const char* filename, hashTableConfig_t* config);
// Read text hash table configuration file
int readHashCfgTxt(const char* filename, hashTableConfig_t* config);

// Helper function to check for really old hash tables without a version string
int isOldHashTable(const char* filename);

#ifdef __cplusplus
}
#endif

#endif
