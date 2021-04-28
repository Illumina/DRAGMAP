// Copyright 2017 Edico Genome Corporation. All rights reserved.
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
#ifndef _gen_hash_table_compress_h_
#define _gen_hash_table_compress_h_

#ifdef __cplusplus
extern "C" {
#endif

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>

#include "hash_table.h"

#define COMP_MAGIC 0xF89B13CD
#define COMP_VERSION 0x00000002
#define COMP_NOLITERAL_RUN_BITS 6
#define COMP_NOLITERAL_RUN_LIMIT ((1 << COMP_NOLITERAL_RUN_BITS) - 1)
#define COMP_EXT_NOLITERAL_RUN_BITS 10
#define COMP_EXT_NOLITERAL_RUN_LIMIT ((1 << COMP_EXT_NOLITERAL_RUN_BITS) - 1)

#define EXTTAB_INDEX_BUCKET_BITS 8

// Information to encode an "automatic" hit record
#define SEED_POP_REC_SLOT_OFFSET_BITS 21
#define SEED_POP_REC_BUCKET_BIN_BITS 22  // Supports EXTTAB_INDEX_BUCKET_BITS >= 7
#define SEED_POP_REC_LIFT_GROUP_BITS 17
typedef struct {
  uint64_t slotOffset : SEED_POP_REC_SLOT_OFFSET_BITS;
  uint64_t bucketBin : SEED_POP_REC_BUCKET_BIN_BITS;
  uint64_t liftGroup : SEED_POP_REC_LIFT_GROUP_BITS;
  uint64_t liftCode : 2;
  uint64_t extTable : 1;
  uint64_t automatic : 1;
} seedPopRec_t;

typedef struct {
  uint32_t start;
  uint32_t length;
  uint8_t  offsetBits;
} extIndexRec_t;

// Context structure for compressing a hash table
typedef struct {
  uint8_t* buf;
  int      bufBytes;
  FILE*    file;
  uint64_t bitLen;
  int64_t  bufStartBit;
  uint8_t  bitsPend;
  uint64_t litPos;
  uint64_t extLitPos;
  uint32_t extTabRecs;
  uint64_t extTabRecsFilePos;
  uint32_t extendIdBits;
  uint64_t autoNulRecs, autoPriRecs, autoSecRecs, specialRecs, literalRecs, chainPtrRecs, chainEndRecs,
      extLitRecs, totalRecs;
  uint64_t autoNulBits, autoPriBits, autoSecBits, specialBits, literalBits, chainPtrBits, chainEndBits,
      extLitBits, totalBits, miscBits;
} writeCompHashTableCtx_t;

// Compression functions, buffers to file - integrate into hash table builder
int writeCompHashTableCtxOpen(writeCompHashTableCtx_t* ctx, char* fileName);
int writeCompHashTableAlignAndFlush(writeCompHashTableCtx_t* ctx);
int writeCompHashTableClose(writeCompHashTableCtx_t* ctx);
int writeCompHashTableHeader(writeCompHashTableCtx_t* ctx, hashTableHeader_t* cfgHdr);
int writeCompHashTableLiterals(
    writeCompHashTableCtx_t* ctx, hashrec_t* recs, uint64_t numRecs, uint8_t* litFlags);
int writeCompHashTableExtTabLiterals(writeCompHashTableCtx_t* ctx, extend_hit_t* recs, uint64_t numRecs);
int writeCompHashTableLiteralsClose(writeCompHashTableCtx_t* ctx);
int writeCompHashTableExtIndex(
    writeCompHashTableCtx_t* ctx, extIndexRec_t* extIndexRecs, uint32_t numExtIndexRecs, uint32_t extTabRecs);
int writeCompHashTableAutoHits(
    writeCompHashTableCtx_t* ctx,
    seedPopRec_t*            popRecs,
    uint32_t                 numPopRecs,
    uint32_t                 maxExtendId,
    extIndexRec_t*           extIndexRecs);

// Decompress function, buffer to buffer.
// Pre-allocated buffer may be provided, else one will get malloc'ed.
// Pass buffers containing compressed hash table (hash_table.cmp) and binary reference (reference.bin).
// Uncompressed hash table digest is calculated unless the digest pointer is NULL.
// Returns NULL on success, or error message.
char* decompHashTable(
    int       threads,
    uint8_t*  compBuf,
    uint64_t  compLen,
    uint8_t*  refBuf,
    uint64_t  refLen,
    uint8_t** decompBuf,
    uint64_t* decompLen,
    uint8_t** extendTableBuf,
    uint64_t* extendTableLen,
    uint32_t* hashDigest,
    uint32_t* extTabDigest);

// Context structure for decompressing a hash table
typedef struct {
  int               threads;
  uint8_t*          buf;
  int64_t           bufBytes;
  int64_t           bufBits;
  int64_t           bitPos;
  int64_t           oldBitPos;
  hashTableHeader_t cfgHdr;
  hashrec_t*        hashTable;
  extend_hit_t*     extendTable;
  extIndexRec_t*    extIndexRecs;
  int64_t           hashRecs;
  int64_t           extendRecs;
  uint8_t*          refBuf2Bit;
  int64_t           refBases;
  uint32_t          extendIdBits;
  void*             priCrcInit;
  void*             secCrcInit;
} decompHashTableCtx_t;

// For testing
void slurpFile(const char* name, uint8_t** buffer, uint64_t* bufLen);
void dumpFile(const char* name, uint8_t* buffer, uint64_t bufLen);
int  decompAndWriteHashTable(
     const char* refPath,
     const char* hashCmpPath,
     const char* hashBinPath,
     const char* extTabPath,
     const int   numThreads);

#ifdef __cplusplus
}
#endif

#endif
