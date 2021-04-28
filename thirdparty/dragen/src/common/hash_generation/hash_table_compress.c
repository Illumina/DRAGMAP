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

#include "hash_table_compress.h"
#include <inttypes.h>
#include <math.h>
#include <pthread.h>
#include <semaphore.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include "crc_hash.h"
#include "gen_hash_table.h"
#include "hash_table.h"
#ifndef LOCAL_BUILD
#include "crc32_hw.h"
#include "watchdog.h"
#ifdef _TARGET_PPC_
#include "crc32_powerpc.h"
#endif
/* ============ end of #include directives ============ */
#endif

// The following filter out these specific compiler warnings
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wformat"
#pragma GCC diagnostic ignored "-Wsequence-point"
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"

#ifdef LOCAL_BUILD
#define DEBUG_TIME 1
#else
#define DEBUG_TIME 0
#endif
#define DEBUG_AUTO 0  // (seedIndex == 0x2D48B61C)
#define DEBUG_FLAG 0  // (initBucket == 172438)
#define DEBUG_THREAD 0

#define PRINTIF(c, x) \
  if (c) {            \
    printf x;         \
    fflush(stdout);   \
  }

#define MAGIC_FILE_START (COMP_MAGIC + 0x00)
#define MAGIC_AFTER_HEADER (COMP_MAGIC + 0x01)
#define MAGIC_LIT_BLOCK (COMP_MAGIC + 0x12)
#define MAGIC_LIT_END (COMP_MAGIC + 0x13)
#define MAGIC_EXT_IDX_START (COMP_MAGIC + 0x31)
#define MAGIC_EXT_IDX_END (COMP_MAGIC + 0x032)
#define MAGIC_AUTO_START (COMP_MAGIC + 0x21)
#define MAGIC_AUTO_BLOCK (COMP_MAGIC + 0x22)
#define MAGIC_AUTO_END (COMP_MAGIC + 0x23)

#define MTHREAD_MAX_THREADS 256
// Shared context for multithreading
typedef struct {
  void*           ctx;
  void**          inpWorkPtrs;
  void**          outWorkPtrs;
  void*           workUnitsBuf;
  char*           errMsg;
  int             numThreads;
  int             nextThreadId;
  int             workBytes;
  int             workUnits;
  int             inpWorkHead;
  int             outWorkHead;
  int             inpWorkTail;
  int             outWorkTail;
  int             inpWorkSemVld;
  int             outWorkSemVld;
  int             inpWorkMutexVld;
  int             outWorkMutexVld;
  int             errMsgMutexVld;
  sem_t           inpWorkSem;
  sem_t           outWorkSem;
  pthread_mutex_t inpWorkMutex;
  pthread_mutex_t outWorkMutex;
  pthread_mutex_t errMsgMutex;
  pthread_t*      threads;
  int8_t          threadsActive[MTHREAD_MAX_THREADS];
} mthreadShare_t;

// Work unit for threads populating blocks of literal records or automatic hit records
typedef struct {
  int64_t startPos;
  int64_t endPos;
  int64_t startBit;
  int64_t endBit;
} blockWorkRec_t;

void mthreadCancel(mthreadShare_t* share)
{
  int i;
  if (!share) return;
  if (share->inpWorkMutexVld) {
    PRINTIF(DEBUG_THREAD, ("Terminating remaining child threads...\n"));
    pthread_mutex_lock(&share->inpWorkMutex);
    memset(share->inpWorkPtrs, 0, share->workUnits * sizeof(void*));
    pthread_mutex_unlock(&share->inpWorkMutex);
    for (i = 0; i < share->numThreads; i++) sem_post(&share->inpWorkSem);
    for (i = 0; i < share->numThreads; i++) {
      if (share->threadsActive[i]) {
        PRINTIF(DEBUG_THREAD, ("  Joining child thread %d\n", i));
        if (pthread_join(share->threads[i], NULL)) {
          PRINTIF(DEBUG_THREAD, ("    Failed, canceling active thread %d\n", i));
          pthread_cancel(share->threads[i]);
        }
      }
    }
    PRINTIF(DEBUG_THREAD, ("  All threads done\n"));
  }
  free(share->workUnitsBuf);
  free(share->inpWorkPtrs);
  free(share->outWorkPtrs);
  free(share->threads);
  if (share->inpWorkSemVld) sem_destroy(&share->inpWorkSem);
  if (share->outWorkSemVld) sem_destroy(&share->outWorkSem);
  if (share->inpWorkMutexVld) pthread_mutex_destroy(&share->inpWorkMutex);
  if (share->outWorkMutexVld) pthread_mutex_destroy(&share->outWorkMutex);
  if (share->errMsgMutexVld) pthread_mutex_destroy(&share->errMsgMutex);
  free(share);
}

int mthreadGetThreadId(mthreadShare_t* share)
{
  pthread_mutex_lock(&share->inpWorkMutex);
  int id = share->nextThreadId++;
  pthread_mutex_unlock(&share->inpWorkMutex);
  return id;
}

void mthreadSetErrMsg(mthreadShare_t* share, char* errMsg)
{
  pthread_mutex_lock(&share->errMsgMutex);
  if (!share->errMsg) share->errMsg = errMsg;
  pthread_mutex_unlock(&share->errMsgMutex);
}

mthreadShare_t* mthreadInit(void* ctx, void* (*threadFunction)(void*), int numThreads, int workBytes)
{
  int i;
  // Constrain parameters
  if (numThreads < 1) numThreads = 1;
  if (numThreads > MTHREAD_MAX_THREADS) numThreads = MTHREAD_MAX_THREADS;
  int workUnits = numThreads * 8;
  PRINTIF(
      DEBUG_THREAD,
      ("Initializing multithreading: ctx=%p, numThreads=%d, workBytes=%d, workUnits=%d\n",
       ctx,
       numThreads,
       workBytes,
       workUnits));
  // Allocate a shared context for parent and child threads
  mthreadShare_t* share = calloc(1, sizeof(mthreadShare_t));
  if (!share) goto mthreadInitError;
  share->ctx        = ctx;
  share->workUnits  = workUnits;
  share->workBytes  = workBytes;
  share->numThreads = numThreads;
  // Allocate arrays inside
  share->workUnitsBuf = calloc(workUnits, workBytes);
  share->inpWorkPtrs  = calloc(workUnits, sizeof(void*));
  share->outWorkPtrs  = calloc(workUnits, sizeof(void*));
  share->threads      = calloc(numThreads, sizeof(pthread_t));
  if (!share->workUnitsBuf || !share->inpWorkPtrs || !share->outWorkPtrs || !share->threads)
    goto mthreadInitError;
  // Initialize semaphores and mutexes
  if (sem_init(&share->inpWorkSem, 0, 0)) goto mthreadInitError;
  share->inpWorkSemVld = 1;
  if (sem_init(&share->outWorkSem, 0, 0)) goto mthreadInitError;
  share->outWorkSemVld = 1;
  if (pthread_mutex_init(&share->inpWorkMutex, 0)) goto mthreadInitError;
  share->inpWorkMutexVld = 1;
  if (pthread_mutex_init(&share->outWorkMutex, 0)) goto mthreadInitError;
  share->outWorkMutexVld = 1;
  if (pthread_mutex_init(&share->errMsgMutex, 0)) goto mthreadInitError;
  share->errMsgMutexVld = 1;
  // Put work batches on the free list
  pthread_mutex_lock(&share->outWorkMutex);
  for (i = 0; i < share->workUnits; i++) {
    PRINTIF(DEBUG_THREAD, ("  Readying work unit %d\n", i));
    share->outWorkPtrs[share->outWorkHead] = &((uint8_t*)share->workUnitsBuf)[i * workBytes];
    share->outWorkHead                     = (share->outWorkHead + 1) % share->workUnits;
    sem_post(&share->outWorkSem);
  }
  pthread_mutex_unlock(&share->outWorkMutex);
  // Start worker threads
  for (i = 0; i < numThreads; i++) {
    PRINTIF(DEBUG_THREAD, ("  Spawning worker thread %d\n", i));
    if (pthread_create(&share->threads[i], NULL, threadFunction, share)) goto mthreadInitError;
    share->threadsActive[i] = 1;
  }
  // Success
  return share;
mthreadInitError:
  mthreadCancel(share);
  return NULL;
}

#define WORK_UNIT_ID (((uint8_t*)work - (uint8_t*)share->workUnitsBuf) / share->workBytes)
void* mthreadParentGetJob(mthreadShare_t* share)
{
  sem_wait(&share->outWorkSem);
  pthread_mutex_lock(&share->outWorkMutex);
  void* work         = share->outWorkPtrs[share->outWorkTail];
  share->outWorkTail = (share->outWorkTail + 1) % share->workUnits;
  pthread_mutex_unlock(&share->outWorkMutex);
  if (!work) mthreadSetErrMsg(share, "Null work pointer returned by decompression thread");
  if (!work || share->errMsg) return NULL;
  PRINTIF(DEBUG_THREAD, ("Opening empty work unit %ld\n", WORK_UNIT_ID));
  return work;
}

void mthreadParentPutJob(mthreadShare_t* share, void* work)
{
  PRINTIF(DEBUG_THREAD && work, ("Posting work unit %ld\n", WORK_UNIT_ID));
  pthread_mutex_lock(&share->inpWorkMutex);
  share->inpWorkPtrs[share->inpWorkHead] = work;
  share->inpWorkHead                     = (share->inpWorkHead + 1) % share->workUnits;
  pthread_mutex_unlock(&share->inpWorkMutex);
  sem_post(&share->inpWorkSem);
}

void* mthreadChildGetJob(int threadId, mthreadShare_t* share)
{
  PRINTIF(DEBUG_THREAD, ("Thread %d waiting for work\n", threadId));
  sem_wait(&share->inpWorkSem);
  pthread_mutex_lock(&share->inpWorkMutex);
  void* work         = share->inpWorkPtrs[share->inpWorkTail];
  share->inpWorkTail = (share->inpWorkTail + 1) % share->workUnits;
  pthread_mutex_unlock(&share->inpWorkMutex);
  if (!work) {
    PRINTIF(DEBUG_THREAD, ("Thread %d got NULL work, should terminate\n", threadId));
    return NULL;
  }
  PRINTIF(DEBUG_THREAD, ("Thread %d starting on work unit %ld\n", threadId, WORK_UNIT_ID));
  return work;
}

void mthreadChildPutJob(int threadId, mthreadShare_t* share, void* work)
{
  PRINTIF(DEBUG_THREAD, ("Thread %d finished work unit %ld\n", threadId, WORK_UNIT_ID));
  PRINTIF(DEBUG_THREAD && share->errMsg, ("  Noticed error message: %s\n", share->errMsg));
  pthread_mutex_lock(&share->outWorkMutex);
  share->outWorkPtrs[share->outWorkHead] = work;
  share->outWorkHead                     = (share->outWorkHead + 1) % share->workUnits;
  pthread_mutex_unlock(&share->outWorkMutex);
  sem_post(&share->outWorkSem);
}

char* mthreadClose(mthreadShare_t* share)
{
  int i;
  PRINTIF(DEBUG_THREAD, ("Closing multi-threading\n"));
  if (!share) return NULL;
  // Post one NULL job per thread to terminate them
  for (i = 0; i < share->numThreads; i++) {
    PRINTIF(DEBUG_THREAD, ("Posting null work unit %d\n", i));
    // Fetch an empty work batch (not used)
    mthreadParentGetJob(share);
    mthreadParentPutJob(share, NULL);
  }
  // Wait for threads to return
  for (i = 0; i < share->numThreads; i++) {
    PRINTIF(DEBUG_THREAD, ("Joining child thread %d\n", i));
    if (!pthread_join(share->threads[i], NULL)) share->threadsActive[i] = 0;
  }
  PRINTIF(DEBUG_THREAD, ("All threads joined\n"));
  char* errMsg = share->errMsg;
  mthreadCancel(share);
  return errMsg;
}

#define BLOCK_LIM 1000000
#define WRITE_BUF_LEN (BLOCK_LIM + 1000)
#define WRITE_BUF_LIM (BLOCK_LIM + 900)

// Works only for len = 1-64
int lclPutBits(writeCompHashTableCtx_t* ctx, uint64_t val, int len)
{
  int cur = (ctx->bitLen & 7);
  val &= 0xFFFFFFFFFFFFFFFFULL >> (64 - len);
  uint8_t ofl                            = val >> (64 - cur);
  val                                    = (val << cur) | ctx->bitsPend;
  int sum                                = cur + len;
  int wbits                              = sum & 0xF8;
  *(uint64_t*)(&ctx->buf[ctx->bufBytes]) = val;
  ctx->bufBytes += wbits >> 3;
  ctx->bitsPend = sum > 64 ? ofl : wbits == 64 ? 0 : val >> wbits;
  ctx->bitLen += len;
  if (ctx->bufBytes < WRITE_BUF_LIM) return 0;
  if (fwrite(ctx->buf, ctx->bufBytes, 1, ctx->file) != 1) return 1;
  ctx->bufBytes    = 0;
  ctx->bufStartBit = ctx->bitLen;
  return 0;
}

#if 0
#define PUT_QWORD(val)                      \
  {                                         \
    if (lclPutBits(ctx, val, 64)) return 1; \
  }
#define PUT_BITS(val, len)                   \
  {                                          \
    if (lclPutBits(ctx, val, len)) return 1; \
  }
#define PUT_FLUSH                                  \
  {                                                \
    if (writeCompHashTableAlign(ctx, 1)) return 1; \
  }
#define PUT_ALIGN                                  \
  {                                                \
    if (writeCompHashTableAlign(ctx, 0)) return 1; \
  }
#else

#define PUT_QWORD(val)                                                \
  {                                                                   \
    if (cacheBits && lclPutBits(ctx, cacheWord, cacheBits)) return 1; \
    if (lclPutBits(ctx, val, 64)) return 1;                           \
    cacheWord = 0;                                                    \
    cacheBits = 0;                                                    \
  }
#define PUT_BITS(val, len)                                                 \
  {                                                                        \
    if ((len) + cacheBits > 64) {                                          \
      if (cacheBits && lclPutBits(ctx, cacheWord, cacheBits)) return 1;    \
      cacheWord = (uint64_t)(val) & ((1ULL << (len)) - 1);                 \
      cacheBits = (len);                                                   \
    } else {                                                               \
      cacheWord |= ((uint64_t)(val) & ((1ULL << (len)) - 1)) << cacheBits; \
      cacheBits += (len);                                                  \
    }                                                                      \
  }
#define PUT_FLUSH                                                     \
  {                                                                   \
    if (cacheBits && lclPutBits(ctx, cacheWord, cacheBits)) return 1; \
    if (writeCompHashTableAlign(ctx, 1)) return 1;                    \
    cacheWord = 0;                                                    \
    cacheBits = 0;                                                    \
  }
#define PUT_ALIGN                                                     \
  {                                                                   \
    if (cacheBits && lclPutBits(ctx, cacheWord, cacheBits)) return 1; \
    if (writeCompHashTableAlign(ctx, 0)) return 1;                    \
    cacheWord = 0;                                                    \
    cacheBits = 0;                                                    \
  }
#endif

int writeCompHashTableCtxOpen(writeCompHashTableCtx_t* ctx, char* fileName)
{
  memset(ctx, 0, sizeof(writeCompHashTableCtx_t));
  ctx->buf = malloc(WRITE_BUF_LEN);
  if (!ctx->buf) return 1;
  if (!fileName) return 1;
  ctx->file = fopen(fileName, "wb");
  if (!ctx->file) {
    free(ctx->buf);
    return 1;
  }
  return 0;
}

int writeCompHashTableAlign(writeCompHashTableCtx_t* ctx, int flush)
{
  int padLen = (8 - (ctx->bitLen & 7)) % 8;
  if (padLen) lclPutBits(ctx, 0, padLen);
  ctx->miscBits += padLen;
  if (!flush || !ctx->bufBytes) return 0;
  if (fwrite(ctx->buf, ctx->bufBytes, 1, ctx->file) != 1) return 1;
  ctx->bufBytes    = 0;
  ctx->bufStartBit = ctx->bitLen;
  return 0;
}

int writeCompHashTableClose(writeCompHashTableCtx_t* ctx)
{
  if (!ctx->file) return 0;
  if (writeCompHashTableAlign(ctx, 1)) return 1;
  free(ctx->buf);
  ctx->buf = NULL;
  // Seek back to write extend table length (not known up front) to its header field
  if (fseeko(ctx->file, ctx->extTabRecsFilePos, SEEK_SET)) return 1;
  if (fwrite(&ctx->extTabRecs, 4, 1, ctx->file) != 1) return 1;
  if (fclose(ctx->file)) return 1;
  ctx->file = NULL;
  return 0;
}

// Writes comressed hash table file header, returning nonzero on error
int writeCompHashTableHeader(writeCompHashTableCtx_t* ctx, hashTableHeader_t* cfgHdr)
{
  uint64_t cacheWord = 0;
  size_t   cacheBits = 0;
  uint8_t* hdr       = (uint8_t*)cfgHdr;
  size_t   i;
  // Magic number for this file type
  PUT_BITS(MAGIC_FILE_START, 32);
  ctx->extTabRecsFilePos += 4;
  // Compressed format version
  PUT_BITS(COMP_VERSION, 32);
  ctx->extTabRecsFilePos += 4;
  // Dummy extTabRecs value, since we don't know the real one yet
  cfgHdr->extTabRecs = 0xFFFFFFFF;
  // Entire configuration header
  for (i = 0; i < sizeof(hashTableHeader_t); i++) PUT_BITS(hdr[i], 8);
  // Remember where the header's extTabRecs field was in the file (in bytes)
  ctx->extTabRecsFilePos += ((uint8_t*)(&cfgHdr->extTabRecs) - (uint8_t*)cfgHdr);
  // Magic number to verify synchronization
  PUT_BITS(MAGIC_AFTER_HEADER, 32);
  // Flush
  PUT_FLUSH;
  ctx->miscBits = ctx->bitLen;
  return 0;
}

// Routine to encode the records from a hash table chunk which have corresponding flags for literal output
int writeCompHashTableLiterals(
    writeCompHashTableCtx_t* ctx, hashrec_t* recs, uint64_t numRecs, uint8_t* litFlags)
{
  uint64_t  cacheWord = 0;
  int       cacheBits = 0;
  uint64_t  i, oldBits;
  uint32_t  opc, con, run = 0;
  hashrec_t rec;
  uint64_t* blockLenPtr   = NULL;
  uint64_t* blockEndPtr   = NULL;
  int64_t   blockStartBit = 0;
  // Flush buffer so all room is available
  PUT_FLUSH;
  for (i = 0; 1; i++) {
    // Insert block boundary when the buffer nears full
    if (!i || ctx->bufBytes + (((ctx->bitLen & 7) + cacheBits) >> 3) >= BLOCK_LIM || i == numRecs) {
      // At end, encode remaining run length
      if (i == numRecs) {
        oldBits = ctx->bitLen + cacheBits;
        if (run) {
          for (; run >= COMP_NOLITERAL_RUN_LIMIT; run -= COMP_NOLITERAL_RUN_LIMIT)
            PUT_BITS(COMP_NOLITERAL_RUN_LIMIT, COMP_NOLITERAL_RUN_BITS);
          PUT_BITS(run, COMP_NOLITERAL_RUN_BITS);
        }
        ctx->miscBits += ctx->bitLen + cacheBits - oldBits;
      }
      // Align to byte boundary
      PUT_ALIGN;
      // Fill in previous block's end record and coded bit length
      if (i) {
        *blockEndPtr = i + ctx->litPos;
        *blockLenPtr = ctx->bitLen - blockStartBit;
      }
      // Flush buffer so all room is available
      PUT_FLUSH;
      // Done?
      if (i == numRecs) break;
      // Magic number for synchronization
      PUT_BITS(MAGIC_LIT_BLOCK, 32);
      ctx->miscBits += 32;
      // Block start record
      PUT_QWORD(i + ctx->litPos);
      ctx->miscBits += 64;
      // Placeholder for block end record
      blockEndPtr = (uint64_t*)&ctx->buf[(ctx->bitLen - ctx->bufStartBit) >> 3];
      PUT_QWORD(0);
      ctx->miscBits += 64;
      // Placeholder for block encoded length
      blockLenPtr = (uint64_t*)&ctx->buf[(ctx->bitLen - ctx->bufStartBit) >> 3];
      PUT_QWORD(0);
      ctx->miscBits += 64;
      blockStartBit = ctx->bitLen;
      // Flag indicating block type is hash table literals
      PUT_BITS(0, 1);
    }
    // Just tally records without literal flags
    if (!((litFlags[i >> 3] >> (i & 7)) & 1)) {
      run++;
      continue;
    }
    // Encode run length before the next literal record.  Repeat the maximum length until the remainder is
    // less.
    oldBits = ctx->bitLen + cacheBits;
    for (; run >= COMP_NOLITERAL_RUN_LIMIT; run -= COMP_NOLITERAL_RUN_LIMIT) {
      PUT_BITS(COMP_NOLITERAL_RUN_LIMIT, COMP_NOLITERAL_RUN_BITS);
    }
    PUT_BITS(run, COMP_NOLITERAL_RUN_BITS);
    run = 0;
    rec = recs[i];
    opc = rec.general.opcode;
    // Non-Chain record: write the whole thing
    if (!HASH_OPC_IS_CHAIN(opc)) {
      PUT_BITS(0, 1);        // Chain flag
      PUT_QWORD(rec.qword);  // Literal record
      if (HASH_OPC_IS_HIT(opc)) {
        ctx->specialBits += ctx->bitLen + cacheBits - oldBits;
        ctx->specialRecs++;
      } else {
        ctx->literalBits += ctx->bitLen + cacheBits - oldBits;
        ctx->literalRecs++;
      }
      continue;
    }
    // Chain terminator: short code
    if (rec.qword == HASH_REC_CHAIN_TERM_QWORD) {
      PUT_BITS(3, 2);  // Chain flag, term flag
      ctx->chainEndBits += ctx->bitLen + cacheBits - oldBits;
      ctx->chainEndRecs++;
      continue;
    }
    // Chain pointer: partial content
    con = HASH_OPC_IS_CON(opc);
    PUT_BITS(1 | (con << 2), 3);        // Chain flag, term flag, continuation flag
    PUT_BITS(rec.chain.chain_ptr, 18);  // Chain pointer
    ctx->chainPtrBits += ctx->bitLen + cacheBits - oldBits;
    ctx->chainPtrRecs++;
  }

  // Update the position in the literal records stream
  ctx->litPos += numRecs;
  return 0;
}

// Routine to encode the records from an extension table chunk which have corresponding flags for literal
// output
int writeCompHashTableExtTabLiterals(writeCompHashTableCtx_t* ctx, extend_hit_t* recs, uint64_t numRecs)
{
  uint64_t  cacheWord = 0;
  int       cacheBits = 0;
  int       recBits;
  uint64_t  i, oldBits;
  uint64_t  rec;
  uint32_t  run           = 0;
  uint64_t* blockLenPtr   = NULL;
  uint64_t* blockEndPtr   = NULL;
  int64_t   blockStartBit = 0;
  // Flush buffer so all room is available
  PUT_FLUSH;
  for (i = 0; 1; i++) {
    // Insert block boundary when the buffer nears full
    if (!i || ctx->bufBytes + (((ctx->bitLen & 7) + cacheBits) >> 3) >= BLOCK_LIM || i == numRecs) {
      // At end, encode remaining run length
      if (i == numRecs) {
        oldBits = ctx->bitLen + cacheBits;
        if (run) {
          for (; run >= COMP_EXT_NOLITERAL_RUN_LIMIT; run -= COMP_EXT_NOLITERAL_RUN_LIMIT)
            PUT_BITS(COMP_EXT_NOLITERAL_RUN_LIMIT, COMP_EXT_NOLITERAL_RUN_BITS);
          PUT_BITS(run, COMP_EXT_NOLITERAL_RUN_BITS);
        }
        ctx->miscBits += ctx->bitLen + cacheBits - oldBits;
      }
      // Align to byte boundary
      PUT_ALIGN;
      // Fill in previous block's end record and coded bit length
      if (i) {
        *blockEndPtr = i + ctx->extLitPos;
        *blockLenPtr = ctx->bitLen - blockStartBit;
      }
      // Flush buffer so all room is available
      PUT_FLUSH;
      // Done?
      if (i == numRecs) break;
      // Magic number for synchronization
      PUT_BITS(MAGIC_LIT_BLOCK, 32);
      ctx->miscBits += 32;
      // Block start record
      PUT_QWORD(i + ctx->extLitPos);
      ctx->miscBits += 64;
      // Placeholder for block end record
      blockEndPtr = (uint64_t*)&ctx->buf[(ctx->bitLen - ctx->bufStartBit) >> 3];
      PUT_QWORD(0);
      ctx->miscBits += 64;
      // Placeholder for block encoded length
      blockLenPtr = (uint64_t*)&ctx->buf[(ctx->bitLen - ctx->bufStartBit) >> 3];
      PUT_QWORD(0);
      ctx->miscBits += 64;
      blockStartBit = ctx->bitLen;
      // Flag indicating block type is extension table literals
      PUT_BITS(1, 1);
    }
    // Just tally records without literal flags
    if (!recs[i].literal) {
      run++;
      continue;
    }
    // Clear the literal flags before possible uncompressed output
    recs[i].literal = 0;
    // Encode run length before the next literal record.  Repeat the maximum length until the remainder is
    // less.
    oldBits = ctx->bitLen + cacheBits;
    for (; run >= COMP_EXT_NOLITERAL_RUN_LIMIT; run -= COMP_EXT_NOLITERAL_RUN_LIMIT) {
      PUT_BITS(COMP_EXT_NOLITERAL_RUN_LIMIT, COMP_EXT_NOLITERAL_RUN_BITS);
    }
    PUT_BITS(run, COMP_EXT_NOLITERAL_RUN_BITS);
    run = 0;
    rec = ((uint64_t*)(recs))[i];
    // Discard lift_group field if lift_code is null
    recBits = (recs[i].lift_code == LIFT_CODE_NONE ? 35 : 63);
    PUT_BITS(rec, recBits);
    ctx->extLitBits += ctx->bitLen + cacheBits - oldBits;
    ctx->extLitRecs++;
  }

  // Update the position in the literal records stream
  ctx->extLitPos += numRecs;
  return 0;
}

// Writes an index of the extension table, indicating the number of extend table records corresponding
// to each bin of hash table buckets, to facilitate using short relative offsets to place extend table records
int writeCompHashTableExtIndex(
    writeCompHashTableCtx_t* ctx, extIndexRec_t* extIndexRecs, uint32_t numExtIndexRecs, uint32_t extTabRecs)
{
  uint64_t cacheWord = 0;
  size_t   cacheBits = 0;
  size_t   i;
  // Remember the extension table length
  ctx->extTabRecs = extTabRecs;
  // Magic number to verify synchronization
  PUT_BITS(MAGIC_LIT_END, 32);
  ctx->miscBits += 32;
  // Magic number to verify synchronization
  PUT_BITS(MAGIC_EXT_IDX_START, 32);
  ctx->miscBits += 32;
  // Length of extend table index
  PUT_BITS(numExtIndexRecs, 32);
  ctx->miscBits += 32;
  // Extend table index length fields
  for (i = 0; i < numExtIndexRecs; i++) PUT_BITS(extIndexRecs[i].length, 32);
  ctx->autoSecBits = numExtIndexRecs * 32;
  // Magic number to verify synchronization
  PUT_BITS(MAGIC_EXT_IDX_END, 32);
  ctx->miscBits += 32;
  // Flush
  PUT_FLUSH;
  return 0;
}

// Routine to encode "automatic" hit records in the compressed file
int writeCompHashTableAutoHits(
    writeCompHashTableCtx_t* ctx,
    seedPopRec_t*            popRecs,
    uint32_t                 numPopRecs,
    uint32_t                 maxExtendIds,
    extIndexRec_t*           extIndexRecs)
{
  uint64_t     cacheWord = 0;
  int          cacheBits = 0;
  seedPopRec_t popRec;
  uint64_t     oldBits;
  uint32_t     i, offset, code, len, extendIdBits, liftGroup, liftGroupBits, liftGroupBitsCode, offsetBits;
  for (extendIdBits = 0; maxExtendIds; maxExtendIds >>= 1, extendIdBits++)
    ;
  ctx->extendIdBits = extendIdBits;
  // Magic number to verify synchronization
  PUT_BITS(MAGIC_AUTO_START, 32);
  ctx->miscBits += 32;
  // Put maximum extend ID bits in the archive
  PUT_BITS(extendIdBits, 32);
  ctx->miscBits += 32;
  uint64_t* blockLenPtr   = NULL;
  uint64_t* blockEndPtr   = NULL;
  int64_t   blockStartBit = 0;
  int       liftCodePri;
  // Flush buffer so all room is available
  PUT_FLUSH;
  // Loop through automatic seed population records
  for (i = 0; 1; i++) {
    // Insert block boundary when the buffer nears full
    if (!i || ctx->bufBytes + (((ctx->bitLen & 7) + cacheBits) >> 3) >= BLOCK_LIM || i == numPopRecs) {
      // Align to byte boundary
      PUT_ALIGN;
      // Fill in previous block's end record and coded bit length
      if (i) {
        *blockEndPtr = i;
        *blockLenPtr = ctx->bitLen - blockStartBit;
      }
      // Flush buffer so all room is available
      PUT_FLUSH;
      // Done?
      if (i == numPopRecs) break;
      // Magic number for synchronization
      PUT_BITS(MAGIC_AUTO_BLOCK, 32);
      ctx->miscBits += 32;
      // Block start record
      PUT_QWORD(i);
      ctx->miscBits += 64;
      // Placeholder for block end record
      blockEndPtr = (uint64_t*)&ctx->buf[(ctx->bitLen - ctx->bufStartBit) >> 3];
      PUT_QWORD(0);
      ctx->miscBits += 64;
      // Placeholder for block encoded length
      blockLenPtr = (uint64_t*)&ctx->buf[(ctx->bitLen - ctx->bufStartBit) >> 3];
      PUT_QWORD(0);
      blockStartBit = ctx->bitLen;
    }
    popRec  = popRecs[i];
    offset  = popRec.slotOffset;
    oldBits = ctx->bitLen + cacheBits;
    // Record type prefix
    if (!popRec.automatic) {
      PUT_BITS(3, 2);  // No record
      ctx->autoNulBits += ctx->bitLen + cacheBits - oldBits;
      ctx->autoNulRecs++;
      continue;
    }
    if (!popRec.extTable) {
      PUT_BITS(0, 1);  // Primary
      // Variable-length offset code:
      // Prefix run of '1's indicating offset bit-length;
      // Offset bits follow, with MSB dropped except for bit-length 1
      //  1-bit offset (0-1)   : length=2,  code=0x
      //  2-bit offset (2-3)   : length=3,  code=10x
      //  3-bit offset (4-7)   : length=5,  code=110xx
      //  4-bit offset (8-15)  : length=7,  code=1110xxx
      //  5-bit offset (16-31) : length=9,  code=11110xxxx
      //  6-bit offset (32-63) : length=11, code=111110xxxxx
      //  ...
      code = offset << 1;
      len  = 1;
      for (offset >>= 1; offset; offset >>= 1) {
        code = (code << 1) | 1;
        len += 2;
      }
      len += (len == 1);  // Retain the MSB for length 1
      PUT_BITS(code, len);
      ctx->autoPriBits += ctx->bitLen + cacheBits - oldBits;
      ctx->autoPriRecs++;
    } else {
      PUT_BITS(1, 2);  // Secondary
      if (popRec.liftCode == LIFT_CODE_NONE) {
        PUT_BITS(0, 1);  // No liftover
      } else {
        PUT_BITS(1, 1);  // Liftover
        liftCodePri = (popRec.liftCode == LIFT_CODE_PRI);
        PUT_BITS(liftCodePri, 1);  // Lift code is primary (not alt)
        // Measure the bit-length of liftGroup
        liftGroup     = popRec.liftGroup;
        liftGroupBits = 0;
        for (; liftGroup; liftGroup = liftGroup >> 1) liftGroupBits++;
        // Don't want to use more than 4 bits to encode the bit-length of liftGroup,
        // so let a code of 15 represent the maximum number of bits
        if (liftGroupBits >= 15) {
          liftGroupBits     = SEED_POP_REC_LIFT_GROUP_BITS;
          liftGroupBitsCode = 15;
        } else {
          liftGroupBitsCode = liftGroupBits;
        }
        PUT_BITS(liftGroupBitsCode, 4);  // LiftGroup bits
        // Drop the liftGroup MSB '1' for lengths 1-14, because it's implied by the length
        liftGroupBits -= (liftGroupBits > 0 && liftGroupBits < 15);
        PUT_BITS(popRec.liftGroup, liftGroupBits);  // LiftGroup
      }
      // Look up the number of offset bits to use for this bucket bin of the extension table
      offsetBits = extIndexRecs[popRec.bucketBin].offsetBits;
      PUT_BITS(offset, offsetBits);  // Extension table offset relative to hash bucket bin
      ctx->autoSecBits += ctx->bitLen + cacheBits - oldBits;
      ctx->autoSecRecs++;
    }
  }
  // Magic number for synchronization
  PUT_BITS(MAGIC_AUTO_END, 32);
  ctx->miscBits += 32;
  PUT_FLUSH;
  // Stats
  ctx->totalRecs = ctx->autoNulRecs + ctx->autoPriRecs + ctx->autoSecRecs + ctx->specialRecs +
                   ctx->literalRecs + ctx->chainPtrRecs + ctx->chainEndRecs;
  ctx->totalBits = ctx->autoNulBits + ctx->autoPriBits + ctx->autoSecBits + ctx->specialBits +
                   ctx->literalBits + ctx->chainPtrBits + ctx->chainEndBits;
  return 0;
}

uint64_t lclGetBits(decompHashTableCtx_t* ctx, int len)
{
  uint64_t val = 0, mask = 0xFFFFFFFFFFFFFFFFULL >> (64 - len);
  int64_t  bitPos = ctx->bitPos;
  int64_t  dist   = ctx->bufBits - ctx->bitPos;
  int      s;
  ctx->bitPos += len;
  if (dist > 64)
    val = *((uint64_t*)&ctx->buf[bitPos >> 3]);
  else if (dist < len)
    return 0;
  else if (dist > 0)
    memcpy(&val, &ctx->buf[bitPos >> 3], (dist + 7) >> 3);
  s = bitPos & 7;
  val >>= s;
  s = 64 - s;
  if (len > s) val |= ((uint64_t)ctx->buf[ctx->bitPos >> 3] << s);

  return val & mask;
}

#define GET_BIT                  \
  (ctx->bitPos < ctx->bufBits && \
   (ctx->oldBitPos = ctx->bitPos++, (ctx->buf[ctx->oldBitPos >> 3] & (1 << (ctx->oldBitPos & 7)))))
//#define GET_BIT (ctx->bitPos<ctx->bufBits&&(ctx->buf[ctx->bitPos>>3]&(1<<(ctx->bitPos&7))), ctx->bitPos++)
#define GET_BITS(len) (lclGetBits(ctx, len))
#define GET_ALIGN (lclGetBits(ctx, (8 - (ctx->bitPos & 7)) & 7))
#define GET_ERROR (ctx->bitPos > ctx->bufBits)

// Routine to initialize the context record for decompressing a hash table
char* decompHashTableCtxInit(
    decompHashTableCtx_t* ctx,
    int                   threads,
    uint8_t*              compBuf,
    int64_t               compLen,
    uint8_t*              refBuf,
    int64_t               refLen,
    uint8_t*              decompBuf,
    int64_t               decompLen,
    uint8_t*              extendTableBuf,
    int64_t               extendTableLen)
{
  memset(ctx, 0, sizeof(decompHashTableCtx_t));
  if (!compBuf) return "Decompress called with null compressed buffer";
  if (compLen <= 0) return "Decompress called with nonpositive length buffer";
  if (!refBuf) return "Decompress called with null reference buffer";
  if (refLen <= 0) return "Decompress called with nonpositive length reference";
  ctx->threads     = threads;
  ctx->buf         = compBuf;
  ctx->bufBytes    = compLen;
  ctx->bufBits     = compLen << 3;
  ctx->hashTable   = (hashrec_t*)decompBuf;
  ctx->hashRecs    = decompLen / HASH_RECORD_BYTES;
  ctx->extendTable = (extend_hit_t*)extendTableBuf;
  ctx->extendRecs  = extendTableLen / sizeof(extend_hit_t);
  ctx->refBuf2Bit  = malloc(refLen / 2);
  if (!ctx->refBuf2Bit) return "Decompress failed to allocate reference buffer";
  ctx->refBases = refLen * 2;
  // Make table for translating 4-bit to 2-bit reference codes
  // (only needs to work for single-base input codes 1, 2, 4, 8)
  uint8_t table4to2[256] = {0};
  int     i, j;
  for (i = 0; i < 16; i++)
    for (j = 0; j < 16; j++)
      table4to2[i | j << 4] =
          (i == 8 ? 3 : i == 4 ? 2 : i == 2 ? 1 : 0) | (j == 8 ? 3 : j == 4 ? 2 : j == 2 ? 1 : 0) << 2;
  // Translate the 4-bit reference input to 2-bit bases
  int64_t  n;
  uint8_t *src = refBuf, *dst = ctx->refBuf2Bit;
  for (n = 0; n < refLen; n += 2) {
    *dst = table4to2[*src++];
    *dst++ |= table4to2[*src++] << 4;
  }
  return NULL;
}

// Routine to decode the header of a compressed hash table
char* decompHashTableHeader(decompHashTableCtx_t* ctx)
{
  uint8_t* p = (uint8_t*)&ctx->cfgHdr;
  int      i;
  // Magic number for this file type
  if (GET_BITS(32) != MAGIC_FILE_START) return "Compressed hash table format wrong at start";
  // Compressed format version
  if (GET_BITS(32) != COMP_VERSION) return "Compressed hash table format version wrong";
  // Entire configuration header
  for (i = 0; i < sizeof(hashTableHeader_t); i++) *p++ = GET_BITS(8);
  // Magic number to verify synchronization
  if (GET_BITS(32) != MAGIC_AFTER_HEADER) return "Compressed hash table format wrong after header";
  // Make sure config header's extTabRecs field isn't uninitialized
  if (ctx->cfgHdr.extTabRecs == 0xFFFFFFFF) return "Decompressed header extTabRecs field is uninitialized";
  ctx->priCrcInit = crcHash64Init(ctx->cfgHdr.priCrcBits, ctx->cfgHdr.priCrcPoly);
  ctx->secCrcInit = crcHash64Init(ctx->cfgHdr.secCrcBits, ctx->cfgHdr.secCrcPoly);
  // Output buffer provided?
  if (ctx->hashTable) {
    // Verify expected hash table length
    if (ctx->cfgHdr.hashTableBytes / HASH_RECORD_BYTES != ctx->hashRecs)
      return "Decompressed hash table buffer length doesn't match expected";
  }
  // Or allocate our own
  else {
    ctx->hashTable = malloc(ctx->cfgHdr.hashTableBytes);
    if (!ctx->hashTable) return "Failed to allocate buffer for uncompressed hash table";
    ctx->hashRecs = ctx->cfgHdr.hashTableBytes / HASH_RECORD_BYTES;
  }
  // Extend table buffer provided?
  if (ctx->extendTable) {
    // Verify expected hash table length
    if (ctx->cfgHdr.extTabRecs != ctx->extendRecs)
      return "Decompressed extend table buffer length doesn't match expected";
  }
  // Or allocate our own
  else {
    ctx->extendTable = malloc((size_t)ctx->cfgHdr.extTabRecs * sizeof(extend_hit_t));
    if (!ctx->extendTable) return "Failed to allocate buffer for uncompressed extend table";
    ctx->extendRecs = ctx->cfgHdr.extTabRecs;
  }
  return GET_ERROR ? "Unexpected end of compressed hash table after header" : NULL;
}

typedef union {
  extend_hit_t extend_hit;
  uint64_t     qword;
} extend_hit_union_t;

// Thread function for populating literal records
void* literalThread(void* sharePtr)
{
  mthreadShare_t*       share     = sharePtr;
  int                   threadId  = mthreadGetThreadId(share);
  decompHashTableCtx_t* masterCtx = share->ctx;
  decompHashTableCtx_t  localCtx  = *masterCtx;
  decompHashTableCtx_t* ctx       = &localCtx;
  blockWorkRec_t*       work      = NULL;
  int64_t               totRun, run;
  hashrec_t             hashRec;
  extend_hit_union_t    extHitUnion;
  uint64_t *            rp, *recEnd, emptyRec;
  uint32_t              blockTypeExt, runBits, runLimit;
  // Flag chain terminators temporarily
  hashrec_t termRec;
  termRec.qword           = HASH_REC_CHAIN_TERM_QWORD;
  termRec.chain.chain_pad = 1;

  while (1) {
    // Get more work
    work = mthreadChildGetJob(threadId, share);
    if (!work) return NULL;
    // Set the bit cursor and end barrier in the local copy of the decompression context
    ctx->bitPos  = work->startBit;
    ctx->bufBits = work->endBit;
    // Block type flag
    blockTypeExt = GET_BIT;
    if (blockTypeExt) {
      // Set the start and end record pointers in the extension table
      rp     = (uint64_t*)(&ctx->extendTable[work->startPos]);
      recEnd = (uint64_t*)(&ctx->extendTable[work->endPos]);
      // Other mode parameters
      emptyRec = 0;
      runBits  = COMP_EXT_NOLITERAL_RUN_BITS;
      runLimit = COMP_EXT_NOLITERAL_RUN_LIMIT;
    } else {
      // Set the start and end record pointers in the hash table
      rp     = (uint64_t*)(&ctx->hashTable[work->startPos]);
      recEnd = (uint64_t*)(&ctx->hashTable[work->endPos]);
      // Other mode parameters
      emptyRec = HASH_REC_EMPTY_QWORD;
      runBits  = COMP_NOLITERAL_RUN_BITS;
      runLimit = COMP_NOLITERAL_RUN_LIMIT;
    }
    while (1) {
      if (GET_ERROR) {
        mthreadSetErrMsg(share, "Unexpected end of compressed block inside literals section");
        break;
      }
      if (rp >= recEnd) break;
      // Decode run length before the next literal record.  The maximum length repeats until the remainder is
      // less.
      totRun = 0;
      do {
        totRun += run = GET_BITS(runBits);
      } while (run == runLimit);
      if (rp + totRun > recEnd) return "Compressed non-literal run overran hash table end";
      // Populate empty records for the run length
      while (totRun--) *rp++ = emptyRec;
      if (rp >= recEnd) break;
      // Decode literal extension table record
      if (blockTypeExt) {
        // First read record without liftGroup field
        extHitUnion.qword = GET_BITS(35);
        // On active liftCode, also read liftGroup
        if (extHitUnion.extend_hit.lift_code != LIFT_CODE_NONE)
          extHitUnion.extend_hit.lift_group = GET_BITS(28);
        // Write into extension table
        *rp++ = extHitUnion.qword;
      }
      // Decode literal hash table record
      else {
        // Chain flag
        if (GET_BIT) {
          // Term flag
          if (GET_BIT) {
            *rp++ = termRec.qword;
            continue;
          }
          // Non-terminator
          hashRec.qword           = 0;
          hashRec.chain.opcode    = GET_BIT ? HASH_OPC_CHAIN_CON_MASK : HASH_OPC_CHAIN_BEG_MASK;
          hashRec.chain.chain_ptr = GET_BITS(18);
          *rp++                   = hashRec.qword;
          continue;
        }
        // Non-chain
        *rp++ = GET_BITS(64);
      }
    }
    // Align to byte boundary
    GET_ALIGN;
    if (ctx->bitPos != work->endBit)
      mthreadSetErrMsg(share, "Decompression populating literals, block ended at wrong bit position");
    // Return work batch
    mthreadChildPutJob(threadId, share, work);
  }
}

// Fetch 32 bases in a uint64_t
#define GET_SEQ_BITS(pos)                                        \
  (((*((uint64_t*)(&refSeq[(pos) >> 2]))) >> (((pos)&3) << 1)) | \
   ((uint64_t)(refSeq[((pos) >> 2) + 8] << (8 - (((pos)&3) << 1)))) << 56)
// Thread function for populating automatic hit records
void* autoHitThread(void* sharePtr)
{
  mthreadShare_t*       share     = sharePtr;
  int                   threadId  = mthreadGetThreadId(share);
  decompHashTableCtx_t* masterCtx = share->ctx;
  decompHashTableCtx_t  localCtx  = *masterCtx;
  decompHashTableCtx_t* ctx       = &localCtx;
  blockWorkRec_t*       work      = NULL;
  uint8_t*              refSeq    = ctx->refBuf2Bit;
  uint32_t              len, sec;
  hashrec_t             rec, *slot, *bucket;
  extend_hit_t          extHit, *extTabPtr, *extendTable = ctx->extendTable;
  extIndexRec_t         extIndexRec, *extIndexRecs = ctx->extIndexRecs;
  uint64_t              pos, fwSeed, rcSeed, hashKey, hash, bucketAddr, bucketIndex;
  uint64_t refBinMask = ctx->cfgHdr.anchorBinBits ? ((0ULL - 1) << ctx->cfgHdr.anchorBinBits) : 0;
  uint32_t seedIndex, extTabPos, msb, liftGroupBits, liftGroupBitsCode;
  double   refSeedInterval = ctx->cfgHdr.refSeedInterval;
  double   squeezeRatio    = (double)ctx->cfgHdr.tableSize64ths / 64;  // This floating point value is exact
  int      seedLen         = ctx->cfgHdr.priSeedBases;
  int      useRc, offset, chain;
  uint64_t seedMask    = ((uint64_t)1 << (seedLen << 1)) - 1;
  uint32_t wrapBuckets = (uint32_t)(MAX_WRAP_BYTES * squeezeRatio) >> HASH_BUCKET_BYTES_LOG2;

  rcSeed = 0;
  PRINTIF(DEBUG_THREAD, ("Auto-hit worker thread %d started\n", threadId));
  // Work batch loop
  while (1) {
    // Get more work
    work = mthreadChildGetJob(threadId, share);
    if (!work) return NULL;
    // Set the bit cursor and end barrier in the local copy of the decompression context
    ctx->bitPos  = work->startBit;
    ctx->bufBits = work->endBit;
    // Loop through automatic seed population records
    for (seedIndex = work->startPos; seedIndex < work->endPos; seedIndex++) {
      if (GET_ERROR) {
        mthreadSetErrMsg(share, "Unexpected end of compressed block inside automatic section");
        break;
      }
      // Secondary or no record?
      if ((sec = GET_BIT) && GET_BIT) continue;
      // Calculate seed position
      pos = floor(seedIndex * refSeedInterval);
      // First hash the primary seed, even for a secondary record
      // Fetch seed
      fwSeed = GET_SEQ_BITS(pos) & seedMask;
      // Compute reverse complement seed
      revComp(&rcSeed, &fwSeed, seedLen);
      // Use whichever is numerically smaller as the seed
      useRc   = (rcSeed < fwSeed);
      hashKey = useRc ? rcSeed : fwSeed;
      // Append the reference bin if applicable
      hashKey |= (pos & refBinMask) << KEY_ANCHOR_OFFSET;
      // Hash the seed
      crcHash64(ctx->priCrcInit, &hashKey, &hash);
      PRINTIF(DEBUG_AUTO, ("----------------------------------------\n"));
      PRINTIF(DEBUG_AUTO, ("bitPos=%lld, pos=%lld, sec=%d\n", ctx->bitPos, pos, sec));
      PRINTIF(
          DEBUG_AUTO, ("fwSeed=0x%llX, rcSeed=0x%llX, useRc=%d, hash=0x%llX\n", fwSeed, rcSeed, useRc, hash));
      // Calculate base bucket address
      bucketAddr  = (hash >> HASH_BYTE_ADDR_START) * squeezeRatio;
      bucketIndex = bucketAddr >> HASH_BUCKET_BYTES_LOG2;
      PRINTIF(
          DEBUG_AUTO,
          ("bucketIndex=0x%llX, byteAddr=0x%llX\n", bucketIndex, bucketIndex << HASH_BUCKET_BYTES_LOG2));
      // Secondary, in seed extension table...
      if (sec) {
        // Fill hit information into extension table record
        extHit.pos = seedIndex;
        extHit.rc  = useRc;
        // No liftover?
        if (!GET_BIT) {
          extHit.lift_code  = LIFT_CODE_NONE;
          extHit.lift_group = 0;
        }
        // Liftover
        else {
          // Lift code
          extHit.lift_code = GET_BIT ? LIFT_CODE_PRI : LIFT_CODE_ALT;
          // Lift group bits code 0-15
          liftGroupBitsCode = GET_BITS(4);
          // Inferred MSB '1' on bits code 1-14
          msb = (liftGroupBitsCode > 0 && liftGroupBitsCode < 15);
          // Bits code 15 means the maximum bit length; otherwise deduct bit for inferred MSB
          liftGroupBits = (liftGroupBitsCode == 15 ? SEED_POP_REC_LIFT_GROUP_BITS : liftGroupBitsCode - msb);
          // Append lower bits
          extHit.lift_group = (msb << liftGroupBits) | (liftGroupBits ? GET_BITS(liftGroupBits) : 0);
          PRINTIF(DEBUG_AUTO, ("liftGroupBitsCode=%u, liftGroupBits=%u\n", liftGroupBitsCode, liftGroupBits));
        }
        PRINTIF(DEBUG_AUTO, ("liftCode=%u, liftGroup=%u\n", extHit.lift_code, extHit.lift_group));
        // Look up the hash bucket's larger bin in the extension table index
        extIndexRec = extIndexRecs[bucketIndex >> EXTTAB_INDEX_BUCKET_BITS];
        // The index entry indicates how many bits long offset values are for this bucket bin.
        // Read an offset of that bit-length, and add to the bucket bin's start position
        // for a net pointer into the extension table.
        offset    = extIndexRec.offsetBits ? GET_BITS(extIndexRec.offsetBits) : 0;
        extTabPos = extIndexRec.start + offset;
        extTabPtr = &extendTable[extTabPos];
        PRINTIF(
            DEBUG_AUTO,
            ("bucketBin=%u, start=%u, length=%u, offsetBits=%u, offset=%d\n",
             bucketIndex >> EXTTAB_INDEX_BUCKET_BITS,
             extIndexRec.start,
             extIndexRec.length,
             extIndexRec.offsetBits,
             offset));
        PRINTIF(DEBUG_AUTO, ("Populating at extTabPos=%u\n", extTabPos));
        if (extTabPtr->pos) {
          mthreadSetErrMsg(
              share, "Decompression populating automatic hit, target extend table slot nonempty");
          break;
        }
        // Populate
        *extTabPtr = extHit;
        continue;
      }
      // Else primary, in hash table...
      //
      // Variable-length offset code:
      // Prefix run of '1's indicating offset bit-length;
      // Offset bits follow, with MSB dropped except for bit-length 1
      //  1-bit offset (0-1)   : length=2,  code=0x
      //  2-bit offset (2-3)   : length=3,  code=10x
      //  3-bit offset (4-7)   : length=5,  code=110xx
      //  4-bit offset (8-15)  : length=7,  code=1110xxx
      //  5-bit offset (16-31) : length=9,  code=11110xxxx
      //  6-bit offset (32-63) : length=11, code=111110xxxxx
      // Could terminate here and drop the '0' for 6-bit offset, but too rare to care, so keep extensible
      for (len = 0; GET_BIT; len++)
        ;
      offset = GET_BITS(len + !len) | ((len != 0) << len);
      PRINTIF(DEBUG_AUTO, ("slotOffset=%u\n", offset));
      // Clear hash record
      rec.qword = 0;
      // Pointer to hash record at bucket address
      bucket = &ctx->hashTable[bucketIndex << HASH_RECORDS_PER_BUCKET_LOG2];
      // Fill in the HIT record
      rec.hit.pos       = seedIndex;
      rec.hit.rc        = useRc;
      rec.hit.hash_bits = hash;
      rec.hit.thread_id = bucketAddr >> ADDR_THREAD_ID_START;
      // Chaining if chain-begin record in the target bucket's last slot
      if (offset >= HASH_RECORDS_PER_BUCKET - 1 &&
          HASH_OPC_IS_BEG(bucket[HASH_RECORDS_PER_BUCKET - 1].general.opcode)) {
        // Deduct initial bucket records before chaining
        offset -= HASH_RECORDS_PER_BUCKET - 1;
        // Chain pointer slot
        slot = bucket + HASH_RECORDS_PER_BUCKET - 1;
        PRINTIF(DEBUG_AUTO, ("Chain rec in slot 7: %016llX, remaining offset=%d\n", slot->qword, offset));
        while (1) {
          // Follow chain pointer
          bucketIndex = (bucketIndex & ((uint64_t)0 - CHAIN_BLOCK_BUCKETS)) | slot->chain.chain_ptr;
          bucket      = &ctx->hashTable[bucketIndex << HASH_RECORDS_PER_BUCKET_LOG2];
          PRINTIF(
              DEBUG_AUTO,
              ("Chained to bucketIndex=0x%llX, byteAddr=0x%llX\n",
               bucketIndex,
               bucketIndex << HASH_BUCKET_BYTES_LOG2));
          // Scan next bucket for chain continuation record
          for (chain = 0; chain < HASH_RECORDS_PER_BUCKET && !HASH_OPC_IS_CON(bucket[chain].general.opcode);
               chain++)
            ;
          if (chain == HASH_RECORDS_PER_BUCKET) {
            PRINTIF(DEBUG_AUTO, ("Missing chain continuation record\n"));
            mthreadSetErrMsg(
                share, "Decompression populating automatic hit, missing chain continuation record");
            goto autoHitThreadWorkFail;
          }
          // Deduct the bucket's chain capacity from remaining offset
          offset -= (HASH_RECORDS_PER_BUCKET - 1) - chain;
          PRINTIF(
              DEBUG_AUTO,
              ("Chain rec in slot %d: %016llX, remaining offset=%d\n", chain, bucket[chain].qword, offset));
          // Negative: stop here
          if (offset < 0) {
            slot = bucket + HASH_RECORDS_PER_BUCKET + offset;
            break;
          }
          // Keep chaining
          slot = &bucket[chain];
          if (slot->chain.chain_pad) {
            PRINTIF(DEBUG_AUTO, ("Unexpected chain terminator\n"));
            mthreadSetErrMsg(share, "Decompression populating automatic hit, unexpected chain terminator");
            goto autoHitThreadWorkFail;
          }
        }
      }
      // Probing
      else {
        while (offset >= HASH_RECORDS_PER_BUCKET) {
          // Deduct size of this bucket
          offset -= HASH_RECORDS_PER_BUCKET;
          // Next bucket, wrapping if necessary
          if (++bucketIndex % wrapBuckets == 0) bucketIndex -= wrapBuckets;
          bucket = &ctx->hashTable[bucketIndex << HASH_RECORDS_PER_BUCKET_LOG2];
          PRINTIF(
              DEBUG_AUTO,
              ("Probed to bucketIndex=0x%llX, byteAddr=0x%llX\n",
               bucketIndex,
               bucketIndex << HASH_BUCKET_BYTES_LOG2));
        }
        slot = bucket + offset;
      }
      // Populate in the target slot
      PRINTIF(
          DEBUG_AUTO,
          ("Populating in bucketIndex=0x%llX, slot %d, byteAddr=0x%llX\n",
           bucketIndex,
           slot - bucket,
           (slot - ctx->hashTable) << HASH_RECORD_BYTES_LOG2));
      if (slot->general.opcode != HASH_OPC_EMPTY) {
        mthreadSetErrMsg(share, "Decompression populating automatic hit, target slot nonempty");
        break;
      }
      *slot = rec;
    }
    GET_ALIGN;
    if (ctx->bitPos != work->endBit)
      mthreadSetErrMsg(share, "Decompression populating automatic hits, block ended at wrong bit position");
  autoHitThreadWorkFail:
    // Return work batch
    mthreadChildPutJob(threadId, share, work);
  }
}

#define BLOCK_TYPE_LITERAL 1
#define BLOCK_TYPE_AUTO 2

// Routine to decode blocks and populate hash table records from the compressed file
char* decompHashTableBlocks(decompHashTableCtx_t* ctx, int blockType)
{
  blockWorkRec_t* work = NULL;
  uint32_t        blockMagic, endMagic, magic;
  int64_t         totalRecs;
  int64_t         bitLen, numPos;
  char*           err            = NULL;
  void* (*threadFunction)(void*) = NULL;

  switch (blockType) {
  case BLOCK_TYPE_LITERAL:
    threadFunction = literalThread;
    totalRecs      = ctx->hashRecs;
    blockMagic     = MAGIC_LIT_BLOCK;
    endMagic       = MAGIC_LIT_END;
    break;
  case BLOCK_TYPE_AUTO:
    threadFunction = autoHitThread;
    totalRecs      = (uint32_t)ceil(ctx->cfgHdr.refSeqLen / ctx->cfgHdr.refSeedInterval);
    blockMagic     = MAGIC_AUTO_BLOCK;
    endMagic       = MAGIC_AUTO_END;
    // Magic number to verify synchronization
    if (GET_BITS(32) != MAGIC_AUTO_START) return "Compressed hash table format wrong at automatic section";
    // Extension ID bits from the archive
    ctx->extendIdBits = GET_BITS(32);
    break;
  default:
    return "Decompress routine called with invalid block type";
  }

  // Initialize shared structure for muti-threading
  mthreadShare_t* share = mthreadInit(ctx, threadFunction, ctx->threads, sizeof(blockWorkRec_t));

  while (1) {
    // Fetch an empty work batch
    work = mthreadParentGetJob(share);
    if (!work || share->errMsg) {
      err = share->errMsg;
      if (!err) err = "Null work pointer returned by decompression thread";
      goto decompHashTableBlocksError;
    }
    // Magic number to verify synchronization of each block
    magic = GET_BITS(32);
    // Done marker
    if (magic == endMagic) break;
    if (magic != blockMagic) return "Compressed hash table format wrong at block boundary";
    // Read the next block header, populating the work unit
    work->startPos = GET_BITS(64);
    work->endPos   = GET_BITS(64);
    bitLen         = GET_BITS(64);
    numPos         = work->endPos - work->startPos;
    work->startBit = ctx->bitPos;
    work->endBit   = ctx->bitPos + bitLen;
    PRINTIF(
        DEBUG_THREAD,
        ("Next block: startPos=%lld, endPos=%lld, numPos=%lld, startBit=%lld, endBit=%lld, bitLen=%lld, byteLen=%g\n",
         work->startPos,
         work->endPos,
         numPos,
         work->startBit,
         work->endBit,
         bitLen,
         bitLen / 8.0));
    if (work->endBit > ctx->bufBits) {
      err = "Block endpoint leaves the compressed buffer";
      goto decompHashTableBlocksError;
    }
    if (work->endPos > totalRecs) {
      err = "Compressed block end position too high";
      goto decompHashTableBlocksError;
    }
    // Skip over the rest of the block
    ctx->bitPos += bitLen;
    // Submit work batch
    mthreadParentPutJob(share, work);
  }

  // Wait for worker threads to terminate
  err   = mthreadClose(share);
  share = NULL;
  if (err) goto decompHashTableBlocksError;
  // Advance to byte boundary
  GET_ALIGN;
  if (GET_ERROR) {
    err = "Unexpected end of compressed hash table inside automatic section";
    goto decompHashTableBlocksError;
  }
decompHashTableBlocksError:
  mthreadCancel(share);
  return err;
}

// Routine to decode the literal records section from a compressed hash table and populate them
char* decompHashTableLiterals(decompHashTableCtx_t* ctx)
{
  return decompHashTableBlocks(ctx, BLOCK_TYPE_LITERAL);
}

// Routine to read the extension table index
char* decompHashTableExtIndex(decompHashTableCtx_t* ctx)
{
  size_t i;
  // Magic number to verify synchronization
  if (GET_BITS(32) != MAGIC_EXT_IDX_START)
    return "Compressed hash table format wrong at extend table index start";
  // Length of extend table index
  uint64_t numExtIndexRecs = GET_BITS(32);
  if ((numExtIndexRecs << (HASH_BUCKET_BYTES_LOG2 + EXTTAB_INDEX_BUCKET_BITS)) != ctx->cfgHdr.hashTableBytes)
    return "Compressed extend table index doesn't match hash table size";
  // Allocate extend table index
  if (!(ctx->extIndexRecs = malloc(numExtIndexRecs * sizeof(extIndexRec_t))))
    return "Unable to allocate extend table index during decompress";
  // Build extend table index
  extIndexRec_t* p     = ctx->extIndexRecs;
  uint32_t       start = 0, length, maxOffset, offsetBits;
  for (i = 0; i < numExtIndexRecs; i++, p++) {
    // Decode length of the extend table section corresponding to this index entry / bin of hash buckets
    length = GET_BITS(32);
    // Derive number of bits needed to represent relative offsets in this bin, less than this length
    maxOffset = length ? length - 1 : 0;
    for (offsetBits = 0; maxOffset; maxOffset = maxOffset >> 1) offsetBits++;
    // Save the index record
    p->start      = start;
    p->length     = length;
    p->offsetBits = offsetBits;
    // Accumulate the start position within the extend table for the next index entry
    start += length;
  }
  // Magic number to verify synchronization
  if (GET_BITS(32) != MAGIC_EXT_IDX_END)
    return "Compressed hash table format wrong at extend table index end";
  return GET_ERROR ? "Unexpected end of compressed hash table after header" : NULL;
}

// Routine to decode "automatic" hit records from the compressed file, regenerate and populate in the hash
// table
char* decompHashTableAutoHits(decompHashTableCtx_t* ctx)
{
  return decompHashTableBlocks(ctx, BLOCK_TYPE_AUTO);
}

// Context for threads populating automatic hit records
typedef struct {
  int                   threadId;
  int64_t               startBucket;
  int64_t               endBucket;
  char*                 errMsg;
  decompHashTableCtx_t* ctx;
} flagThreadCtx_t;

// Routine to annotate last flags and chain filters into the decompressed hash table
void* flagsThread(void* ctxPtr)
{
  flagThreadCtx_t*      threadCtx   = (flagThreadCtx_t*)ctxPtr;
  int64_t               startBucket = threadCtx->startBucket;
  int64_t               endBucket   = threadCtx->endBucket;
  decompHashTableCtx_t* ctx         = threadCtx->ctx;
  threadCtx->errMsg                 = NULL;
  hashrec_t* trail[4096];
  int        trailLen, i, j, probeDist;
  double     squeezeRatio = (double)ctx->cfgHdr.tableSize64ths / 64;  // This floating point value is exact
  uint32_t   wrapBuckets  = (uint32_t)(MAX_WRAP_BYTES * squeezeRatio) >> HASH_BUCKET_BYTES_LOG2;
  uint32_t   threadMask, chainMask, chainList;
  uint8_t*   chainListBytes = (uint8_t*)&chainList;
  int        srcThreadId, thread, listVal, maskVal, chainListLen, listDiffs, listTemp, opcode;
#define SOURCE_THREAD_MASK (((1 << THREAD_ID_BITS) - 1) & ~BUCKET_THREAD_MASK)

  // Loop through starting buckets
  int64_t    bucketIndex, initBucket;
  hashrec_t *bucket, *chainSlot, *slot;
  for (initBucket = startBucket; initBucket < endBucket; initBucket++) {
    // Chase the probing or chaining trail
    bucketIndex = initBucket;
    bucket      = &ctx->hashTable[bucketIndex << HASH_RECORDS_PER_BUCKET_LOG2];
    chainSlot   = &bucket[HASH_RECORDS_PER_BUCKET - 1];
    trailLen    = 0;
    PRINTIF(DEBUG_FLAG, ("----------------------------------------\n"));
    PRINTIF(DEBUG_FLAG, ("Setting flags starting at bucketIndex=0x%llX\n", bucketIndex));
    // Chaining if chain-begin record in the target bucket's last slot
    if (HASH_OPC_IS_BEG(chainSlot->general.opcode)) {
      PRINTIF(DEBUG_FLAG, ("Chain record in slot %d: 0x%016llX\n", chainSlot - bucket, chainSlot->qword));
      // Copy whole first bucket's slot pointers into our trail
      for (i = 0; i < HASH_RECORDS_PER_BUCKET; i++) {
        PRINTIF(DEBUG_FLAG, ("Saving trail record from slot %d: 0x%016llX\n", i, bucket[i].qword));
        trail[trailLen++] = &bucket[i];
      }
      // Follow chain pointers until terminator flag in a chain record
      while (!chainSlot->chain.chain_pad) {
        // Follow chain pointer
        bucketIndex = (bucketIndex & ((uint64_t)0 - CHAIN_BLOCK_BUCKETS)) | chainSlot->chain.chain_ptr;
        bucket      = &ctx->hashTable[bucketIndex << HASH_RECORDS_PER_BUCKET_LOG2];
        PRINTIF(
            DEBUG_FLAG,
            ("Chained to bucketIndex=0x%llX, byteAddr=0x%llX\n",
             bucketIndex,
             bucketIndex << HASH_BUCKET_BYTES_LOG2));
        // Scan next bucket for chain continuation record
        for (i = 0; i < HASH_RECORDS_PER_BUCKET && !HASH_OPC_IS_CON(bucket[i].general.opcode); i++)
          ;
        if (i == HASH_RECORDS_PER_BUCKET) {
          threadCtx->errMsg = "Decomressed hash table bucket missing chain continuation record";
          return NULL;
        }
        chainSlot = &bucket[i++];
        PRINTIF(DEBUG_FLAG, ("Chain record in slot %d: 0x%016llX\n", chainSlot - bucket, chainSlot->qword));
        // Copy subsequent slot pointers into our trail
        for (; i < HASH_RECORDS_PER_BUCKET && bucket[i].general.opcode != HASH_OPC_EMPTY; i++) {
          PRINTIF(DEBUG_FLAG, ("Saving trail record from slot %d: 0x%016llX\n", i, bucket[i].qword));
          trail[trailLen++] = &bucket[i];
        }
        // Chain slot follows them in the trail
        PRINTIF(
            DEBUG_FLAG,
            ("Saving trail record from slot %d: 0x%016llX\n", chainSlot - bucket, chainSlot->qword));
        trail[trailLen++] = chainSlot;
        PRINTIF(DEBUG_FLAG, ("Copied further slots into trail\n"));
        if (trailLen > 4000) {
          threadCtx->errMsg = "Decompressed hash table has improperly long chain";
          return NULL;
        }
      }
    }
    // Probing
    else {
      srcThreadId = (bucketIndex << (THREAD_ID_BITS - HASH_EXTRA_BITS)) & SOURCE_THREAD_MASK;
      probeDist   = 0;
      PRINTIF(DEBUG_FLAG, ("Probing, srcThreadId = 0x%X\n", srcThreadId));
      while (1) {
        PRINTIF(
            DEBUG_FLAG,
            ("Probing step %d: bucketIndex=0x%llX, byteAddr=0x%llX\n",
             probeDist,
             bucketIndex,
             bucketIndex << HASH_BUCKET_BYTES_LOG2));
        // Loop through bucket slots
        for (i = 0; i < HASH_RECORDS_PER_BUCKET; i++) {
          // Quit probing at empty or chain continuation record
          opcode = bucket[i].general.opcode;
          if (opcode == HASH_OPC_EMPTY || HASH_OPC_IS_CON(opcode)) {
            PRINTIF(DEBUG_FLAG, ("Stopping at slot %d record: 0x%016llX\n", i, bucket[i].qword));
            break;
          }
          // Copy non-chain records with threadIDs matching the initial bucket
          if ((bucket[i].matchable.thread_id & SOURCE_THREAD_MASK) == srcThreadId &&
              !HASH_OPC_IS_CHAIN(opcode)) {
            PRINTIF(DEBUG_FLAG, ("Saving trail record from slot %d: 0x%016llX\n", i, bucket[i].qword));
            trail[trailLen++] = &bucket[i];
          }
        }
        if (i < HASH_RECORDS_PER_BUCKET || ++probeDist >= MAX_PROBES) break;
        // Next bucket, wrapping if necessary
        if (++bucketIndex % wrapBuckets == 0) bucketIndex -= wrapBuckets;
        bucket = &ctx->hashTable[bucketIndex << HASH_RECORDS_PER_BUCKET_LOG2];
      }
    }
    PRINTIF(DEBUG_FLAG, ("Trail complete, %d records\n", trailLen));
    chainMask    = 0;
    threadMask   = 0;
    chainList    = 0;
    chainListLen = 0;
    // Walk the trail backwards
    for (i = trailLen - 1; i >= 0; i--) {
      slot = trail[i];
      PRINTIF(DEBUG_FLAG, ("trail[%d] = 0x%016llX\n", i, slot->qword));
      // Special processing for chain records
      if (HASH_OPC_IS_CHAIN(slot->general.opcode)) {
        // Clear the terminal flag we stashed
        slot->chain.chain_pad = 0;
        // Fill in filter list if 1-4 listed so far, else filter mask
        if (chainListLen > 0 && chainListLen <= 4) {
          slot->chain.opcode |= HASH_OPC_CHAIN_LIST_FLAG;
          slot->chain.filter = chainList;
        } else {
          slot->chain.filter = chainMask;
        }
        continue;
      }
      thread  = slot->matchable.thread_id & BUCKET_THREAD_MASK;
      maskVal = slot->matchable.hash_bits & CHAIN_MASK_HASH_MASK;
      listVal = slot->matchable.hash_bits & CHAIN_LIST_HASH_MASK;
      PRINTIF(
          DEBUG_FLAG,
          ("  thread=%d, maskVal=%d, listVal=%d, threadMask=0x%02X, chainMask=0x%08X, chainListLen=%d, chainList=%d,%d,%d,%d\n",
           thread,
           maskVal,
           listVal,
           threadMask,
           chainMask,
           chainListLen,
           chainListBytes[0],
           chainListBytes[1],
           chainListBytes[2],
           chainListBytes[3]));
      // Flag the last instance of each thread ID (first seen walking backwards)
      if (!(threadMask & (1 << thread))) {
        threadMask |= 1 << thread;
        slot->matchable.lf = 1;
        PRINTIF(DEBUG_FLAG, ("  Setting last flag => 0x%016llX\n", slot->qword));
      }
      // Flag filter mask bits required by records ahead (seen so far)
      chainMask |= ((uint32_t)1 << maskVal);
      // List up to 4 filter list values required by records ahead (seen so far)
      if (chainListLen <= 4) {
        for (listDiffs = 0; (listDiffs < chainListLen) && (listVal != chainListBytes[listDiffs]);
             listDiffs++) {
          // Keep list sorted
          if (listVal < chainListBytes[listDiffs]) {
            listTemp                  = listVal;
            listVal                   = chainListBytes[listDiffs];
            chainListBytes[listDiffs] = listTemp;
            continue;
          }
        }
        if (listDiffs == chainListLen) {
          for (j = chainListLen++; j < 4; j++) chainListBytes[j] = listVal;
        }
      }
    }
  }
  return NULL;
}

// Routine to annotate last flags and chain filters into the decompressed hash table
char* decompHashTableSetFlags(decompHashTableCtx_t* ctx)
{
  char*           errMsg        = NULL;
  int             numThreads    = ctx->threads, i;
  int64_t         numBuckets    = ctx->hashRecs >> HASH_RECORDS_PER_BUCKET_LOG2;
  int64_t         threadBuckets = numBuckets / numThreads;
  pthread_t       threads[numThreads];
  flagThreadCtx_t threadCtx[numThreads];

  // Start worker threads
  for (i = 0; i < numThreads; i++) {
    PRINTIF(DEBUG_THREAD, ("Spawning worker thread %d\n", i));
    threadCtx[i].threadId    = i;
    threadCtx[i].startBucket = threadBuckets * i;
    threadCtx[i].endBucket   = i == numThreads - 1 ? numBuckets : threadBuckets * (i + 1);
    threadCtx[i].ctx         = ctx;
    if (pthread_create(&threads[i], NULL, flagsThread, &threadCtx[i])) break;
  }
  // Cancel successful threads if any failed
  if (i < numThreads) {
    for (; i >= 0; i--) pthread_cancel(threads[i]);
    return "Decompressor unable to create flags worker thread";
  }
  // Wait for threads to return
  for (i = 0; i < numThreads; i++) {
    PRINTIF(DEBUG_THREAD, ("Joining worker thread %d\n", i));
    pthread_join(threads[i], NULL);
    if (threadCtx[i].errMsg && !errMsg) errMsg = threadCtx[i].errMsg;
  }
  PRINTIF(DEBUG_THREAD, ("All threads joined\n"));
  return NULL;
}

uint32_t decompHashTableDigest(decompHashTableCtx_t* ctx)
{
  uint64_t* p      = (uint64_t*)ctx->hashTable;
  uint64_t* e      = (uint64_t*)ctx->hashTable + ctx->hashRecs;
  uint64_t  digest = 0;
  for (; p != e; p++) {
#if defined(LOCAL_BUILD) && defined(__x86_64__)
    __asm__ __volatile__(
        "crc32q\t"
        "(%1), %0"
        : "=r"(digest)
        : "r"(p), "0"(digest));
#elif !defined(LOCAL_BUILD)
    digest = crc32c_hw(digest, (const unsigned char*)p, 8);
#endif
  }
  return (uint32_t)digest;
}

uint32_t decompExtendTableDigest(decompHashTableCtx_t* ctx)
{
  uint64_t* p      = (uint64_t*)ctx->extendTable;
  uint64_t* e      = (uint64_t*)ctx->extendTable + ctx->extendRecs;
  uint64_t  digest = 0;
  for (; p != e; p++) {
#if defined(LOCAL_BUILD) && defined(__x86_64__)
    __asm__ __volatile__(
        "crc32q\t"
        "(%1), %0"
        : "=r"(digest)
        : "r"(p), "0"(digest));
#elif !defined(LOCAL_BUILD)
    digest = crc32c_hw(digest, (const unsigned char*)p, 8);
#endif
  }
  return (uint32_t)digest;
}

#if DEBUG_TIME
#include <sys/time.h>
double dblTime()
{
  struct timeval currentTime;
  gettimeofday(&currentTime, NULL);
  return currentTime.tv_sec + currentTime.tv_usec / 1000000.0;
}
#else
#define dblTime() 0
#endif

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
    uint32_t* extTabDigest)
{
  char*                errMsg = NULL;
  double               t0, t1;
  decompHashTableCtx_t ctx;

  PRINTIF(DEBUG_TIME, ("decompHashTableCtxInit...\n"));
  t0 = dblTime();
  if ((errMsg = decompHashTableCtxInit(
           &ctx,
           threads,
           compBuf,
           compLen,
           refBuf,
           refLen,
           *decompBuf,
           *decompLen,
           *extendTableBuf,
           *extendTableLen)))
    goto decompHashTableErr;
  t1 = dblTime();
  PRINTIF(DEBUG_TIME, ("  %.3f seconds\n", t1 - t0));
  t0 = t1;
  PRINTIF(DEBUG_TIME, ("decompHashTableHeader...\n"));
  if ((errMsg = decompHashTableHeader(&ctx))) goto decompHashTableErr;
  t1 = dblTime();
  PRINTIF(DEBUG_TIME, ("  %.3f seconds\n", t1 - t0));
  t0 = t1;
  PRINTIF(DEBUG_TIME, ("decompHashTableLiterals...\n"));
  if ((errMsg = decompHashTableLiterals(&ctx))) goto decompHashTableErr;
  t1 = dblTime();
  PRINTIF(DEBUG_TIME, ("  %.3f seconds\n", t1 - t0));
  t0 = t1;
  PRINTIF(DEBUG_TIME, ("decompHashTableExtIndex...\n"));
  if ((errMsg = decompHashTableExtIndex(&ctx))) goto decompHashTableErr;
  t1 = dblTime();
  PRINTIF(DEBUG_TIME, ("  %.3f seconds\n", t1 - t0));
  t0 = t1;
  PRINTIF(DEBUG_TIME, ("decompHashTableAutoHits...\n"));
  if ((errMsg = decompHashTableAutoHits(&ctx))) goto decompHashTableErr;
  if (ctx.bitPos != ctx.bufBits) {
    errMsg = "Unexpected additional bytes after end of compressed hash table";
    goto decompHashTableErr;
  }
  t1 = dblTime();
  PRINTIF(DEBUG_TIME, ("  %.3f seconds\n", t1 - t0));
  t0 = t1;
  PRINTIF(DEBUG_TIME, ("decompHashTableSetFlags...\n"));
  if ((errMsg = decompHashTableSetFlags(&ctx))) goto decompHashTableErr;
  t1 = dblTime();
  PRINTIF(DEBUG_TIME, ("  %.3f seconds\n", t1 - t0));
  t0 = t1;
  if (hashDigest) {
    PRINTIF(DEBUG_TIME, ("decompHashTableDigest...\n"));
    *hashDigest = decompHashTableDigest(&ctx);
    t1          = dblTime();
    PRINTIF(DEBUG_TIME, ("  %.3f seconds\n", t1 - t0));
    t0 = t1;
  }
  if (extTabDigest) {
    PRINTIF(DEBUG_TIME, ("decompExtendTableDigest...\n"));
    *extTabDigest = decompExtendTableDigest(&ctx);
    t1            = dblTime();
    PRINTIF(DEBUG_TIME, ("  %.3f seconds\n", t1 - t0));
    t0 = t1;
  }
  free(ctx.extIndexRecs);
  ctx.extIndexRecs = NULL;
  free(ctx.refBuf2Bit);
  ctx.refBuf2Bit = NULL;
  PRINTIF(DEBUG_TIME, ("finished decompress\n"));
  if (!*decompBuf) {
    *decompBuf = (uint8_t*)ctx.hashTable;
    *decompLen = ctx.hashRecs * HASH_RECORD_BYTES;
  }
  if (!*extendTableBuf) {
    *extendTableBuf = (uint8_t*)ctx.extendTable;
    *extendTableLen = ctx.extendRecs * sizeof(extend_hit_t);
  }
  free(ctx.priCrcInit);
  free(ctx.secCrcInit);
  return NULL;
decompHashTableErr:
  //  free(ctx.extIndexRecs);
  //  free(ctx.refBuf2Bit);
  //  if (!*decompBuf) free(ctx.hashTable);
  //  if (!*extendTableBuf) free(ctx.extendTable);
  return errMsg;
}

void slurpFile(const char* name, uint8_t** buffer, uint64_t* bufLen)
{
  printf("Slurping file %s...\n", name);
  fflush(stdout);
  struct stat s;
  if (stat(name, &s) < 0) {
    fprintf(stderr, "Cannot stat file '%s'\n", name);
    exit(1);
  }
  uint64_t len = s.st_size;
  uint8_t* buf = malloc(len);
  if (!buf) {
    fprintf(stderr, "Cannot allocate %lld bytes to slurp '%s'\n", len, name);
    exit(1);
  }
  FILE* inp = fopen(name, "rb");
  if (!inp) {
    fprintf(stderr, "Cannot open '%s'\n", name);
    exit(1);
  }
  if (fread(buf, 1, len, inp) != len) {
    fprintf(stderr, "Cannot read %lld bytes from '%s'\n", len, name);
    exit(1);
  }
  fclose(inp);
  printf("  Done, %lld bytes\n", len);
  fflush(stdout);
  *buffer = buf;
  *bufLen = len;
}

void dumpFile(const char* name, uint8_t* buffer, uint64_t bufLen)
{
  printf("Dumping file %s...\n", name);
  fflush(stdout);
  FILE* out = fopen(name, "wb");
  if (!out) {
    fprintf(stderr, "Cannot open '%s'\n", name);
    exit(1);
  }
  if (fwrite(buffer, 1, bufLen, out) != bufLen) {
    fprintf(stderr, "Cannot write %lld bytes to '%s'\n", bufLen, name);
    exit(1);
  }
  printf("  Done, %lld bytes\n", bufLen);
  fflush(stdout);
  fclose(out);
}

// Called by build_hash_table and dragen
int decompAndWriteHashTable(
    const char* refPath,
    const char* hashCmpPath,
    const char* hashBinPath,
    const char* extTabPath,
    const int   numThreads)
{
  uint64_t refBinLen = 0, hashCmpLen = 0, decompLen = 0, extendTableLen = 0;
  uint8_t *refBin = NULL, *hashCmp = NULL, *decompBuf = NULL, *extendTableBuf = NULL;
  slurpFile(refPath, &refBin, &refBinLen);
  slurpFile(hashCmpPath, &hashCmp, &hashCmpLen);
  printf("Decompressing hash table...\n");
  fflush(stdout);

  uint32_t hashDigest = 0, extTabDigest = 0;
  char*    err = decompHashTable(
      numThreads,
      hashCmp,
      hashCmpLen,
      refBin,
      refBinLen,
      &decompBuf,
      &decompLen,
      &extendTableBuf,
      &extendTableLen,
      &hashDigest,
      &extTabDigest);

  if (err) {
    fprintf(stderr, "ERROR: %s\n", err);
    exit(1);
  }
  printf("  hash_table.bin digest:   0x%08lX\n", hashDigest);
  printf("  extend_table.bin digest: 0x%08lX\n", extTabDigest);
  fflush(stdout);
  if (hashBinPath) dumpFile(hashBinPath, decompBuf, decompLen);
  if (extTabPath) dumpFile(extTabPath, extendTableBuf, extendTableLen);
  free(decompBuf);
  free(extendTableBuf);
  return 1;
}
