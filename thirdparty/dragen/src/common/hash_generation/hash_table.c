// Copyright 2013-2018 Edico Genome Corporation. All rights reserved.
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

#include "hash_table.h"
#include <inttypes.h>
#include <math.h>
#include <pthread.h>
#include <semaphore.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include "crc_hash.h"
#include "hash_table_compress.h"
#ifndef LOCAL_BUILD
#include "crc32_hw.h"
#include "watchdog.h"
#ifdef _TARGET_PPC_
#include "crc32_powerpc.h"
#endif
/* ============ end of #include directives ============ */
#endif

// The following filter out these specific compiler warnings, which are numerous in
// this file
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wformat"
#pragma GCC diagnostic ignored "-Wmissing-braces"
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#pragma GCC diagnostic ignored "-Wparentheses"
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"

#if defined(_TARGET_PPC_)
extern unsigned int crc32_vpmsum(unsigned int crc, unsigned char* p, unsigned long len);
#endif

//#define DEBUG
#define DEBUG_EXTEND 0
#define DEBUG_CHAIN 0
#define DEBUG_PACK 0
#define PRINTIF(c, x) \
  if (c) {            \
    printf x;         \
    fflush(stdout);   \
  }

// This array contains values 0-15 for upper/lowercase IUB nucleotide codes, and -1 for other characters.
// E.g. ENCODE_BASE['a'] = 1, ENCODE_BASE['R'] = 5, ENCODE_BASE[' '] = -1.
// Tilde ('~') maps to BASE_PAD = 0.
// clang-format off
const char ENCODE_BASE[256] = {
  -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
  -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
  -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
  -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
  -1, 1,14, 2,13,-1,-1, 4,11,-1,-1,12,-1, 3,15,-1,
  -1,-1, 5, 6, 8,-1, 7, 9,-1,10,-1,-1,-1,-1,-1,-1,
  -1, 1,14, 2,13,-1,-1, 4,11,-1,-1,12,-1, 3,15,-1,
  -1,-1, 5, 6, 8,-1, 7, 9,-1,10,-1,-1,-1,-1, 0,-1,
  -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
  -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
  -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
  -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
  -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
  -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
  -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
  -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
// clang-format on

// Convert 4-bit code to base value 0-3 for ACGT or 0 otherwise
const uint8_t baseCodeBase[16] = {0, 0, 1, 0, 2, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0};
// Number of '1' bits in 4-bit codes 0-15
const uint8_t baseCodeNumBases[16] = {0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4};
// Number of '1' bits in 4-bit codes 0-15
const uint8_t baseCodeNumBasesMin1[16] = {1, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4};
// List of positions of '1' bits in 4-bit codes 0-15
const uint8_t baseCodeBaseList[16][4] = {
    {0},          // 0 = 0000
    {0},          // 1 = 0001
    {1},          // 2 = 0010
    {0, 1},       // 3 = 0011
    {2},          // 4 = 0100
    {0, 2},       // 5 = 0101
    {1, 2},       // 6 = 0110
    {0, 1, 2},    // 7 = 0111
    {3},          // 8 = 1000
    {0, 3},       // 9 = 1001
    {1, 3},       // A = 1010
    {0, 1, 3},    // B = 1011
    {2, 3},       // C = 1100
    {0, 2, 3},    // D = 1101
    {1, 2, 3},    // E = 1110
    {0, 1, 2, 3}  // F = 1111
};

void printHistogram(FILE* file, uint64_t* hist, int bins, int indent, int lastPlus)
{
  char     binStr[20], valStr[20], pctStr[20], binLine[120] = "", valLine[120] = "", pctLine[120] = "";
  int      i, x, len = 0, first = 1, max, last = 0;
  uint64_t total = 0;
  double   pct;
  for (i = 0; i < bins; i++) {
    total += hist[i];
    if (hist[i]) last = i;
  }
  for (i = 0; i < bins; i++) {
    if (!binLine[0]) {
      for (x = 0; x < indent; x++) {
        strcat(valLine, " ");
        strcat(binLine, " ");
        strcat(pctLine, " ");
        len++;
      }
      if (total == 0) {
        fprintf(file, "%s-\n", binLine);
        fprintf(file, "%s-\n", valLine);
        fprintf(file, "%s-\n", pctLine);
        fflush(file);
        return;
      }
    }
    if (hist[i]) {
      pct = (double)hist[i] / total * 100;
      if (pct >= 10)
        sprintf(pctStr, "%.0f%%", pct);
      else if (pct > 1)
        sprintf(pctStr, "%.1f%%", pct);
      else if (pct > .005)
        sprintf(pctStr, "%.2f%%", pct);
      else
        sprintf(pctStr, "-");
      sprintf(binStr, "%u", i);
      sprintf(valStr, "%llu", hist[i]);
      if (lastPlus && i == bins - 1) strcat(binStr, "+");
      max = strlen(valStr);
      max = (max < strlen(pctStr) ? strlen(pctStr) : max);
      max = (max < strlen(binStr) ? strlen(binStr) : max);
      while (strlen(binStr) < max) strcat(binStr, " ");
      while (strlen(valStr) < max) strcat(valStr, " ");
      while (strlen(pctStr) < max) strcat(pctStr, " ");
      strcat(valLine, " ");
      strcat(binLine, " ");
      strcat(pctLine, " ");
      strcat(valLine, valStr);
      strcat(binLine, binStr);
      strcat(pctLine, pctStr);
      len += strlen(binStr) + 1;
    }
    if ((len > 92 && i != last) || i == bins - 1) {
      if (!first) {
        for (x = 0; x < indent; x++) fprintf(file, " ");
        for (; x < 100; x++) fprintf(file, "-");
        fprintf(file, "\n");
      }
      fprintf(file, "%s\n", binLine);
      fprintf(file, "%s\n", valLine);
      fprintf(file, "%s\n", pctLine);
      binLine[0] = 0;
      valLine[0] = 0;
      pctLine[0] = 0;
      first      = 0;
      len        = 0;
    }
  }
  fflush(file);
}

char* bytesReadable(uint64_t bytes)
{
  static char s[256][64], units[][4] = {"B", "KB", "MB", "GB", "TB", "PB", "EB"};
  static int  index = 0;
  int         unit  = 0;
  index             = (index + 1) % 256;
  if (bytes == 0) {
    sprintf(s[index], "0 B");
    return s[index];
  }
  while ((bytes & 0x3FF) == 0) {
    bytes >>= 10;
    unit++;
  }
  if (bytes < 1024) {
    sprintf(s[index], "%llu %s", bytes, units[unit]);
    return s[index];
  }
  double f = bytes;
  while (f >= 1024) {
    f /= 1024;
    unit++;
  }
  f = (double)((uint64_t)(f * 100000)) / 100000;
  sprintf(s[index], "%g %s", f, units[unit]);
  return s[index];
}

// clang-format off
const uint8_t rcTable[256] = {
  0xFF,0xBF,0x7F,0x3F,0xEF,0xAF,0x6F,0x2F,0xDF,0x9F,0x5F,0x1F,0xCF,0x8F,0x4F,0x0F,
  0xFB,0xBB,0x7B,0x3B,0xEB,0xAB,0x6B,0x2B,0xDB,0x9B,0x5B,0x1B,0xCB,0x8B,0x4B,0x0B,
  0xF7,0xB7,0x77,0x37,0xE7,0xA7,0x67,0x27,0xD7,0x97,0x57,0x17,0xC7,0x87,0x47,0x07,
  0xF3,0xB3,0x73,0x33,0xE3,0xA3,0x63,0x23,0xD3,0x93,0x53,0x13,0xC3,0x83,0x43,0x03,
  0xFE,0xBE,0x7E,0x3E,0xEE,0xAE,0x6E,0x2E,0xDE,0x9E,0x5E,0x1E,0xCE,0x8E,0x4E,0x0E,
  0xFA,0xBA,0x7A,0x3A,0xEA,0xAA,0x6A,0x2A,0xDA,0x9A,0x5A,0x1A,0xCA,0x8A,0x4A,0x0A,
  0xF6,0xB6,0x76,0x36,0xE6,0xA6,0x66,0x26,0xD6,0x96,0x56,0x16,0xC6,0x86,0x46,0x06,
  0xF2,0xB2,0x72,0x32,0xE2,0xA2,0x62,0x22,0xD2,0x92,0x52,0x12,0xC2,0x82,0x42,0x02,
  0xFD,0xBD,0x7D,0x3D,0xED,0xAD,0x6D,0x2D,0xDD,0x9D,0x5D,0x1D,0xCD,0x8D,0x4D,0x0D,
  0xF9,0xB9,0x79,0x39,0xE9,0xA9,0x69,0x29,0xD9,0x99,0x59,0x19,0xC9,0x89,0x49,0x09,
  0xF5,0xB5,0x75,0x35,0xE5,0xA5,0x65,0x25,0xD5,0x95,0x55,0x15,0xC5,0x85,0x45,0x05,
  0xF1,0xB1,0x71,0x31,0xE1,0xA1,0x61,0x21,0xD1,0x91,0x51,0x11,0xC1,0x81,0x41,0x01,
  0xFC,0xBC,0x7C,0x3C,0xEC,0xAC,0x6C,0x2C,0xDC,0x9C,0x5C,0x1C,0xCC,0x8C,0x4C,0x0C,
  0xF8,0xB8,0x78,0x38,0xE8,0xA8,0x68,0x28,0xD8,0x98,0x58,0x18,0xC8,0x88,0x48,0x08,
  0xF4,0xB4,0x74,0x34,0xE4,0xA4,0x64,0x24,0xD4,0x94,0x54,0x14,0xC4,0x84,0x44,0x04,
  0xF0,0xB0,0x70,0x30,0xE0,0xA0,0x60,0x20,0xD0,0x90,0x50,0x10,0xC0,0x80,0x40,0x00};
// clang-format on

void copyBases(void* dst, void* src, int len, int srcOffset)
{
  srcOffset &= 3;  // Source offset may be a large absolute position, but we just want the base within a byte
  int      topByte = (len - 1) >> 2, i, shift = srcOffset << 1, rshift = 8 - shift;
  uint8_t *d = dst, *s = src;
  for (i = 0; i < topByte; i++) d[i] = (s[i] >> shift) | (s[i + 1] << rshift);
  d[topByte] = s[topByte] >> shift;
  // Get more bits from the next source byte if still missing some (conditional to avoid memory access error)
  if ((len - (topByte << 2)) > (4 - srcOffset)) d[topByte] |= s[topByte + 1] << rshift;
  // Mask unused bits from destination
  d[topByte] &= (0xFF >> (((4 - (len & 3)) & 3) << 1));
}

void revComp(void* dst, void* src, int len)
{
  int      topByte = (len - 1) >> 2, i, shift = ((4 - (len & 3)) & 3) << 1;
  uint8_t  buf[topByte + 1];
  uint8_t *d = dst, *s = src;
  // Copy source aside first if working in-place
  if (dst == src) {
    memcpy(buf, src, topByte + 1);
    s = buf;
  }
  for (i = 0; i <= topByte; i++) d[i] = rcTable[s[topByte - i]];
  if (shift) {
    if (len <= 32) {
      *((uint64_t*)d) >>= shift;
    } else {
      for (i = 0; i < topByte; i++) d[i] = (d[i] >> shift) | (d[i + 1] << (8 - shift));
      d[topByte] = (d[topByte] >> shift);
    }
  }
}

uint64_t qwordExtractBits(uint64_t qword, int start, int len)
{
  uint64_t val;
  if (start < 0)
    val = qword << -start;
  else
    val = (qword >> start);
  return val & (((uint64_t)1 << len) - 1);
}

uint64_t hashRecExtractHashBits(hashrec_t rec, int start, int len)
{
  uint64_t hash = rec.matchable.hash_bits;
  return qwordExtractBits(hash, start, len);
}

// Hash record sorting comparison function - first by match bits, then by refPosition, then by RC flag (for
// palindromes)
int hashRecCompareHash(const void* a, const void* b)
{
  uint32_t aHashBits = ((hashrec_t*)a)->match_bits.match_bits;
  uint32_t bHashBits = ((hashrec_t*)b)->match_bits.match_bits;
  uint32_t aRefPos   = ((hashrec_t*)a)->hit.pos;
  uint32_t bRefPos   = ((hashrec_t*)b)->hit.pos;
  int      aRc       = ((hashrec_t*)a)->hit.rc;
  int      bRc       = ((hashrec_t*)b)->hit.rc;
  // clang-format off
  return  aHashBits < bHashBits ? -1
        : aHashBits > bHashBits ? 1
        : aRefPos < bRefPos ? -1
        : aRefPos > bRefPos ? 1
        : aRc < bRc ? -1
        : aRc > bRc ? 1
        : 0;
  // clang-format on
}

// Hash record sorting comparison function.  Custom sort keys should be placed
// as 3 additional qwords following each record.  The first sort key is the
// last qword, and the last sort key is the hash record qword.
int hashRecCompareCustom(const void* a, const void* b)
{
  uint64_t* A = (uint64_t*)a;
  uint64_t* B = (uint64_t*)b;
  // clang-format off
  return  A[3] < B[3] ? -1
        : A[3] > B[3] ? 1
        : A[2] < B[2] ? -1
        : A[2] > B[2] ? 1
        : A[1] < B[1] ? -1
        : A[1] > B[1] ? 1
        : A[0] < B[0] ? -1
        : A[0] > B[0] ? 1
        : 0;
  // clang-format on
}

// Simple qword comparison function
int qwordCompare(const void* a, const void* b)
{
  // clang-format off
  return  *(uint64_t*)a < *(uint64_t*)b ? -1
        : *(uint64_t*)a > *(uint64_t*)b ? 1
        : 0;
  // clang-format on
}

// Sortable record for seed extension analysis
typedef struct {
  uint8_t   extension[MAX_NET_SEED_EXTENSION];
  hashrec_t rec;
  int       masked;
  int       maxLen;
  int       liftCode;
  uint32_t  liftIdx;
  uint32_t  liftGroup;
} kmerToExtend_t;
#define LIFT_NON_ALT 0xFFFFFFFF
#define LIFT_NONE 0xFFFFFFFE
#define LIFT_UNIQ_ID_INVALID 0xFFFFFFFF

// Center-symmetric-lexicographical comparison of kmers to extend, for sorting
int kmerToExtendCompare(const void* a, const void* b)
{
  const kmerToExtend_t* A = a;
  const kmerToExtend_t* B = b;
  // Extension bases already toggle from each end, for linear comparison
  int x = memcmp(A->extension, B->extension, MAX_NET_SEED_EXTENSION);
  if (x) return x;
  // Make the sort unique by using reference genome position as a further key
  uint32_t aRefPos = A->rec.hit.pos;
  uint32_t bRefPos = B->rec.hit.pos;
  return aRefPos < bRefPos ? -1 : aRefPos > bRefPos ? 1 : 0;
}

const seedPopRec_t seedPopRecNull = {0};

typedef struct extendTreeRec {
  double   cost;
  uint32_t frequency;
  uint16_t liftExtraFreq;
  uint8_t  extendLen;
  uint8_t  anyAlt : 1;
  uint8_t  masked : 1;
} extendTreeRec;

typedef struct extendSummaryRec {
  uint32_t intervalBeg;
  uint32_t intervalEnd;
  uint16_t extraLifts;
  uint8_t  extendLen;
  uint8_t  depth;
} extendSummaryRec;

// clang-format off
const uint8_t thinningPlacementOrder[THINNING_MAX_PERIOD][THINNING_MAX_PERIOD] = {
  {1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
  {1,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
  {1,  2,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
  {1,  3,  2,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
  {1,  4,  2,  3,  5,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
  {1,  3,  5,  2,  4,  6,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
  {1,  4,  6,  2,  5,  3,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0},
  {1,  5,  3,  7,  2,  6,  4,  8,  0,  0,  0,  0,  0,  0,  0,  0},
  {1,  6,  4,  8,  2,  5,  7,  3,  9,  0,  0,  0,  0,  0,  0,  0},
  {1,  6,  8,  4,  9,  2,  5,  7,  3, 10,  0,  0,  0,  0,  0,  0},
  {1,  8,  6,  4, 10,  2,  7,  9,  3,  5, 11,  0,  0,  0,  0,  0},
  {1,  6, 10,  4,  8, 11,  2,  5,  9,  3,  7, 12,  0,  0,  0,  0},
  {1,  6, 10,  4, 12,  7,  2,  9,  5, 11,  3,  8, 13,  0,  0,  0},
  {1, 10,  6, 13,  4,  8, 12,  2,  9,  5, 11,  3,  7, 14,  0,  0},
  {1, 10, 12,  7,  4,  9, 14,  2,  6, 11, 13,  3,  8,  5, 15,  0},
  {1, 12,  6, 14,  4, 10,  8, 15,  2, 11,  5, 13,  3,  9,  7, 16}};
// clang-format on

const char BASE_CHARS[5]          = "ACGT";
const char COMP_CHARS[5]          = "TGCA";
const char MAYBE_COMP_CHARS[2][5] = {"ACGT", "TGCA"};

const char CODE_CHARS[17]               = ".ACMGRSVTWYHKDBN";
const char CODE_COMP_CHARS[17]          = ".TGKCYSBAWRDMHVN";
const char CODE_MAYBE_COMP_CHARS[2][17] = {".ACMGRSVTWYHKDBN", ".TGKCYSBAWRDMHVN"};

// Routine to add a hash record to a bucket, inreasing allocation if necessary
int bucketAddRecord(hashrec_t** bucketPtr, uint32_t* bucketAllocPtr, uint32_t* bucketCountPtr, hashrec_t rec)
{
  // Check allocation in the bucket.  Increase in powers of two minus one (1,3,7,15,...) because
  // 64-bit memory management appears to work in blocks of 16 bytes, but steals 8 bytes for metadata.
  if (*bucketCountPtr == *bucketAllocPtr) {
    // If the allocation is the physical size, the bucket pointer still references physical record space;
    // redirect it to a larger virtual bucket, and copy the contents there
    if (*bucketAllocPtr == HASH_RECORDS_PER_BUCKET) {
      *bucketAllocPtr = HASH_RECORDS_PER_BUCKET * 2 - 1;
      hashrec_t* tp   = malloc((uint64_t)(*bucketAllocPtr) * HASH_RECORD_BYTES);
      if (!tp) return 1;
      memcpy(tp, *bucketPtr, HASH_BUCKET_BYTES);
      *bucketPtr = tp;
    }
    // Otherwise reallocate the virtual bucket larger
    else {
      *bucketAllocPtr = (*bucketAllocPtr + 1) * 2 - 1;
      hashrec_t* tp   = realloc(*bucketPtr, (uint64_t)(*bucketAllocPtr) * HASH_RECORD_BYTES);
      if (!tp) return 1;
      *bucketPtr = tp;
    }
  }
  // Populate the record
  (*bucketPtr)[*bucketCountPtr] = rec;
  (*bucketCountPtr)++;
  return 0;
}

// Routine to free a virtual bucket if it has shrunk back below the physical bucket size
void shrinkBucket(hashrec_t** bucketPtr, uint32_t* bucketAllocPtr, uint32_t bucketCount, hashrec_t* physPtr)
{
  if (*bucketAllocPtr > HASH_RECORDS_PER_BUCKET && bucketCount <= HASH_RECORDS_PER_BUCKET) {
    memcpy(physPtr, *bucketPtr, bucketCount * HASH_RECORD_BYTES);
    free(*bucketPtr);
    *bucketPtr      = physPtr;
    *bucketAllocPtr = HASH_RECORDS_PER_BUCKET;
  }
}

// Bucket hash record insertion for simple use in each routine
#define BUCKET_ADD_REC(n, r)                                              \
  if (bucketAddRecord(&bucket[n], &bucketAlloc[n], &bucketCount[n], r)) { \
    sprintf(errMsg, "Bucket allocation failure");                         \
    goto buildHashTablesError;                                            \
  }
#define MAIN_BUCKET_ADD_REC(c, n, r)                                               \
  if (bucketAddRecord(&bucket[c][n], &bucketAlloc[c][n], &bucketCount[c][n], r)) { \
    sprintf(errMsg, "Bucket allocation failure");                                  \
    goto mainBuildHashTablesError;                                                 \
  }
#define THREAD_BUCKET_ADD_REC(c, n, r) \
  (fail = bucketAddRecord(&bucket[c][n], &bucketAlloc[c][n], &bucketCount[c][n], r))
// Same for the special bucket
#define SPECIAL_ADD_REC(r)                                                                                 \
  if (prevSpecialCount = specialCount, bucketAddRecord(&specialBucket, &specialAlloc, &specialCount, r)) { \
    sprintf(errMsg, "Special bucket allocation failure");                                                  \
    goto buildHashTablesError;                                                                             \
  }
// Bucket shrink
#define BUCKET_SHRINK(n) \
  shrinkBucket(&bucket[n], &bucketAlloc[n], bucketCount[n], (physRecords + (n)*HASH_RECORDS_PER_BUCKET))
#define BUCKET_CLEAR(n) \
  shrinkBucket(&bucket[n], &bucketAlloc[n], (bucketCount[n] = 0), (physRecords + (n)*HASH_RECORDS_PER_BUCKET))
// Fetch one base
#define GET_BASE_BITS(pos) (((refSeq[(pos) >> 2]) >> (((pos)&3) << 1)) & 3)
#define GET_BASE_CHAR(pos) BASE_CHARS[GET_BASE_BITS(pos)]
// Fetch one mask bit
#define GET_BASE_MASK(pos) (((refMask[(pos) >> 3]) >> ((pos)&7)) & 1)
// Fetch one base 4-bit code
#define GET_BASE_CODE(pos) (((refCode[(pos) >> 1]) >> (((pos)&1) << 2)) & 0xF)
// Fetch 32 bases in a uint64_t
#define GET_SEQ_BITS(pos)                                        \
  (((*((uint64_t*)(&refSeq[(pos) >> 2]))) >> (((pos)&3) << 1)) | \
   ((uint64_t)(refSeq[((pos) >> 2) + 8] << (8 - (((pos)&3) << 1)))) << 56)
// Fetch 64 mask bits in a uint64_t
#define GET_SEQ_MASK(pos)                                  \
  (((*((uint64_t*)(&refMask[(pos) >> 3]))) >> ((pos)&7)) | \
   ((uint64_t)(refMask[((pos) >> 3) + 8] << (8 - ((pos)&7)))) << 56)
// Test for a mask flag in the next len <= 64 positions
#define SEQ_MASK_FAIL(pos, len) ((GET_SEQ_MASK(pos) << (64 - (len))) != 0)
// Check if the parent thread says we need to quit
#define CHECK_ABORT            \
  if (*ctx->abort) {           \
    goto buildHashTablesError; \
  }

typedef struct {
  uint64_t kmerFreqHist[64], kmerLogFreqHist[32];
  uint64_t kmerFreqSumSquare, probeLenSum, chainLenSum, extendLenWtdSum, extendIncrWtdSum;
  uint64_t kmerExtendFreqSum, kmerHitFreqSum, extendStepsWtdSum, secKmerHitFreqSum;
  uint64_t countRawKmer, countUniqKmer, countHitRec;
  uint64_t passExtendSeed, passChains, countSeedHits[2];
  uint64_t countExtendRec, countIntervalRec, countExtendKmer, countChainRec, countEmptyRec,
      countKmerExtendIncr;
  uint64_t countChains, countChainBuckets, countThinnedKmers;
  uint64_t extendFreqHist[64], hitFreqHist[2][64];
  uint64_t bucketLevelRawHist[64], bucketLevelEscHist[64], bucketLevelMapHist[64];
  uint64_t extendLengthHist[MAX_EXTENDED_LENGTH + 1], extendIncrHist[64], extendStepsHist[64];
  uint64_t probeLengthHist[64], chainLengthHist[64];
  uint64_t chainOrProbeLengthHist[64], probeDistChainedHist[64];
  uint64_t extendIdBinPctHist[101];
  uint64_t priSeedAltContig, priSeedAltNoLiftover, priSeedAltLiftPresent, priSeedAltLiftAbsent;
  uint64_t rawNoLiftoverCountHist[64], popNoLiftoverCountHist[64], rawLiftoverGroupSizeHist[64],
      popLiftoverGroupSizeHist[64];
  uint64_t rawLiftoverGroupCount, popLiftoverGroupCount, rawLiftoverGroupSizeSum, popLiftoverGroupSizeSum;
  uint64_t liftoverPriHitsMatching, liftoverPriHitsAdded, liftoverDummyPriSeeds;
  uint64_t cyclesOverhead, cyclesBucketOverhead, cyclesExtendPrep, cyclesExtendSort, cyclesLiftover,
      cyclesPriSort;
  uint64_t cyclesExtendFreq, cyclesExtendDynProg, cyclesExtendConstruct, cyclesBucketSort,
      cyclesBucketOrganize;
  uint64_t cyclesBucketChain, cyclesBucketWrite, cyclesBucketCompress, cyclesExtendIntervals;
} buildThreadStats_t;

typedef struct {
  hashTableConfig_t* config;
  uint8_t*           refSeq;
  uint8_t*           refMask;
  uint8_t*           refCode;
  uint8_t*           litFlags;

  seedPopRec_t*      seedPopRecs;
  extend_hit_t*      extendHitRecs;
  extIndexRec_t*     extIndexRecs;
  uint32_t           extIndexBase;
  uint32_t           numExtendHitRecs;
  int                threadId;
  int                chunkNum;
  uint32_t           bucketStart;
  uint32_t*          bucketAlloc;
  uint32_t*          bucketCount;
  hashrec_t**        bucket;
  hashrec_t*         physRecords;
  uint32_t           numBuckets;
  uint32_t           maxExtendIds;
  buildThreadStats_t stats;
  char               errMsg[256];
  int*               abort;
} buildThreadCtx_t;

#if defined(LOCAL_BUILD) && !defined(_TARGET_PPC_)
static inline uint64_t RDTSC()
{
  uint32_t hi, lo;
  __asm__ volatile("rdtsc" : "=a"(lo), "=d"(hi));
  return ((uint64_t)hi << 32) | lo;
}
#define RESET_CYCLES (tsc0 = RDTSC())
#define COUNT_CYCLES(name) (tsc1 = RDTSC(), stats.name += tsc1 - tsc0, tsc0 = tsc1)
#else
#define RESET_CYCLES
#define COUNT_CYCLES(name)
#endif

// Returns NULL on success, or error message
//void __attribute__((optimize("O0")))*buildHashTableThread(void *ctxPtr) {
void* buildHashTableThread(void* ctxPtr)
{
  // Extract arguments
  uint64_t tsc0, tsc1;
  (void)tsc0;
  (void)tsc1;
  RESET_CYCLES;
  buildThreadCtx_t*  ctx          = (buildThreadCtx_t*)ctxPtr;
  hashTableConfig_t* config       = ctx->config;
  hashTableHeader_t* configHeader = config->hdr;
  uint8_t*           refSeq       = ctx->refSeq;
  uint8_t*           refMask      = ctx->refMask;
  uint8_t*           refCode      = ctx->refCode;
  uint8_t*           litFlags     = ctx->litFlags;

  seedPopRec_t*      seedPopRecs  = ctx->seedPopRecs;
  extIndexRec_t*     extIndexRecs = ctx->extIndexRecs;
  uint32_t           extIndexBase = ctx->extIndexBase;
  uint32_t*          bucketAlloc  = ctx->bucketAlloc;
  uint32_t*          bucketCount  = ctx->bucketCount;
  uint32_t           numBuckets   = ctx->numBuckets;
  hashrec_t**        bucket       = ctx->bucket;
  hashrec_t*         physRecords  = ctx->physRecords;
  hashrec_t          rec;
  char*              errMsg = ctx->errMsg;
  buildThreadStats_t stats  = {0};
  ctx->stats                = stats;  // Local copy to avoid indirection on all counter increments
  errMsg[0]                 = 0;
  uint64_t* tempBucket      = NULL;

  uint32_t allocExtendHitRecs = 32;
  ctx->numExtendHitRecs       = 0;
  ctx->extendHitRecs          = malloc(allocExtendHitRecs * sizeof(extend_hit_t));

  uint32_t        freq, numAlt = 0;
  int32_t         i, j, k, n, m, len, VAL, LOG, FRQ;
  uint8_t *       bucketProbeDist = NULL, *bucketOccupancy = NULL;
  kmerToExtend_t* kmerBuf                              = NULL;
  uint8_t*        kmerLiftHitFlag                      = NULL;
  uint32_t*       kmerLiftHitList                      = NULL;
  uint32_t        kmerLiftHitCount                     = 0;
  uint32_t        kmerBufLen                           = 0;
  uint32_t        thinningMaxFreq[THINNING_MAX_PERIOD] = {0xFFFFFFFFul, 0};

  // Compute frequency caps for various seed positions modulo the thinning period
  for (i = 1; i < configHeader->thinningPeriod; i++)
    thinningMaxFreq[i] = (uint32_t)round(
        configHeader->thinningFreqCap * configHeader->thinningPeriod /
        thinningPlacementOrder[configHeader->thinningPeriod][i]);

  // Grab config params
  int      addrBits        = configHeader->tableAddrBits;
  double   refSeedInterval = configHeader->refSeedInterval;
  uint32_t maxSeedFreq     = configHeader->maxSeedFreq;
  uint32_t refAltSeed      = configHeader->refAltSeed;
  int      writeCompFile   = config->writeCompFile;
  // Disable alt-aware seed handling for anchored mapping
  if (configHeader->anchorBinBits) refAltSeed = 0xFFFFFFFF;

  // Derive hash record layout parameters
  int      REC_BYTES    = HASH_RECORD_BYTES;
  double   squeezeRatio = (float)configHeader->tableSize64ths / 64;  // This floating point value is exact
  uint32_t wrapBytes    = MAX_WRAP_BYTES * squeezeRatio;
  uint32_t wrapRecords  = wrapBytes >> HASH_RECORD_BYTES_LOG2;
  uint32_t wrapBuckets  = wrapBytes >> HASH_BUCKET_BYTES_LOG2;
  uint32_t numWraps     = numBuckets / wrapBuckets;

  int   byteAddrStart = HASH_BYTE_ADDR_START;
  int   probeLimit    = MAX_PROBES;
  void* secCrcInit    = crcHash64Init(configHeader->secCrcBits, configHeader->secCrcPoly);
  if (!secCrcInit) return "Failure initializing secondary hash function\n";

  // Dynamic length arrays on the stack
  uint32_t  threadCounts[HASH_THREADS];
  uint8_t   threadFlags[HASH_THREADS];
  uint64_t  bucketProbeDistHist[MAX_PROBES + 1];
  uint64_t  bucketOccupancyHist[HASH_RECORDS_PER_BUCKET + 1];
  uint8_t   recUsed[wrapRecords];
  uint8_t   recsPlaced[wrapBuckets];
  int32_t   wrapOffsets[wrapRecords];
  int8_t    wrapSpecials[wrapRecords];
  hashrec_t physRecs[wrapRecords];

  // Some global stuff
  int seedLen            = configHeader->priSeedBases;
  int seedEven           = !(seedLen & 1);
  int extIdHashBits      = configHeader->secCrcBits - SEC_CRC_BITS_MINUS_EXT_ID_HASH_BITS;
  extIdHashBits          = (extIdHashBits < 0 ? 0 : extIdHashBits);
  uint32_t  numExtIDs    = (uint32_t)1 << extIdHashBits;
  uint32_t* extensionIDs = calloc(numExtIDs, 4);
  if (!extensionIDs || !ctx->extendHitRecs) {
    sprintf(errMsg, "Failed to allocate %d bytes for seed extension IDs\n", numExtIDs);
    goto buildHashTablesError;
  }
  extendTreeRec*    extendTreeRecs[MAX_EXTENDED_LENGTH + 1] = {NULL};
  extendSummaryRec* extendSummaryRecs                       = NULL;
  int               extendTreeRecsAlloc                     = 0;

  // Initialize a 64-bit CRC for use generating pseudo-random numbers
  void* randCrcInit = crcHash64Init(64, CRC_POLYS[64][0]);
  if (!randCrcInit) return "Failure initializing CRC random function\n";

  // Allocate a special "bucket" for HIT records needing special treatment in the compressed format
  uint32_t   specialAlloc = 1023, specialCount = 0, prevSpecialCount = 0;
  hashrec_t* specialBucket = malloc(specialAlloc * HASH_RECORD_BYTES);
  if (!specialBucket) {
    sprintf(errMsg, "Failed to allocate special bucket\n");
    goto buildHashTablesError;
  }

  // Allocate bucket attributes
  bucketProbeDist = malloc(numBuckets);
  bucketOccupancy = malloc(numBuckets);
  if (!bucketProbeDist || !bucketOccupancy) {
    sprintf(errMsg, "Failed to allocate bucket attributes\n");
    goto buildHashTablesError;
  }
#ifndef LOCAL_BUILD
  unsigned int watchdogID = WatchDogRegister("_build_hash_tbl_thread");
#endif
  // Extend or reject high frequency K-mers
  PRINTIF(DEBUG_EXTEND, ("  Extending or rejecting high frequency seeds\n"));
  uint64_t pos;
  uint32_t beg, end = 0, num, priRecs;
  double   cost, minCost;
  int      rc, maxLen, bestLen, nextLen;
  uint32_t maxBucketCount = 0;
  for (n = 0; n < numBuckets; n++)
    if (bucketCount[n] > maxBucketCount) maxBucketCount = bucketCount[n];
  tempBucket = malloc(maxBucketCount * 2 * 8);
  if (!tempBucket) {
    sprintf(errMsg, "Failed to allocate temporary buckets for sort");
    goto buildHashTablesError;
  }
  COUNT_CYCLES(cyclesOverhead);
  for (n = 0; n < numBuckets; n++) {
    CHECK_ABORT;
    // Count un-extended records in the bucket (extended ones will come after)
    // ... also quit if we see an un-extended non-HIT opcode, an EXTEND record
    // we generated earlier while processing previous HIT records in this bucket
    for (priRecs = 0; priRecs < bucketCount[n]; priRecs++)
      if (!HASH_OPC_IS_NORM(bucket[n][priRecs].general.opcode) || bucket[n][priRecs].hit.ex) break;
    // Stats
    stats.bucketLevelRawHist[bucketCount[n] > 63 ? 63 : bucketCount[n]]++;
    stats.countRawKmer += priRecs;
    // Sort the primary portion of the bucket by hash, to group records from identical k-mers
    qsort(bucket[n], priRecs, REC_BYTES, &hashRecCompareHash);
    // Scan for groups of matching hashBits
    for (beg = 0; beg < bucketCount[n]; beg = end) {
      // Quit if we struck an extended record
      if (!HASH_OPC_IS_NORM(bucket[n][beg].general.opcode) || bucket[n][beg].hit.ex) break;
      for (end = beg; end < bucketCount[n]; end++)
        if (!HASH_OPC_IS_NORM(bucket[n][end].general.opcode) ||
            bucket[n][end].match_bits.match_bits != bucket[n][beg].match_bits.match_bits)
          break;
      // Now we have a group of matching hashBits from beg to end-1.  Skip groups of 1.
      num = end - beg;
      // Thin out high frequency seeds.  Convert each record's reference position back to the sequence number
      // by dividing by refSeedInterval, take the remainder modulo thinningPeriod, and look up the maximum
      // frequency for that remainder.  Delete the record if this K-mer group exceeds the maximum.
      if (num > configHeader->thinningFreqCap && configHeader->thinningPeriod > 1) {
        for (i = j = beg; i < end; i++) {
          if (num <= thinningMaxFreq[((uint32_t)bucket[n][i].hit.pos) % configHeader->thinningPeriod])
            bucket[n][j++] = bucket[n][i];
        }
        // Shift further records into the freed space
        stats.countThinnedKmers += end - j;
        stats.countSeedHits[0] -= end - j;
        end = j;
        num = end - beg;
        for (; i < bucketCount[n]; i++) bucket[n][j++] = bucket[n][i];
        bucketCount[n] = j;
        BUCKET_SHRINK(n);
      }
      // Count ALT contig seeds in the group
      numAlt = 0;
      for (i = beg; i < end; i++) numAlt += (bucket[n][i].hit.pos >= refAltSeed);
      // Skip further processing for singleton
      if (num + numAlt <= 1) {
        // Stats
        stats.countSeedHits[0]++;
        stats.kmerFreqHist[1]++;
        stats.hitFreqHist[0][1]++;
        stats.kmerLogFreqHist[0]++;
        stats.kmerFreqSumSquare++;
        stats.kmerHitFreqSum++;
        stats.countUniqKmer++;
        continue;
      }
      // Make sure the k-mer buffer is big enough
      uint32_t kmerBufCurLen = num + numAlt + 1;
      if (kmerBufLen < kmerBufCurLen) {
        free(kmerBuf);
        free(kmerLiftHitFlag);
        free(kmerLiftHitList);
        kmerBufLen      = kmerBufCurLen;
        kmerBuf         = malloc(kmerBufLen * sizeof(kmerToExtend_t));
        kmerLiftHitFlag = malloc(kmerBufLen);
        kmerLiftHitList = malloc(kmerBufLen * 4);
        if (!kmerBuf || !kmerLiftHitFlag || !kmerLiftHitList) {
          sprintf(errMsg, "K-mer buffer allocation failure");
          goto buildHashTablesError;
        }
        memset(kmerLiftHitFlag, 0, kmerBufLen);
      }
      // Make sure the seed extension table is big enough
      uint32_t reserveExtendHitRecs = ctx->numExtendHitRecs + kmerBufCurLen + 32;
      if (allocExtendHitRecs < reserveExtendHitRecs) {
        allocExtendHitRecs = reserveExtendHitRecs + reserveExtendHitRecs / 8 + 4096;
        extend_hit_t* newp = realloc(ctx->extendHitRecs, allocExtendHitRecs * sizeof(extend_hit_t));
        if (!newp) {
          sprintf(errMsg, "K-mer buffer allocation failure");
          goto buildHashTablesError;
        }
        ctx->extendHitRecs = newp;
      }
      extend_hit_t* extendHitRecs = ctx->extendHitRecs;
      PRINTIF(
          DEBUG_EXTEND, ("\nbucketCount[%d]=%d, beg=%d, end=%d, num=%d\n", n, bucketCount[n], beg, end, num));
      maxLen = seedLen + MAX_NET_SEED_EXTENSION;
      maxLen = (maxLen > configHeader->maxSeedBases ? configHeader->maxSeedBases : maxLen);
      maxLen -= ((maxLen - seedLen) & 1);  // Even extension lengths
      COUNT_CYCLES(cyclesBucketOverhead);
      // Each record in the group
      for (i = 0; i < num; i++) {
        j = i + beg;
        // Reference seed position and orientation
        rc             = bucket[n][j].hit.rc;
        pos            = floor(bucket[n][j].hit.pos * refSeedInterval);
        uint64_t left  = pos - 1;
        uint64_t right = pos + seedLen;
        // Sorting record to populate
        kmerToExtend_t* kmer      = &kmerBuf[i];
        uint8_t*        extension = kmer->extension;
        int             maxExt    = maxLen - seedLen;
        int             unmasked  = 4;
        // Alternate extension bases from each end.  Start on the right if forward,
        // or on the left if reverse complemented, complementing bases as we go.
        // Detect when extension would be blocked on either end by a mask bit,
        // and add 4 to each unmasked base, so extension-limited seeds sort
        // before longer extensions they otherwise match.
        if (rc) {
          for (k = 0; k < MAX_NET_SEED_EXTENSION - 1;) {
            if ((GET_BASE_MASK(left) || GET_BASE_MASK(right)) && unmasked) {
              maxExt   = k;
              unmasked = 0;
            }
            extension[k++] = (GET_BASE_BITS(left) ^ 3) + unmasked;
            extension[k++] = (GET_BASE_BITS(right) ^ 3) + unmasked;
            left--;
            right++;
          }
        } else {
          for (k = 0; k < MAX_NET_SEED_EXTENSION - 1;) {
            if ((GET_BASE_MASK(left) || GET_BASE_MASK(right)) && unmasked) {
              maxExt   = k;
              unmasked = 0;
            }
            extension[k++] = GET_BASE_BITS(right) + unmasked;
            extension[k++] = GET_BASE_BITS(left) + unmasked;
            left--;
            right++;
          }
        }
        // Copy the hash record too, and note the maximum extension
        kmer->rec       = bucket[n][j];
        kmer->masked    = SEQ_MASK_FAIL(pos, seedLen);
        kmer->maxLen    = maxExt + seedLen;
        kmer->liftIdx   = LIFT_NON_ALT;
        kmer->liftCode  = LIFT_CODE_NONE;
        kmer->liftGroup = LIFT_UNIQ_ID_INVALID;
      }
      COUNT_CYCLES(cyclesExtendPrep);
      // Sort the 64-base extensions center-symmetric-lexicographically
      qsort(kmerBuf, num, sizeof(kmerToExtend_t), &kmerToExtendCompare);
      COUNT_CYCLES(cyclesExtendSort);
      // Lift ALT contig seeds to the primary assembly
      uint32_t priSeedIdx, numWithPri = num, nextLiftGroup = 0;
      int      dummyPriUsed = 0, dummyPriCount = 0;
      uint32_t liftIdx, groupSize              = 1;

      if (DEBUG_EXTEND >= 2) {
        // Display the sorted pile
        int extLimit = maxLen - seedLen;
        for (i = 0; i < numWithPri; i++) {
          rc  = kmerBuf[i].rec.hit.rc;
          pos = floor(kmerBuf[i].rec.hit.pos * refSeedInterval);
          printf("%6d: ", i);
          if (i < num) {
            int maxExt = kmerBuf[i].maxLen - seedLen;
            int dir    = (rc ? -1 : 1);
            if (rc)
              pos += seedLen - 1 + maxExt / 2;
            else
              pos -= maxExt / 2;
            for (k = extLimit; k > maxExt; k -= 2) {
              putchar('.');
            }
            for (; k > 0; k -= 2) {
              putchar(
                  refCode ? CODE_MAYBE_COMP_CHARS[rc][GET_BASE_CODE(pos)]
                          : MAYBE_COMP_CHARS[rc][GET_BASE_BITS(pos)]);
              pos += dir;
            }
            putchar('|');
            for (k = 0; k < seedLen; k++) {
              putchar(
                  refCode ? CODE_MAYBE_COMP_CHARS[rc][GET_BASE_CODE(pos)]
                          : MAYBE_COMP_CHARS[rc][GET_BASE_BITS(pos)]);
              pos += dir;
            }
            putchar('|');
            for (k = 0; k < maxExt; k += 2) {
              putchar(
                  refCode ? CODE_MAYBE_COMP_CHARS[rc][GET_BASE_CODE(pos)]
                          : MAYBE_COMP_CHARS[rc][GET_BASE_BITS(pos)]);
              pos += dir;
            }
            for (; k < extLimit; k += 2) {
              putchar('.');
            }
          } else {
            printf("Liftover PRI seed");
          }
          printf(" - pos=0x%08llX, rc=%d", (uint64_t)floor(kmerBuf[i].rec.hit.pos * refSeedInterval), rc);
          if (kmerBuf[i].liftIdx != LIFT_NON_ALT)
            printf(", ALT liftover to #%lu (group #%d)", kmerBuf[i].liftIdx, kmerBuf[i].liftGroup);
          printf("\n");
        }
        fflush(stdout);
      }
      // Allocate a 2D array of extend records for use optimizing the seed extension tree.
      if (num > extendTreeRecsAlloc) {
        while (num > extendTreeRecsAlloc)
          extendTreeRecsAlloc = extendTreeRecsAlloc ? extendTreeRecsAlloc * 2 : 1024;
        PRINTIF(DEBUG_EXTEND, ("Allocating extension summary array to %u\n", extendTreeRecsAlloc));
        extendSummaryRec* newp =
            realloc(extendSummaryRecs, (uint64_t)extendTreeRecsAlloc * sizeof(extendSummaryRec));
        if (!newp) {
          sprintf(errMsg, "Extension summary record allocation failure");
          goto buildHashTablesError;
        }
        extendSummaryRecs = newp;
        for (i = 1; i <= maxLen; i++) {
          // Only allocate even extension lengths starting at the primary seed length, and a summary array at
          // zero
          if ((i < seedLen) || ((i - seedLen) & 1)) continue;
          PRINTIF(DEBUG_EXTEND, ("Allocating extension tree arrays to %u: len=%d\n", extendTreeRecsAlloc, i));
          extendTreeRec* newp =
              realloc(extendTreeRecs[i], (uint64_t)extendTreeRecsAlloc * sizeof(extendTreeRec));
          if (!newp) {
            sprintf(errMsg, "Extension record allocation failure");
            goto buildHashTablesError;
          }
          extendTreeRecs[i] = newp;
        }
      }
      COUNT_CYCLES(cyclesBucketOverhead);
      // First measure the seed frequencies as extended to various lengths
      for (len = seedLen; len <= maxLen; len += 2) {
        int clump = 0;
        PRINTIF(DEBUG_EXTEND >= 2, ("Measuring seed frequencies at len=%d\n", len));
        for (i = 0, j = 1; j <= num; j++) {
          kmerToExtend_t* kmerI = &kmerBuf[i];
          kmerToExtend_t* kmerJ = &kmerBuf[j];
          // At the end, or if the next k-mer (at j) differs from the previous interval (starting at i),
          // write the frequency (j-i) into the first record (i).  Also break each k-mer with maximum
          // extension length greater than the current length ('masked' e.g. by N positions) into its own
          // frequency 1 group.
          int masked = (len > kmerI->maxLen);
          if (j == num || masked || memcmp(kmerI->extension, kmerJ->extension, len - seedLen)) {
            freq                             = j - i;
            uint32_t freqWithLift            = freq;
            int      anyAlt                  = 0;
            extendTreeRecs[len][i].frequency = freq;
            extendTreeRecs[len][i].masked    = masked;

            uint32_t liftExtraFreq               = freqWithLift - freq;
            extendTreeRecs[len][i].liftExtraFreq = liftExtraFreq > 0xFFFF ? 0xFFFF : liftExtraFreq;
            extendTreeRecs[len][i].anyAlt        = anyAlt;
            PRINTIF(
                DEBUG_EXTEND >= 2,
                ("  %d: freq=%d, freqWithLift=%d, masked=%d\n", i, freq, freqWithLift, masked));
            // Clump flag for any seed interval exceeding target frequency
            // (at primary seed length, flag when frequency is sufficient to allow extension)
            if (configHeader->maxMultBaseSeeds > 0) {
              clump |= (freqWithLift > configHeader->targetSeedFreq);
              if (len == seedLen && freq < configHeader->minFreqToExtend &&
                  freqWithLift <= configHeader->maxSeedFreq)
                clump = 0;
            } else {
              clump |= (len == seedLen) ? (freqWithLift >= configHeader->minFreqToExtend)
                                        : (freqWithLift > configHeader->targetSeedFreq);
            }
            // Advance past this interval
            i = j;
          }
        }
        // Reduce maxLen to this length and exit the loop if all seeds are at or below the target frequency
        // (never any point to extending longer)
        if (!clump) {
          maxLen = len;
        }
      }
      COUNT_CYCLES(cyclesExtendFreq);
      // Dynamic programming to generate an optimal seed extension tree.  In the column of extendTreeRecs for
      // a given seed extension length, only the first record of each matching-extension interval is used.
      // Each such valid record represents the possibility that the extension tree has a node at its extension
      // length, for its interval's seeds.
      // This record contains:
      //  .frequency:     The interval's length, which among other uses allows skipping to
      //                  the next valid record at the beginning of the following interval
      //  .liftExtraFreq: Number of additional nonmatching liftover hits needed for the interval
      //  .cost:          Minimum cost of seeds in this node's interval, considering all possible
      //                  extension sub-trees starting at this node
      //  .extendLen:     Optimal next extension length beyond this node, corresponding with .cost,
      //                  or zero if no further extension is optimal.  The 'children' of this node
      //                  are all the various sub-intervals at this next extension length.
      //  .masked:        Flag indicating an interval of a single seed that cannot extend to this length
      //  .anyAlt:        Flag indicating interval contains at least one ALT contig hit
      // We fill in columns starting with the maximum extension length.  To fill in a node, all possible next
      // extension lengths are considered, child intervals at that length being identified, and their costs
      // summed, and added to an additional cost term for another extension step.  The cost of populating
      // seeds at the current length without further extension is also calculated, and the lowest cost option
      // is chosen.
      PRINTIF(DEBUG_EXTEND >= 2, ("Dynamic programming for seed extension tree...\n"));
      for (len = maxLen; len >= seedLen; len -= 2) {
        // Ramp the maximum seed frequency from priMaxSeedFreq at the primary seed length (seedLen)
        // to maxSeedFreq at extended seed length maxSeedFreqLen
        int maxFreq = configHeader->maxSeedFreq;
        if (len < configHeader->maxSeedFreqLen && configHeader->priMaxSeedFreq > 0)
          maxFreq = (configHeader->maxSeedFreq - configHeader->priMaxSeedFreq) * (len - seedLen) /
                        (configHeader->maxSeedFreqLen - seedLen) +
                    configHeader->priMaxSeedFreq;
        for (i = 0; i < num; i += freq) {
          freq                  = extendTreeRecs[len][i].frequency;
          uint32_t freqWithLift = freq + extendTreeRecs[len][i].liftExtraFreq;
          // If this interval (always a single seed) is flagged 'masked', meaning its seed cannot be extended
          // to this length, just assign it a large constant cost.  The effect of this will be to select an
          // extension length that minimizes the number of seeds lost to failed extensions, with top priority.
          // But an extension length with one or more failed extensions may still be chosen if it is the only
          // way to get other non-failing seeds under the maximum hit frequency.
          if (extendTreeRecs[len][i].masked) {
            extendTreeRecs[len][i].cost      = 1000000 * configHeader->seedLenCost;
            extendTreeRecs[len][i].extendLen = 0;
            PRINTIF(
                DEBUG_EXTEND >= 2,
                ("extendTreeRecs[%d][%d] (%d):  masked=1, cost=%f\n",
                 len,
                 i,
                 freq,
                 extendTreeRecs[len][i].cost));
            continue;
          }
          // Initialize cost minimization
          minCost = 1e99;
          bestLen = 0;
          // Seed extension always terminates at this length if already at or below the target frequency,
          // or if the frequency doesn't shrink further by the maximum extension length
          int term =
              ((freq <= configHeader->targetSeedFreq) || (freq == extendTreeRecs[maxLen][i].frequency));
          // If frequency here is under the maximum, we have the option to terminate seed extension here;
          // compute this cost based on the current extended seed length and frequency
          if (freqWithLift <= maxFreq || term)
            minCost = freqWithLift *
                      (len * configHeader->seedLenCost +
                       fabs(freqWithLift - configHeader->targetSeedFreq) * configHeader->seedFreqCost);
          if (!term) {
            // Extending cost starts with a factor representing the performance impact of yet another
            // extension step for the number of seeds in this interval; add an extra cost factor for the first
            // step
            double baseCost =
                freqWithLift * (configHeader->extStepCost + (len == seedLen) * configHeader->extensionCost);
            // Calculate the longest next extension length legal from here
            int maxNextLen = len + configHeader->maxExtIncrement;
            maxNextLen     = (maxNextLen > maxLen ? maxLen : maxNextLen);
            // Consider all possible next extension lengths
            for (nextLen = len + 2; nextLen <= maxNextLen; nextLen += 2) {
              cost = baseCost;
              // Add the cost of all child nodes we would get, corresponding to the subintervals at this next
              // extension length
              for (j = 0; j < freq; j += extendTreeRecs[nextLen][i + j].frequency)
                cost += extendTreeRecs[nextLen][i + j].cost;
              // Capture any lower cost found
              if (cost < minCost) {
                minCost = cost;
                bestLen = nextLen;
              }
            }
          }
          // Add a cost component per tree node for the added INTERVAL & EXTEND records
          minCost += (bestLen ? 2 : 1) * configHeader->extRecCost;
          // Record the optimal result
          extendTreeRecs[len][i].cost      = minCost;
          extendTreeRecs[len][i].extendLen = bestLen;
          PRINTIF(
              DEBUG_EXTEND >= 2,
              ("extendTreeRecs[%d][%d] (%d):  cost=%f, extendLen=%d\n", len, i, freq, minCost, bestLen));
        }
      }
      COUNT_CYCLES(cyclesExtendDynProg);

      // The root node is at the primary seed length; leave the primary hash records in this bucket if no
      // extension from there
      bestLen     = extendTreeRecs[seedLen][0].extendLen;
      int maxFreq = configHeader->maxSeedFreq;
      if (configHeader->priMaxSeedFreq > 0 && configHeader->priMaxSeedFreq < maxFreq)
        maxFreq = configHeader->priMaxSeedFreq;
      uint32_t numFinal    = numWithPri + dummyPriUsed;
      int      primaryHits = !bestLen && numFinal <= maxFreq;
      PRINTIF(DEBUG_EXTEND, ("bestLen = %d\n", bestLen));

      // For even seed lengths, detect palindromes, and flag the reverse-complemented copy
      if (writeCompFile && seedEven) {
        // First pass: flag positions with forward hits
        for (i = 0; i < num; i++) {
          rec = kmerBuf[i].rec;
          if (!rec.hit.rc) seedPopRecs[rec.hit.pos].slotOffset = 1;
        }
        // Second pass: detect reverse hits with already-flagged positions
        for (m = i = 0; i < num; i++) {
          rec = kmerBuf[i].rec;
          if (rec.hit.rc && seedPopRecs[rec.hit.pos].slotOffset) {
            // Add to the special bucket for compressed output
            SPECIAL_ADD_REC(rec);
            // Replace the position with a reference to the special bucket
            kmerBuf[i].rec.hit.pos = HASH_REC_POS_SPECIAL | prevSpecialCount;
            // Remember we did this
            m = 1;
          }
        }
        // Third pass: clear the flags
        for (i = 0; i < num; i++) {
          rec = kmerBuf[i].rec;
          if (!rec.hit.rc) seedPopRecs[rec.hit.pos].slotOffset = 0;
        }
        // Fourth pass: copy back into the bucket if we're not using the kmerBuf further
        if (m && primaryHits && !numAlt)
          for (i = 0; i < num; i++) bucket[n][beg + i] = kmerBuf[i].rec;
      }

      // No extension
      if (primaryHits) {

        // Stats
        FRQ = (numFinal > 63 ? 63 : numFinal);
        stats.kmerFreqHist[FRQ] += num;
        for (LOG = 0, VAL = num >> 1; VAL; VAL >>= 1, LOG++)
          ;
        stats.kmerLogFreqHist[LOG] += num;
        stats.kmerFreqSumSquare += (uint64_t)numFinal * num;
        if (numFinal <= maxSeedFreq) {
          stats.hitFreqHist[0][FRQ] += num;
          stats.countSeedHits[0] += num;
          stats.kmerHitFreqSum += (uint64_t)numFinal * num;
        }
        // Done with this [beg,end) group
        COUNT_CYCLES(cyclesBucketOverhead);
        continue;
      }
      // Go ahead and delete all the original HIT records from the bucket.
      // We have them saved in the k-mer buffers, and will redistribute into other buckets after extension.
      for (i = end; i < bucketCount[n]; i++) bucket[n][i - num] = bucket[n][i];
      bucketCount[n] -= num;
      BUCKET_SHRINK(n);
      end = beg;
      // If no seed extension is allowed at all, quit without adding anything to the extension table
      // or generating INTERVAL records
      if (configHeader->maxSeedBases <= seedLen) continue;
      // Stats for EXTEND primary seeds
      FRQ = (numWithPri > 63 ? 63 : numWithPri);
      stats.extendFreqHist[FRQ] += num;
      // Allocate memory for seed extension IDs at various extension lengths
      uint32_t extendIds[maxLen + 1];
      int      extendInc[maxLen + 1];
      // Populate the primary seed length into the summary records for all seeds
      for (i = 0; i < num; i++) {
        extendSummaryRecs[i].extendLen = seedLen;
        extendSummaryRecs[i].depth     = 0;
      }
      // To prepare the seed extension table, pre-walk the optimal seed extension tree, depth first
      int32_t top = i;
      PRINTIF(DEBUG_EXTEND >= 2, ("Seed extension table entries:\n"));
      for (i = 0; i < num; i = top) {
        len       = extendSummaryRecs[i].extendLen;
        freq      = extendTreeRecs[len][i].frequency;
        top       = i + freq;
        int depth = extendSummaryRecs[i].depth;
        // Walk deeper as long as there is a child (another extension step)
        while (1) {
          // Exit deeper-walk unless another extension step
          nextLen = extendTreeRecs[len][i].extendLen;
          if (!nextLen) break;
          // Push all seeds in the interval to their next extension length
          for (j = i; j < top; j++) {
            extendSummaryRecs[j].extendLen = nextLen;
            extendSummaryRecs[j].depth     = depth;
          }
          // Update extension length and seed frequency to the first child
          len  = nextLen;
          freq = extendTreeRecs[nextLen][i].frequency;
          top  = i + freq;
        }
        // Now we have reached a leaf node, with some terminal seed interval
        uint32_t        numExtendHitRecs = ctx->numExtendHitRecs;
        extend_hit_t    eh;
        kmerToExtend_t* kp = &kmerBuf[i];
        // Record the interval start in the summary record for the leaf node's first seed
        extendSummaryRecs[i].intervalBeg = numExtendHitRecs;
        PRINTIF(DEBUG_EXTEND >= 2, ("---- Leaf node i=%d, freq=%d, len=%d, depth=%d\n", i, freq, len, depth));
        // Populate hits for this leaf node into the extension table interval
        for (j = i; j < top; j++) {
          liftIdx = kp->liftIdx;
          // Locate the counter for the number of extension table entries in the current bin of 256 buckets
          uint32_t  extTabBucketBin  = extIndexBase + (n >> EXTTAB_INDEX_BUCKET_BITS);
          uint32_t* binExtendHitRecs = &(extIndexRecs[extTabBucketBin].length);

          rec = kp->rec;
          // Replace with real record from the special bucket if applicable
          int special = HASH_OPC_IS_SPEC(rec.general.opcode);
          if (special) rec.hit.pos = specialBucket[rec.hit.pos & HASH_REC_SPECIAL_MASK].hit.pos;
          eh.pos        = rec.hit.pos;
          eh.rc         = rec.hit.rc;
          eh.lift_code  = kp->liftCode;
          eh.lift_group = kp->liftCode == LIFT_CODE_NONE ? 0 : kp->liftGroup;
          if (configHeader->maxMultBaseSeeds > 0) {
            // Flag special-bucket hits (palindromes) and masked seed positions
            // (multi-base codes) for literal encoding
            eh.literal = writeCompFile && (special || kp->masked);
          } else {
            eh.literal = writeCompFile && special;
          }
          PRINTIF(
              DEBUG_EXTEND >= 2,
              ("Extension Hit [%u]: pos=%u, rc=%d, lift_code=%d",
               numExtendHitRecs,
               eh.pos,
               eh.rc,
               eh.lift_code));
          PRINTIF(DEBUG_EXTEND >= 2 && eh.lift_code, (", lift_group=%d", eh.lift_group));
          PRINTIF(DEBUG_EXTEND >= 2, ("  (#%d)\n", j));
          // For extension table compression, save a record describing how to locate and populate this
          // extended hit record (not for extension table records flagged for literal encoding)
          if (((configHeader->maxMultBaseSeeds > 0) && writeCompFile && !eh.literal) ||
              ((configHeader->maxMultBaseSeeds == 0) && writeCompFile && !special)) {
            uint32_t lg = eh.lift_group, liftGroupBits = 0;
            for (; lg; lg = lg >> 1) liftGroupBits++;
            seedPopRec_t* popPtr = &seedPopRecs[eh.pos];
            // Force literal encoding if the slotOffset or liftGroup fields don't have enough bits
            if ((eh.lift_group >> SEED_POP_REC_LIFT_GROUP_BITS) ||
                (*binExtendHitRecs >> SEED_POP_REC_SLOT_OFFSET_BITS)) {
              eh.literal        = 1;
              popPtr->automatic = 0;
            } else {
              popPtr->automatic  = 1;
              popPtr->extTable   = 1;
              popPtr->liftCode   = eh.lift_code;
              popPtr->slotOffset = *binExtendHitRecs;
              popPtr->bucketBin  = extTabBucketBin;
              popPtr->liftGroup  = eh.lift_group;
            }
          }
          extendHitRecs[numExtendHitRecs++] = eh;
          if (writeCompFile) (*binExtendHitRecs)++;
          kp++;
        }
        // Erase hit flags written
        while (kmerLiftHitCount) kmerLiftHitFlag[kmerLiftHitList[--kmerLiftHitCount]] = 0;
        // Record the interval end in the summary record for the leaf node's last seed
        extendSummaryRecs[top - 1].intervalEnd = numExtendHitRecs;
        ctx->numExtendHitRecs                  = numExtendHitRecs;
      }
      uint32_t numHits = 0, numMasked = 0;
      uint32_t extIdBin;
      uint64_t bucketIndex;
      // Reset to the primary seed length in the summary records for all seeds
      for (i = 0; i < num; i++) {
        extendSummaryRecs[i].extendLen = seedLen;
        extendSummaryRecs[i].depth     = 0;
      }
      // Array to stash base offsets that extension intervals are relative to, for each seed length
      uint32_t intervalBase[MAX_EXTENDED_LENGTH + 1];
      intervalBase[seedLen] = 0;
      COUNT_CYCLES(cyclesExtendIntervals);
      // Walk the optimal seed extension tree, depth first
      PRINTIF(DEBUG_EXTEND, ("Tree of seed extensions: (%u seeds)\n", num));
      for (i = 0; i < num; i += freq) {
        len               = extendSummaryRecs[i].extendLen;
        freq              = extendTreeRecs[len][i].frequency;
        int      depth    = extendSummaryRecs[i].depth;
        int      incr     = 0;
        uint32_t extendId = 0;
        // Walk deeper as long as there is a child (another extension step).  Terminate immediately
        // if a masked node is reached, extending further than associated reference positions can handle.
        while (!extendTreeRecs[len][i].masked) {
          // If at the root, EXTEND records go in the present bucket,
          // using the same match bits as the HIT records we are processing
          if (depth == 0) {
            bucketIndex               = n;
            rec.qword                 = 0;
            rec.match_bits.match_bits = kmerBuf[0].rec.match_bits.match_bits;
            // Take low hash bits as a new extension ID bin
            extIdBin = (kmerBuf[0].rec.hit.hash_bits) & (numExtIDs - 1);
          } else {
            // Retrieve the extension increment and ID for this length
            extendId = extendIds[len];
            incr     = extendInc[len];
            // Form a secondary hash key, starting with the extension ID shifted up 24 bits
            uint64_t secKey = (uint64_t)extendId << (MAX_SEED_EXTENSION_INCR << 1);
            // Grab up to 14 bases for this incremental extension from the K-mer record,
            // which are interleaved from right and left sides of the primary seed,
            // and insert them growing left and right of the bit 14:13 boundary.
            int extStart = len - incr - seedLen;
            rc           = kmerBuf[i].rec.hit.rc;
            for (j = 0; j < incr; j += 2) {
              secKey |= ((uint64_t)kmerBuf[i].extension[extStart + j] - 4) << (MAX_SEED_EXTENSION_INCR + j);
              secKey |= ((uint64_t)kmerBuf[i].extension[extStart + j + 1] - 4)
                        << (MAX_SEED_EXTENSION_INCR - j - 2);
            }
            // Hash the key
            uint64_t secHash;
            crcHash64(secCrcInit, &secKey, &secHash);
            // Compute target bucket address and hash thread ID
            uint64_t bucketAddr = qwordExtractBits(secHash, byteAddrStart, addrBits) * squeezeRatio;
            uint32_t threadId   = qwordExtractBits(bucketAddr, ADDR_THREAD_ID_START, THREAD_ID_BITS);
            bucketIndex         = bucketAddr >> HASH_BUCKET_BYTES_LOG2;
            // Prep a general format hash record with this information
            rec.qword               = 0;
            rec.matchable.ex        = 1;
            rec.matchable.hash_bits = secHash;
            rec.matchable.thread_id = threadId;
            // Take low hash bits as a new extension ID bin
            extIdBin = secHash & (numExtIDs - 1);
          }
          // Generate an EXTEND record if another extension step
          nextLen = extendTreeRecs[len][i].extendLen;
          if (nextLen) {
            // Seed extensions are given unique IDs to hash with extended bases, because the ID can be fewer
            // bits than the original seed bases we would otherwise have to re-hash, enabling a longer
            // extension. Part of each ID are low bits from the base seed hash (stored in the hashBits field
            // of the primary record), and some are a separate ID field in the EXTEND record. To keep IDs
            // unique, we have to put a different ID field for each combination of low hash bits. So we
            // consider the low hash bits to define a bin, and track the ID fields issued so far in each bin.
            // Above, we populated extIdBin with the proper number of low hash bits indicating the bin.
            // Increment the seed extension ID count in the bin corresponding to the low hash bits.
            if (++extensionIDs[extIdBin] >> HASH_RECORD_EXT_ID_BITS != 0) {
              stats.extendIdBinPctHist[100] = 1;
              sprintf(
                  errMsg,
                  "Too many seed extensions occurred.  Please increase ht-seed-len and/or ht-max-seed-freq.");
              goto buildHashTablesError;
            }
            // Construct the full extension ID, concatenating the low hash bits in extIdBin
            // with the ID field in the EXTEND record
            uint32_t extendId = (extIdBin << HASH_RECORD_EXT_ID_BITS) | extensionIDs[extIdBin];
            // Save as the extension ID for the next extension length
            extendIds[nextLen] = extendId;
            // Extension increment
            int incr           = nextLen - len;
            extendInc[nextLen] = incr;
            // Generate an iterative EXTEND record in the hash table for this seed interval
            if (DEBUG_EXTEND) {
              for (k = seedLen; k < len; k++) printf(" ");
            }
            PRINTIF(
                DEBUG_EXTEND, ("%d: #%u-%u EXTEND #%d to %d\n", len, i, i + freq - 1, depth + 1, nextLen));
            // Generate an iterative EXTEND record in the hash table for this seed interval
            rec.extend.extend_id  = extendId;
            rec.extend.extend_len = incr;
            rec.extend.opcode     = HASH_OPC_EXTEND;
            rec.extend.rs         = 0;
            rec.extend.al         = 0;
            rec.extend.rf         = 0;
            BUCKET_ADD_REC(bucketIndex, rec);
            // Stats
            stats.extendIncrWtdSum += incr * freq;
            stats.countKmerExtendIncr += freq;
            stats.extendIncrHist[incr] += freq;
          }
          // Generate 1-3 INTERVAL records referencing this node's match interval
          // within the seed extension table for this hash table chunk.
          // Retrieve the match interval endpoints from the summary records for the
          // first and last seeds in this node.
          uint32_t intervalBeg = extendSummaryRecs[i].intervalBeg;
          uint32_t intervalEnd = extendSummaryRecs[i + freq - 1].intervalEnd;
          uint32_t intervalLen = intervalEnd - intervalBeg;
          // Save the beginning as the interval base for the next seed length
          intervalBase[nextLen] = intervalBeg;
          // Make the current interval relative to the current base
          uint32_t relativeBeg = intervalBeg - intervalBase[len];
          uint32_t extraLifts  = 0;

        // For primary seeds (depth 0) we need maximal room to encode the interval start,
        // because later a large offset will be added to adjust for previous hash table chunks
          uint32_t adjustBeg = (depth == 0 ? 0x7FFFFFFF : relativeBeg);
          // Decide what combination of 1-3 INTERVAL formats to use, out of 5 available, or one HIT record
          int use_sle = 0, use_sl0 = 0, use_sl1 = 0, use_st = 0, use_ln = 0, use_hit = 0;
          if (intervalLen == 1) {
            use_hit = 1;
          } else if (extraLifts) {
            use_sle = 1;
            use_st  = (adjustBeg > 0xFF);
          } else if (adjustBeg <= 0x7FFF && intervalLen <= 0x1FF) {
            use_sl0 = 1;
          } else {
            use_st = 1;
            if (intervalLen <= 0xFFFF) {
              use_sl1 = 1;
            } else {
              use_ln  = 1;
              use_sle = (adjustBeg > 0xFFFFFF || intervalLen > 0xFFFFFF);
            }
          }
          // Generate the planned INTERVAL records or HIT record
          if (use_sl0) {
            rec.interval_sl0.opcode = HASH_OPC_INTERVAL_SL;
            rec.interval_sl0.start  = relativeBeg;
            rec.interval_sl0.length = intervalLen;
            rec.interval_sl0.fmt    = 0;
            BUCKET_ADD_REC(bucketIndex, rec);
          }
          if (use_sl1) {
            rec.interval_sl1.opcode = HASH_OPC_INTERVAL_SL;
            rec.interval_sl1.start  = relativeBeg >> 24;
            rec.interval_sl1.length = intervalLen;
            rec.interval_sl1.fmt    = 1;
            BUCKET_ADD_REC(bucketIndex, rec);
          }
          if (use_sle) {
            rec.interval_sle.opcode  = HASH_OPC_INTERVAL_SLE;
            rec.interval_sle.start   = relativeBeg >> (use_st * 24);
            rec.interval_sle.length  = (intervalLen - extraLifts) >> (use_ln * 24);
            rec.interval_sle.exlifts = extraLifts;
            rec.interval_sle.msb     = use_st;
            BUCKET_ADD_REC(bucketIndex, rec);
          }
          if (use_st) {
            rec.interval_st.opcode = HASH_OPC_INTERVAL_S;
            rec.interval_st.start  = relativeBeg;
            rec.interval_st.carry  = 0;
            BUCKET_ADD_REC(bucketIndex, rec);
          }
          if (use_ln) {
            rec.interval_ln.opcode = HASH_OPC_INTERVAL_L;
            rec.interval_ln.length = intervalLen;
            rec.interval_ln.rsvd   = 0;
            BUCKET_ADD_REC(bucketIndex, rec);
          }
          if (use_hit) {
            // All extended seeds are populated in the seed extension table, with "automatic" format
            // used for compression.  When a tree leaf node corresponds to just a single hit, we could have
            // used an INTERVAL record to point to that hit in the extension table, but in the same space
            // as an INTERVAL record, we instead populate a single HIT record directly in the hash table,
            // duplicating the same hit also appearing in the extension table.  This is preferable for
            // run time performance, because the mapper avoids an addional access to the extension table.
            // But since this hit is already "automatic"-coded into the extension table, this HIT record
            // needs to be encoded as a "literal" in compression.  To indicate that, move this record
            // into the special bucket if not already there, replacing it with an indirect reference.
            rec.hit.pos = kmerBuf[i].rec.hit.pos;
            rec.hit.rc  = kmerBuf[i].rec.hit.rc;
            if (!HASH_OPC_IS_SPEC(rec.general.opcode)) {
              SPECIAL_ADD_REC(rec);
              rec.hit.pos = HASH_REC_POS_SPECIAL | prevSpecialCount;
            }
            BUCKET_ADD_REC(bucketIndex, rec);
          }
          // Quit tree descent if no further extension
          if (!nextLen) break;
          // Push all seeds in the interval to their next extension length
          depth++;
          for (j = 0; j < freq; j++) {
            extendSummaryRecs[i + j].extendLen = nextLen;
            // Track extension count
            extendSummaryRecs[i + j].depth = depth;
          }
          // Update extension length and seed frequency to the first child
          len  = nextLen;
          freq = extendTreeRecs[nextLen][i].frequency;
        }
        // Now we have reached a leaf node, with some terminal seed interval.
        uint32_t freqWithLift = freq + extendTreeRecs[len][i].liftExtraFreq;
        if (extendTreeRecs[len][i].masked) {
          if (DEBUG_EXTEND) {
            for (k = seedLen; k < len; k++) printf(" ");
            printf("%d: #%u-%u MASKED (%d)\n", len, i, i + freq - 1, freq);
            numMasked += freq;
          }
        } else {
          if (DEBUG_EXTEND) {
            for (k = seedLen; k < len; k++) printf(" ");
            printf("%d: #%u-%u HITs (%d)\n", len, i, i + freq - 1, freq);
            numHits += freq;
          }
          // Stats
          if (len > seedLen) {
            FRQ = freqWithLift > 63 ? 63 : freqWithLift;
            stats.countExtendKmer += freq;
            stats.kmerExtendFreqSum += (uint64_t)freqWithLift * freq;
            stats.passExtendSeed += freq;
            stats.extendLenWtdSum += len * freq;
            stats.extendLengthHist[len] += freq;
            stats.extendStepsWtdSum += (uint64_t)depth * freq;
            stats.extendStepsHist[depth] += freq;
            if (freqWithLift <= maxSeedFreq) {
              stats.hitFreqHist[1][FRQ] += freq;
              stats.countSeedHits[1] += freq;
              stats.secKmerHitFreqSum += (uint64_t)freqWithLift * freq;
            }
          }
        }
      }
      PRINTIF(DEBUG_EXTEND, ("numHits=%u, numMasked=%u\n", numHits, numMasked));
      COUNT_CYCLES(cyclesExtendConstruct);
    }
#ifndef LOCAL_BUILD
    WatchDogCheckin(watchdogID);
#endif
  }
  PRINTIF(DEBUG_EXTEND, ("    Processed %lu buckets - %lu seeds extended\n", n, stats.passExtendSeed));

  PRINTIF(DEBUG_CHAIN, ("  Converting long probing into hash chains\n"));
  PRINTIF(DEBUG_CHAIN, ("Analyzing bucket occupancy and probe distances\n"));
  // Analyze bucket occupancy and probe distances
  int     thread;
  int32_t residue, maxResidue;
  maxBucketCount = 0;
  for (n = 0; n < numBuckets; n++)
    if (bucketCount[n] > maxBucketCount) maxBucketCount = bucketCount[n];
  free(tempBucket);
  tempBucket = malloc(maxBucketCount * 4 * 8);
  if (!tempBucket) {
    sprintf(errMsg, "Failed to allocate temporary buckets for sort");
    goto buildHashTablesError;
  }
  COUNT_CYCLES(cyclesOverhead);
  for (n = 0; n < numBuckets; n++) {
    CHECK_ABORT;
    // Hash records are assigned to one of N 'threads', where N is the number of records per bucket,
    // by the log2(N) hash bits just below the bucket address.
    // Count how many records in each thread
    memset(threadCounts, 0, HASH_THREADS << 2);
    for (i = 0; i < bucketCount[n]; i++) {
      thread = bucket[n][i].hit.thread_id & BUCKET_THREAD_MASK;
      threadCounts[thread]++;
    }
    // Make a temporary copy of the bucket with records extended by 2 qwords holding sort keys.
    // In each triplet, QW[2] = primary sort key, QW[1] = secondary sort key, QW[0] = hash record.
    int      liftInc = 0, seenAlt = 0, seenExt = 0;
    uint32_t liftGroup = 0;
    for (i = j = 0; i < bucketCount[n]; i++) {
      hashrec_t myRec = bucket[n][i];
      // Replace with real record from the special bucket if applicable
      if (HASH_OPC_IS_SPEC(myRec.general.opcode))
        myRec.hit.pos = specialBucket[myRec.hit.pos & HASH_REC_SPECIAL_MASK].hit.pos;
      // The sorting goal is:
      //    First groups of EXTEND records, and associated match interval references:
      //        (EXTEND, INTERVAL..), (EXTEND, INTERVAL..), ...
      //    Then normal hits:
      //        HIT, HIT, HIT, ...
      //  Within each class, groups are ordered from shortest thread to longest thread.
      //  Hash record "match bits" are used to group associated EXTEND, and INTERVAL records.
      //  Each "liftover group" of hits (a sequence of 1+ ALT hits followed by 1 PRI hit) is kept together
      //  with the PRI last. The CRC32's of normal HIT positions are used to randomize their order.  The
      //  reason for this randomization is that iterative seed extension repair may obtain more than the
      //  maximum number of HIT records, and when we cut off further hits, we want which hits are discarded to
      //  be unbiased.
      //
      // Qword 0: Hash record (original)
      tempBucket[j++] = bucket[n][i].qword;
      // Qword 1: Hash record (replaced)
      tempBucket[j++] = myRec.qword;
      // Qword 2: RC flag, position hash, opcode type, non-ALT flag, liftover group
      // TODO: remove uint32_t myHash = 0; // Don't need to randomize order anymore
      int myOpc = myRec.general.opcode;
      seenExt |= (myOpc == HASH_OPC_EXTEND);
      int isHit  = HASH_OPC_IS_HIT(myOpc);
      int myType = isHit ? 15 : (myOpc & 0xF) - 1;
      int altHit = isHit && (myRec.hit.pos >= refAltSeed);
      int priHit = isHit && (myRec.hit.pos < refAltSeed);
      // Increment liftover group on first ALT in a row, or after PRI following ALTs
      if (liftInc || (altHit && !seenAlt)) {
        liftGroup++;
        seenAlt = 0;
        liftInc = 0;
      }
      liftInc = priHit & seenAlt;
      seenAlt |= altHit;
      tempBucket[j++] = ((uint64_t)liftGroup << 39) | ((uint64_t)priHit << 38) | ((uint64_t)myType << 33) |
                        (uint64_t)myRec.hit.rc;
      // Qword 3: match bits, thread count, and seen-extend flag to sort (EXTEND, INTERVAL..) groups first
      thread          = bucket[n][i].hit.thread_id & BUCKET_THREAD_MASK;
      tempBucket[j++] = ((uint64_t)seenExt << 63) | ((uint64_t)threadCounts[thread] << 30) |
                        (uint64_t)myRec.match_bits.match_bits;
    }
    // Re-sort the buckets using the sort keys we attached
    qsort(tempBucket, bucketCount[n], 32, hashRecCompareCustom);
    // Replace records in original bucket
    for (i = j = 0; i < bucketCount[n]; i++, j += 4) bucket[n][i].qword = tempBucket[j];
    COUNT_CYCLES(cyclesBucketSort);
    // Scan in reverse and flag the last (first seen) record for each thread.
    memset(threadFlags, 0, HASH_THREADS);
    for (i = bucketCount[n] - 1; i >= 0; i--) {
      thread = bucket[n][i].hit.thread_id & BUCKET_THREAD_MASK;
      if (threadFlags[thread]) {
        bucket[n][i].hit.lf = 0;
      } else {
        threadFlags[thread] = 1;
        bucket[n][i].hit.lf = 1;
      }
    }
    // Stats
    VAL = bucketCount[n];
    VAL = (VAL > 63 ? 63 : VAL);
    stats.bucketLevelEscHist[VAL]++;
    // Scan forward, and for each thread, move the first record up front.
    // By policy, if there are any records in a bucket's thread, at least one must
    // appear in the target bucket, not get pushed to other buckets.  That way if
    // no record of a thread being searched for is present, probing can stop immediately.
    hashrec_t tempRec;
    memset(threadFlags, 0, HASH_THREADS);
    for (i = j = 0; i < bucketCount[n]; i++) {
      thread = bucket[n][i].hit.thread_id & BUCKET_THREAD_MASK;
      if (threadFlags[thread]) continue;
      threadFlags[thread] = 1;
      // Move up front
      tempRec = bucket[n][i];
      for (k = i; k > j; k--) bucket[n][k] = bucket[n][k - 1];
      bucket[n][j] = tempRec;
      j++;
    }
    // Scan forward from the current bucket to compute the maximum probe distance
    residue = 0;
    for (i = 0; i < probeLimit; i++) {
      // Probing wraps at a certain chunk size, e.g. 4KB
      m = ((n + i) % wrapBuckets) + (n - (n % wrapBuckets));
      residue += bucketCount[m] - HASH_RECORDS_PER_BUCKET;
      if (residue <= 0) {
        bucketProbeDist[n] = i;
        break;
      }
    }
    if (i == probeLimit) bucketProbeDist[n] = i;
    // Scan backward from the current bucket to compute the expected occupancy here
    maxResidue = residue = bucketCount[n];
    for (i = 1; i < probeLimit; i++) {
      // Probing wraps at a certain chunk size, e.g. 4KB
      m = ((wrapBuckets + n - i) % wrapBuckets) + (n - (n % wrapBuckets));
      residue += bucketCount[m] - HASH_RECORDS_PER_BUCKET;
      maxResidue = (residue > maxResidue ? residue : maxResidue);
    }
    bucketOccupancy[n] = (maxResidue > HASH_RECORDS_PER_BUCKET ? HASH_RECORDS_PER_BUCKET : maxResidue);

#ifndef LOCAL_BUILD
    WatchDogCheckin(watchdogID);
#endif
    COUNT_CYCLES(cyclesBucketOrganize);
  }
  if (DEBUG_CHAIN) {
    printf("Bucket probe distances histogram:\n");
    printHistogram(stdout, bucketProbeDistHist, probeLimit + 1, 2, 1);
    printf("Bucket occupancy histogram:\n");
    printHistogram(stdout, bucketOccupancyHist, HASH_RECORDS_PER_BUCKET + 1, 2, 0);
  }
  // Loop through 'chain blocks', regions corresponding to the reach of chain pointers
  uint32_t chainBlockBuckets = numBuckets < CHAIN_BLOCK_BUCKETS ? numBuckets : CHAIN_BLOCK_BUCKETS;
  uint32_t numChainBlocks    = (numBuckets + (chainBlockBuckets - 1)) / chainBlockBuckets;
  uint32_t chainBlock;
  stats.chainLengthHist[0] = numBuckets;
  stats.probeLengthHist[0] = numBuckets;
  for (chainBlock = 0; chainBlock < numChainBlocks; chainBlock++) {
    uint32_t bucketStart = chainBlock * chainBlockBuckets, bucketEnd = bucketStart + chainBlockBuckets;
    bucketEnd = (bucketEnd > numBuckets ? numBuckets : bucketEnd);
    // Count how many buckets of each probe distance and occupancy within this chain block
    memset(bucketProbeDistHist, 0, (MAX_PROBES + 1) * 8);
    memset(bucketOccupancyHist, 0, (HASH_RECORDS_PER_BUCKET + 1) * 8);
    for (n = bucketStart; n < bucketEnd; n++) {
      bucketProbeDistHist[bucketProbeDist[n]]++;
      bucketOccupancyHist[bucketOccupancy[n]]++;
    }
    // Resolve long probes by chaining from long-probing buckets to low-occupancy buckets,
    // using special escape records to hold chain pointers
    int      targetProbeDist, targetOccupancy, chainRoom, chainLength;
    uint32_t chainBucket = bucketStart;
    for (targetOccupancy = 0; bucketOccupancyHist[targetOccupancy] == 0; targetOccupancy++)
      ;
    // Probe distances 1+ can potentially benefit from chaining
    for (targetProbeDist = MAX_PROBES; targetProbeDist >= 1; targetProbeDist--) {
      PRINTIF(
          DEBUG_CHAIN,
          ("Target probe distance = %u, target occupancy = %u\n", targetProbeDist, targetOccupancy));
      // Skip unrepresented probe distances
      if (bucketProbeDistHist[targetProbeDist] == 0) continue;
      // Hunt for buckets probing this far
      for (n = bucketStart; n < bucketEnd; n++) {
        if (bucketProbeDist[n] != targetProbeDist) continue;
        PRINTIF(
            DEBUG_CHAIN,
            ("  Found source bucket %u, count %u, probe distance = %u\n",
             n,
             bucketCount[n],
             targetProbeDist));
        // Advance to more heavily occupied chain bucket targets if the chain wouldn't fit
        for (; targetOccupancy < HASH_RECORDS_PER_BUCKET - 1; targetOccupancy++, chainBucket = bucketStart) {
          chainRoom   = HASH_RECORDS_PER_BUCKET - targetOccupancy - 1;
          chainLength = (bucketCount[n] - (HASH_RECORDS_PER_BUCKET - 1) + chainRoom - 1) / chainRoom;
          PRINTIF(
              DEBUG_CHAIN,
              ("    targetOccupancy=%d, chainRoom=%u, chainLength=%u\n",
               targetOccupancy,
               chainRoom,
               chainLength));
          if (chainLength <= bucketOccupancyHist[targetOccupancy]) break;
        }
        if (targetOccupancy >= HASH_RECORDS_PER_BUCKET - 1) break;
        // Unless we are forced to chain, check if this bucket would benefit from chaining,
        // meaning a chain comprising buckets with the current room available would be no longer.
        // If the same length, tilt toward chaining, because the filter should improve performance.
        if (targetProbeDist < MAX_PROBES && targetProbeDist < chainLength) {
          PRINTIF(DEBUG_CHAIN, ("    No chaining benefit for this bucket\n"));
          continue;
        }
        i = HASH_RECORDS_PER_BUCKET - 1;
        PRINTIF(DEBUG_CHAIN, ("    Exporting records %u-%u to a chain\n", i, bucketCount[n] - 1));
        hashrec_t *linkFromRec = &bucket[n][i], *linkToRec;
        uint32_t   chainMask, chainList, chainListLen, listVal, listTemp, maskVal, useList, listDiffs;
        uint8_t*   chainListBytes = (uint8_t*)&chainList;
        int        chainFirst     = 1;
        // Find chain buckets for all the overflow records
        for (; i < bucketCount[n]; chainFirst = 0) {
          // Find the next chain bucket
          for (; chainBucket < bucketEnd; chainBucket++) {
            if (bucketOccupancy[chainBucket] == targetOccupancy) break;
          }
          PRINTIF(
              DEBUG_CHAIN, ("    Found chain bucket %u, count %u\n", chainBucket, bucketCount[chainBucket]));
          // Plan a chain link record in the target bucket
          j = bucketCount[chainBucket];
          PRINTIF(DEBUG_CHAIN, ("    Chain link record at bucket[%u][%d]\n", chainBucket, j));
          linkToRec        = &bucket[chainBucket][j++];
          linkToRec->qword = HASH_REC_CHAIN_TERM_QWORD;
          stats.chainLenSum++;
          // for hash table compression, stash count of previously chained records
          int prevChained            = i - (HASH_RECORDS_PER_BUCKET - 1);
          linkToRec->chain.chain_pad = prevChained > 63 ? 63 : prevChained;
          // Follow with expatriated records
          useList      = 1;
          chainMask    = 0;
          chainList    = 0;
          chainListLen = 0;
          int availCount =
              HASH_RECORDS_PER_BUCKET - (bucketOccupancy[chainBucket] - bucketCount[chainBucket]);
          for (k = i; k < bucketCount[n]; k++) {
            if (j < availCount) {
              PRINTIF(DEBUG_CHAIN, ("    Moving bucket[%u][%d] to bucket[%u][%d]\n", n, k, chainBucket, j));
              bucket[chainBucket][j] = bucket[n][k];
              j++;
              i++;
            }
            // Keep track of hash bit sub-segments seen in remaining records to chain,
            // to encode into the previous chain record to filter accesses.  We remember
            // a list of up to 4 distinct 8-bit segments, and in case that overflows, a mask
            // representing 32 5-bit segments.
            maskVal = hashRecExtractHashBits(bucket[n][k], 0, CHAIN_MASK_HASH_BITS);
            if (useList) {
              listVal = hashRecExtractHashBits(bucket[n][k], 0, CHAIN_LIST_HASH_BITS);
              for (listDiffs = 0; listDiffs < chainListLen; listDiffs++) {
                // Keep list sorted
                if (listVal < chainListBytes[listDiffs]) {
                  listTemp                  = listVal;
                  listVal                   = chainListBytes[listDiffs];
                  chainListBytes[listDiffs] = listTemp;
                  continue;
                }
                if (listVal == chainListBytes[listDiffs]) break;
              }
              if (listDiffs == chainListLen) {
                if (chainListLen < 4)
                  chainListBytes[chainListLen++] = listVal;
                else
                  useList = 0;
              }
            }
            chainMask |= ((uint32_t)1 << maskVal);
          }
          bucketCount[chainBucket] = j;
          PRINTIF(DEBUG_CHAIN, ("    Now bucketCount[%u]=%d\n", chainBucket, bucketCount[chainBucket]));
          // Fill in the chain link record
          // On the first chain step, construct the beginning of chain record (we didn't do this earlier
          // because this slot was occupied by a hit record until we just finished moving it)
          if (chainFirst) {
            linkFromRec->qword        = 0;
            linkFromRec->chain.opcode = (useList ? HASH_OPC_CHAIN_BEG_LIST : HASH_OPC_CHAIN_BEG_MASK);
          } else {
            linkFromRec->chain.opcode = (useList ? HASH_OPC_CHAIN_CON_LIST : HASH_OPC_CHAIN_CON_MASK);
          }
          // Put the list or mask filter into the previous record
          if (useList) {
            // If less than 4 listed, copy the last value to the remaining entries
            listVal = chainListBytes[chainListLen - 1];
            for (k = chainListLen; k < 4; k++) chainListBytes[k] = listVal;
            linkFromRec->chain.filter = chainList;
          } else {
            linkFromRec->chain.filter = chainMask;
          }
          // Chain the previous link record here
          PRINTIF(DEBUG_CHAIN, ("    Linking previous to bucket %u\n", chainBucket));
          linkFromRec->chain.chain_ptr = chainBucket & (CHAIN_BLOCK_BUCKETS - 1);
          linkFromRec                  = linkToRec;
          // Tally this bucket done
          bucketOccupancyHist[targetOccupancy]--;
          PRINTIF(
              DEBUG_CHAIN,
              ("    Buckets left with targetOccupancy=%d: %u\n",
               targetOccupancy,
               bucketOccupancyHist[targetOccupancy]));
          chainBucket++;
        }
        // Forget the relocated records
        bucketCount[n] = HASH_RECORDS_PER_BUCKET;
        PRINTIF(DEBUG_CHAIN, ("    Now bucketCount[%u]=%d\n", n, bucketCount[n]));
        // Stats
        VAL = (chainLength > 63 ? 63 : chainLength);
        stats.chainLengthHist[VAL]++;
        stats.chainLengthHist[0]--;
        VAL = (targetProbeDist > 63 ? 63 : targetProbeDist);
        stats.probeDistChainedHist[VAL]++;
        stats.passChains++;
        stats.countChains++;
        stats.countChainBuckets += chainLength;
      }
    }
#ifndef LOCAL_BUILD
    WatchDogCheckin(watchdogID);
#endif
  }
  COUNT_CYCLES(cyclesBucketChain);
  PRINTIF(DEBUG_CHAIN, ("Done with chaining\n"));
  PRINTIF(DEBUG_CHAIN, ("    Processed %u buckets - built %u hash chains\n", numBuckets, stats.passChains));

  // Compression stuff
  int           literal, chainOffset, popOffset, isChain;
  uint8_t       litByte = 0;
  seedPopRec_t* popPtr  = NULL;

  PRINTIF(DEBUG_PACK, ("  Mapping records to physical buckets\n"));
  // Consider one "wrap block" at a time
  uint32_t wrap;
  for (wrap = 0; wrap < numWraps; wrap++) {
    CHECK_ABORT;
    PRINTIF(DEBUG_PACK, ("\nPacking wrap %u\n", wrap));
    uint32_t bucketStart = wrap * wrapBuckets, bucketEnd = bucketStart + wrapBuckets;
    memset(recUsed, 0, wrapRecords);
    memset(wrapSpecials, 0, wrapRecords);
    // Initialize physical record buffer to empty records
    for (i = 0; i < wrapRecords; i++) physRecs[i].qword = HASH_REC_EMPTY_QWORD;
    // Now we need to smear the remaining records into their physical buckets
    // or successive buckets, within the probing limit.
    int dist;
    // Increment a distance value from zero to the limit, and at each stage
    // place remaining records exactly that distance ahead.  There is no
    // inter-bucket contention with this method.  It probably could be
    // optimized by giving threads near their probing limit higher placement
    // priority, but then there is complex contention to work out.
    memset(recsPlaced, 0, wrapBuckets);
    for (dist = 0; dist <= probeLimit; dist++) {
      for (n = bucketStart; n < bucketEnd; n++) {
        CHECK_ABORT;
        m = ((n + dist) % wrapBuckets) * HASH_RECORDS_PER_BUCKET;
        for (j = 0, i = recsPlaced[n % wrapBuckets]; i < bucketCount[n];) {
          // If a chain escape record is found, place it and successive records
          // at the end of the physical bucket.  The number of local and probing
          // entries here has been precomputed, and should fit before the chain escape.
          int opc = bucket[n][i].general.opcode;
          if (HASH_OPC_IS_CHAIN(opc)) {
            j = HASH_RECORDS_PER_BUCKET - (bucketCount[n] - i);
            if (j < 0 || recUsed[m + j]) {
              sprintf(
                  errMsg,
                  "Failed packing chain records: n=%u, opc=0x%X, bucketCount=%u, i=%d, j=%d",
                  n,
                  opc,
                  bucketCount[n],
                  i,
                  j);
              goto buildHashTablesError;
            }
          }
          for (; j < HASH_RECORDS_PER_BUCKET; j++) {
            k = m + j;
            if (!recUsed[k]) {
              hashrec_t rec = bucket[n][i];
              // Replace with real record from the special bucket if applicable
              if (HASH_OPC_IS_SPEC(rec.general.opcode)) {
                rec.hit.pos = specialBucket[rec.hit.pos & HASH_REC_SPECIAL_MASK].hit.pos;
                // Flag that the special bucket was used (this flag is the whole purpose of the special
                // bucket)
                wrapSpecials[k] = 1;
              }
              // Also flag HIT records special for literal encoding if their seeds
              // overlapped masked reference positions (multi-base codes)
              else if ((configHeader->maxMultBaseSeeds > 0) && HASH_OPC_IS_HIT(rec.general.opcode)) {
                pos = floor(rec.hit.pos * refSeedInterval);
                wrapSpecials[k] |= SEQ_MASK_FAIL(pos, seedLen);
              }
              physRecs[k] = rec;
              recUsed[k]  = 1;
              // Record probe distance, for compression
              wrapOffsets[k] =
                  ((dist << HASH_RECORDS_PER_BUCKET_LOG2) | (k & ((1 << HASH_RECORDS_PER_BUCKET_LOG2) - 1))) +
                  1;
              PRINTIF(
                  DEBUG_PACK >= 2,
                  ("    Source bucket 0x%X, record %u (0x%016llX): placed in physical bucket 0x%X, offset %d\n",
                   n,
                   i,
                   bucket[n][i].qword,
                   m / HASH_RECORDS_PER_BUCKET + bucketStart,
                   j));
              j++;
              i++;
              // Stats
              if (i == bucketCount[n]) {
                stats.probeLenSum += dist;
                stats.probeLengthHist[dist]++;
                stats.probeLengthHist[0]--;
              }
              break;
            }
          }
          if (j == HASH_RECORDS_PER_BUCKET) break;
        }
        PRINTIF(
            DEBUG_PACK && (recsPlaced[n % wrapBuckets] != i || !dist),
            ("  Bucket %u: cnt=%u:  put %d in bkt %u, dist %d\n",
             n,
             bucketCount[n],
             i - recsPlaced[n % wrapBuckets],
             m / HASH_RECORDS_PER_BUCKET + bucketStart,
             dist));
        recsPlaced[n % wrapBuckets] = i;
        if (dist == probeLimit && i != bucketCount[n]) {
          sprintf(errMsg, "Hash table needs to be larger.  Failed to place all records for bucket %u.", n);
          goto buildHashTablesError;
        }
      }
    }
    for (i = 0; i < wrapBuckets; i++) {
      CHECK_ABORT;
      int occ = HASH_RECORDS_PER_BUCKET;
      for (j = 0; j < HASH_RECORDS_PER_BUCKET; j++) {
        int opc = physRecs[(i << HASH_RECORDS_PER_BUCKET_LOG2) + j].general.opcode;
        if (HASH_OPC_IS_HIT(opc))
          stats.countHitRec++;
        else if (HASH_OPC_IS_CHAIN(opc))
          stats.countChainRec++;
        else if (opc == HASH_OPC_EMPTY) {
          stats.countEmptyRec++;
          occ--;
        } else if (HASH_OPC_IS_INTVL(opc))
          stats.countIntervalRec++;
        else if (opc == HASH_OPC_EXTEND)
          stats.countExtendRec++;
        else {
          sprintf(errMsg, "Unrecognized hash record opcode 0x%02X", opc);
          goto buildHashTablesError;
        }
      }
      stats.bucketLevelMapHist[occ]++;
    }
    // Free any virtual buckets
    for (n = bucketStart; n < bucketEnd; n++) BUCKET_CLEAR(n);
    // Copy the constructed wrap back into the physical record space
    memcpy(bucket[bucketStart], physRecs, wrapBytes);
    // In the copy, clear chain record offset fields that we overloaded for compression
    for (i = 0; i < wrapRecords; i++) {
      hashrec_t* rp = bucket[bucketStart] + i;
      if (HASH_OPC_IS_CHAIN(rp->general.opcode)) rp->chain.chain_pad = 0;
    }
    COUNT_CYCLES(cyclesBucketWrite);
    // Analyze wrap of records for compression
    if (writeCompFile) {
      for (i = 0; 1; i++) {
        if (!(i & 7) && i) {
          *litFlags++ = litByte;
          litByte     = 0;
          if (i == wrapRecords) break;
        }
        if (!(i & ((1 << HASH_RECORDS_PER_BUCKET_LOG2) - 1)))
          chainOffset = 0;
        else if (chainOffset)
          chainOffset++;
        hashrec_t rec = physRecs[i];
        if (rec.general.opcode == HASH_OPC_EMPTY) continue;
        isChain = HASH_OPC_IS_CHAIN(rec.general.opcode);
        if (isChain)
          chainOffset = rec.chain.chain_pad == 63 ? 100 : rec.chain.chain_pad + HASH_RECORDS_PER_BUCKET - 1;
        // For compression, most HIT records can be encoded in "automatic" mode, which is cheap (less so for
        // extended seeds). Non-HITs must be encoded as literal records, which is more expensive.  Also, some
        // HITs violating the assumptions of automatic encoding must be literals: those flagged as "special"
        // via temporary storage in the special bucket above (non-matching liftover hits, reversed
        // palindromes), and null-position hits used to terminate liftover groups.
        literal =
            !HASH_OPC_IS_HIT(rec.general.opcode) || wrapSpecials[i] || !rec.hit.pos || chainOffset >= 100;
        if (!literal) {
          popOffset = chainOffset ? chainOffset : wrapOffsets[i];
          popPtr    = &seedPopRecs[rec.hit.pos];
          if (popPtr->automatic)
            literal = 1;
          else {
            popPtr->automatic  = 1;
            popPtr->slotOffset = popOffset - 1;
          }
        }
        if (literal) litByte |= 1 << (i & 7);
      }
      COUNT_CYCLES(cyclesBucketCompress);
    }
#ifndef LOCAL_BUILD
    WatchDogCheckin(watchdogID);
#endif
  }

buildHashTablesError:
  // Return largest seed extension ID count in any bin for compression
  ctx->maxExtendIds = 0;
  for (i = 0; i < numExtIDs; i++)
    if (extensionIDs[i] > ctx->maxExtendIds) ctx->maxExtendIds = extensionIDs[i];
  // Null-pad 256 bytes onto the seed extension table (we allocated enough for this)
  extend_hit_t nullExtendHit = {0};
  for (i = 0; i < 32; i++) ctx->extendHitRecs[ctx->numExtendHitRecs + i] = nullExtendHit;
  // Finalize stats
  for (i = 0; i < numExtIDs; i++)
    stats.extendIdBinPctHist[(int)(100 * (double)extensionIDs[i] / (1 << HASH_RECORD_EXT_ID_BITS) + 0.5)]++;
  stats.chainOrProbeLengthHist[0] = numBuckets;
  for (i = 1; i < 64; i++) {
    stats.chainOrProbeLengthHist[i] = stats.chainLengthHist[i] + stats.probeLengthHist[i];
    stats.chainOrProbeLengthHist[0] -= (stats.chainLengthHist[i] + stats.probeLengthHist[i]);
  }
  ctx->stats = stats;

  // Cleanup
  free(specialBucket);
  free(bucketOccupancy);
  free(bucketProbeDist);
  free(extensionIDs);
  free(kmerBuf);
  free(kmerLiftHitFlag);
  free(kmerLiftHitList);
  free(secCrcInit);
  free(randCrcInit);
  free(tempBucket);
  free(extendSummaryRecs);
  for (n = 0; n <= MAX_EXTENDED_LENGTH; n++) free(extendTreeRecs[n]);
  for (n = 0; n < numBuckets; n++) BUCKET_CLEAR(n);

#ifndef LOCAL_BUILD
  WatchDogDeregister(watchdogID);
#endif

  COUNT_CYCLES(cyclesOverhead);
  if (errMsg[0])
    return errMsg;
  else
    return NULL;
}

typedef struct {
  int          threadId;
  int          numThreads;
  int          seedLen;
  int          addrBits;
  int          binBits;
  int          numChunks;
  int          allocFail;
  void*        crcInit;
  uint64_t     maxPos;
  uint32_t     numHashes;
  uint32_t     numPalindromes;
  uint32_t     numRecs;
  uint32_t     numMultBasePos;
  uint32_t     numMultBaseSeeds;
  uint32_t     maxMultBaseSeeds;
  uint32_t     liftMatchSeedInt;
  uint8_t*     refSeq;
  uint8_t*     refMask;
  uint8_t*     refCode;
  uint8_t*     altMatches;
  uint64_t     refAltStart;
  double       refSeedInterval;
  double       squeezeRatio;
  uint64_t     chunkBuckets;
  uint64_t     passBucketStart;
  uint64_t     passBucketEnd;
  uint32_t**   bucketAlloc;
  uint32_t**   bucketCount;
  uint8_t**    bucketLocks;
  hashrec_t*** bucket;
} hashThreadCtx_t;

// Worker thread to hash blocks of reference seeds
void* hashThreadMultBase(void* arg)
{
  // Extract arguments
  hashThreadCtx_t* ctx              = arg;
  int              threadId         = ctx->threadId;
  int              numThreads       = ctx->numThreads;
  int              seedLen          = ctx->seedLen;
  int              addrBits         = ctx->addrBits;
  int              binBits          = ctx->binBits;
  int              numChunks        = ctx->numChunks;
  void*            crcInit          = ctx->crcInit;
  uint64_t         maxPos           = ctx->maxPos;
  uint8_t*         refSeq           = ctx->refSeq;
  uint8_t*         refMask          = ctx->refMask;
  uint8_t*         refCode          = ctx->refCode;
  uint8_t*         altMatches       = ctx->altMatches;
  uint64_t         refAltStart      = ctx->refAltStart;
  uint32_t         maxMultBaseSeeds = ctx->maxMultBaseSeeds;
  uint32_t         liftMatchSeedInt = ctx->liftMatchSeedInt;
  double           refSeedInterval  = ctx->refSeedInterval;
  double           squeezeRatio     = ctx->squeezeRatio;
  uint64_t         chunkBuckets     = ctx->chunkBuckets;
  uint64_t         passBucketStart  = ctx->passBucketStart;
  uint64_t         passBucketEnd    = ctx->passBucketEnd;
  uint32_t**       bucketAlloc      = ctx->bucketAlloc;
  uint32_t**       bucketCount      = ctx->bucketCount;
  uint8_t**        bucketLocks      = ctx->bucketLocks;
  hashrec_t***     bucket           = ctx->bucket;

  uint64_t bucketIndex;
  uint16_t chunkIndex;

  uint64_t priSeedMask = ((uint64_t)1 << (seedLen * 2)) - 1;
  uint32_t priMaskMask = (uint32_t)(((uint64_t)1 << seedLen) - 1);
  uint64_t refBinMask  = binBits ? ((0ULL - 1) << binBits) : 0;
  uint64_t fwSeed, rcSeed = 0, hashKey, hash;
  int      i, j, k;

  int      multNumPos, code;
  uint8_t  multNum[32], multIdx[32];
  uint32_t multNumSeeds, seqMask, num;
  uint32_t altMatchVec;
  uint64_t altPos;
  uint64_t multFwSeeds[32][4], multRcSeeds[32][4], fwBaseSeed, rcBaseSeed;

  uint32_t numHashes        = 0;
  uint32_t numPalindromes   = 0;
  uint32_t numRecs          = 0;
  uint32_t numMultBasePos   = 0;
  uint32_t numMultBaseSeeds = 0;
  uint32_t index            = threadId;
  uint64_t pos              = floor(index * refSeedInterval);
#define SAVED_RECS 16
  int       savedChunkIndex[SAVED_RECS];
  uint64_t  savedBucketIndex[SAVED_RECS];
  hashrec_t savedRec[SAVED_RECS];
  int       savedCount = 0;
  int       fail       = 0;
  hashrec_t rec;
  rec.qword = 0;

#ifndef LOCAL_BUILD
  unsigned long watchdogID = WatchDogRegister("_hash_thread");
#endif
  // Loop through all positions in our window
  for (; pos <= maxPos; pos = floor((index += numThreads) * refSeedInterval)) {
    // Grab the k-mer, and its reverse complement
    fwBaseSeed = GET_SEQ_BITS(pos) & priSeedMask;
    revComp(&rcBaseSeed, &fwBaseSeed, seedLen);
    // Skip alt contig positions where the entire primary seed matches its liftover sequence exactly
    if (pos >= refAltStart) {
      altPos      = pos - refAltStart;
      altMatchVec = (*(uint64_t*)(&altMatches[altPos >> 3]) >> (altPos & 7)) & priMaskMask;
      if (altMatchVec == priMaskMask && (index % liftMatchSeedInt)) continue;
    }
    // Initialize multi-base seed processing if any masked position is hit
    seqMask = (uint32_t)GET_SEQ_MASK(pos) & priMaskMask;
    if (seqMask) {
      // Skip position if multi-base codes not supported
      if (refCode == NULL) continue;
      multNumPos   = 0;
      multNumSeeds = 1;
      for (i = 0, j = seedLen - 1; seqMask; i++, j--, seqMask = seqMask >> 1) {
        // Process masked positions, which encode a number of nucleotides other than 1
        if ((seqMask & 1)) {
          // Grab the 4-bit code, containing a '1' for each matching base
          code = GET_BASE_CODE(pos + i);
          // Number of '1's
          num = baseCodeNumBases[code];
          // Multiply the total number of seeds
          multNumSeeds *= num;
          // Early exit if we hit a zero-base code, or if the total number of seeds will be too great
          if (multNumSeeds == 0 || multNumSeeds > maxMultBaseSeeds) break;
          // Initialize the iterator for this multi-base position
          multNum[multNumPos] = num;
          multIdx[multNumPos] = 0;
          // Erase this position from the base seeds
          fwBaseSeed &= ~((uint64_t)3 << (i << 1));
          rcBaseSeed &= ~((uint64_t)3 << (j << 1));
          // List the nucleotides in this code, shifted into their seed positions
          for (k = 0; k < num; k++) {
            multFwSeeds[multNumPos][k] = (uint64_t)baseCodeBaseList[code][k] << (i << 1);
            multRcSeeds[multNumPos][k] = ((uint64_t)baseCodeBaseList[code][k] ^ 3) << (j << 1);
          }
          multNumPos++;
        }
      }
      // Skip position if number of seeds will be zero or too many
      if (multNumSeeds == 0 || multNumSeeds > maxMultBaseSeeds) continue;
      numMultBaseSeeds += multNumSeeds;
      numMultBasePos++;
    }
    // Only 1-base codes: initialize trivial iteration
    else {
      multNumPos = 0;
    }
    // Iteration to generate seeds combining choices for multi-base codes
    while (1) {
      // Grab the baseline k-mer, and its reverse complement
      fwSeed = fwBaseSeed;
      rcSeed = rcBaseSeed;
      // Inject current choices for multi-base positions
      for (i = 0; i < multNumPos; i++) {
        fwSeed |= multFwSeeds[i][multIdx[i]];
        rcSeed |= multRcSeeds[i][multIdx[i]];
      }
      // Use whichever is numerically smaller as the seed
      int useRc = (rcSeed < fwSeed);
      hashKey   = useRc ? rcSeed : fwSeed;
      // Append the reference bin if applicable
      hashKey |= (pos & refBinMask) << KEY_ANCHOR_OFFSET;
      // Hash the seed
      crcHash64(crcInit, &hashKey, &hash);
      numHashes++;
      // Extract the bucket byte address from the left end of the hash,
      // scaled down by the squeeze ratio for non-power-of-two table sizes
      uint64_t bucketAddr  = qwordExtractBits(hash, HASH_BYTE_ADDR_START, addrBits) * squeezeRatio;
      uint64_t bucketIndex = bucketAddr >> HASH_BUCKET_BYTES_LOG2;
      int      palindrome  = (rcSeed == fwSeed);
      numPalindromes += palindrome;
      // Skip if not in the address range for this pass
      if (bucketIndex < passBucketStart || bucketIndex >= passBucketEnd) continue;
      int chunkIndex = (bucketIndex / chunkBuckets) % numChunks;
      bucketIndex %= chunkBuckets;
      rec.hit.pos       = index;
      rec.hit.hash_bits = hash;
      rec.hit.thread_id = qwordExtractBits(bucketAddr, ADDR_THREAD_ID_START, THREAD_ID_BITS);
      // Generate seeds in both orientations for palindromes
      do {
        rec.hit.rc                   = useRc;
        savedChunkIndex[savedCount]  = chunkIndex;
        savedBucketIndex[savedCount] = bucketIndex;
        savedRec[savedCount]         = rec;
        savedCount++;
        __builtin_prefetch(&bucket[chunkIndex][bucketIndex]);
        __builtin_prefetch(&bucketAlloc[chunkIndex][bucketIndex]);
        __builtin_prefetch(&bucketCount[chunkIndex][bucketIndex]);
        __builtin_prefetch(&bucketLocks[chunkIndex][bucketIndex]);
        if (savedCount == SAVED_RECS) {
          for (j = 0; j < savedCount; j++) {
            chunkIndex  = savedChunkIndex[j];
            bucketIndex = savedBucketIndex[j];
            while (!__sync_bool_compare_and_swap(&bucketLocks[chunkIndex][bucketIndex], 0, 1))
              ;
            THREAD_BUCKET_ADD_REC(chunkIndex, bucketIndex, savedRec[j]);
#if defined(_TARGET_PPC_)
            //
            // If you run into an issue on x86 where the hash tables don't match when built multiple times
            // likely may happen when we build with more threads per core
            // uncomment out this line
            // __asm__ __volatile__ ("lfence ::: "memory");
            __asm__ __volatile__("sync" ::: "memory");
#endif
            __atomic_clear(&bucketLocks[chunkIndex][bucketIndex], __ATOMIC_SEQ_CST);
            if (fail) {
              ctx->allocFail = 0;
#ifndef LOCAL_BUILD
              WatchDogDeregister(watchdogID);
#endif
              return NULL;
            }
          }
          savedCount = 0;
        }
        numRecs++;
      } while (palindrome && !useRc++);
#ifndef LOCAL_BUILD
      WatchDogCheckin(watchdogID);
#endif
      // Advance iteration to next combination of choices for multi-base positions
      for (i = 0; i < multNumPos; i++) {
        if (++multIdx[i] < multNum[i]) {
          break;
        } else {
          multIdx[i] = 0;
        }
      }
      // Done when no more combinations were found
      if (i == multNumPos) break;
    }
  }
  for (j = 0; j < savedCount; j++) {
    chunkIndex  = savedChunkIndex[j];
    bucketIndex = savedBucketIndex[j];
    while (!__sync_bool_compare_and_swap(&bucketLocks[chunkIndex][bucketIndex], 0, 1))
      ;
    THREAD_BUCKET_ADD_REC(chunkIndex, bucketIndex, savedRec[j]);
#if defined(_TARGET_PPC_)
    __asm__ __volatile__("sync" ::: "memory");
#endif
    __atomic_clear(&bucketLocks[chunkIndex][bucketIndex], __ATOMIC_SEQ_CST);
  }
  // Report counts
  ctx->numHashes        = numHashes;
  ctx->numPalindromes   = numPalindromes;
  ctx->numRecs          = numRecs;
  ctx->numMultBasePos   = numMultBasePos;
  ctx->numMultBaseSeeds = numMultBaseSeeds;
  ctx->allocFail        = 0;
#ifndef LOCAL_BUILD
  WatchDogDeregister(watchdogID);
#endif
  return NULL;
}

// Worker thread to hash blocks of reference seeds
// TODO This calls the multi base version of the function when (maxMultBaseSeeds > 0). The next time the hash
// table version increases, remove this version of the function and only use the multi base version.
void* hashThread(void* arg)
{
  // Extract arguments
  hashThreadCtx_t* ctx = arg;

  uint32_t maxMultBaseSeeds = ctx->maxMultBaseSeeds;
  if (maxMultBaseSeeds > 0) {
    return hashThreadMultBase(arg);
  }

  int          threadId        = ctx->threadId;
  int          numThreads      = ctx->numThreads;
  int          seedLen         = ctx->seedLen;
  int          addrBits        = ctx->addrBits;
  int          binBits         = ctx->binBits;
  int          numChunks       = ctx->numChunks;
  void*        crcInit         = ctx->crcInit;
  uint64_t     maxPos          = ctx->maxPos;
  uint8_t*     refSeq          = ctx->refSeq;
  uint8_t*     refMask         = ctx->refMask;
  double       refSeedInterval = ctx->refSeedInterval;
  double       squeezeRatio    = ctx->squeezeRatio;
  uint64_t     chunkBuckets    = ctx->chunkBuckets;
  uint64_t     passBucketStart = ctx->passBucketStart;
  uint64_t     passBucketEnd   = ctx->passBucketEnd;
  uint32_t**   bucketAlloc     = ctx->bucketAlloc;
  uint32_t**   bucketCount     = ctx->bucketCount;
  uint8_t**    bucketLocks     = ctx->bucketLocks;
  hashrec_t*** bucket          = ctx->bucket;

  uint64_t bucketIndex;
  uint16_t chunkIndex;

  uint64_t priSeedMask = ((uint64_t)1 << (seedLen * 2)) - 1;
  uint64_t refBinMask  = binBits ? ((0ULL - 1) << binBits) : 0;
  uint64_t fwSeed, rcSeed = 0, hashKey, hash;
  int      j;

  uint32_t numHashes      = 0;
  uint32_t numPalindromes = 0;
  uint32_t numRecs        = 0;
  uint32_t index          = threadId;
  uint64_t pos            = floor(index * refSeedInterval);
#define SAVED_RECS 16
  int       savedChunkIndex[SAVED_RECS];
  uint64_t  savedBucketIndex[SAVED_RECS];
  hashrec_t savedRec[SAVED_RECS];
  int       savedCount = 0;
  int       fail       = 0;
  hashrec_t rec;
  rec.qword = 0;

#ifndef LOCAL_BUILD
  unsigned long watchdogID = WatchDogRegister("_hash_thread");
#endif
  // Loop through all positions in our window
  for (; pos <= maxPos; pos = floor((index += numThreads) * refSeedInterval)) {
    // Skip if not all hashable bases
    if (SEQ_MASK_FAIL(pos, seedLen)) continue;
    // Grab the k-mer, and its reverse complement
    fwSeed = GET_SEQ_BITS(pos) & priSeedMask;
    revComp(&rcSeed, &fwSeed, seedLen);
    // Use whichever is numerically smaller as the seed
    int useRc = (rcSeed < fwSeed);
    hashKey   = useRc ? rcSeed : fwSeed;
    // Append the reference bin if applicable
    hashKey |= (pos & refBinMask) << KEY_ANCHOR_OFFSET;
    // Hash the seed
    crcHash64(crcInit, &hashKey, &hash);
    numHashes++;
    // Extract the bucket byte address from the left end of the hash,
    // scaled down by the squeeze ratio for non-power-of-two table sizes
    uint64_t bucketAddr  = qwordExtractBits(hash, HASH_BYTE_ADDR_START, addrBits) * squeezeRatio;
    uint64_t bucketIndex = bucketAddr >> HASH_BUCKET_BYTES_LOG2;
    int      palindrome  = (rcSeed == fwSeed);
    numPalindromes += palindrome;
    // Skip if not in the address range for this pass
    if (bucketIndex < passBucketStart || bucketIndex >= passBucketEnd) continue;
    int chunkIndex = (bucketIndex / chunkBuckets) % numChunks;
    bucketIndex %= chunkBuckets;
    rec.hit.pos       = index;
    rec.hit.hash_bits = hash;
    rec.hit.thread_id = qwordExtractBits(bucketAddr, ADDR_THREAD_ID_START, THREAD_ID_BITS);
    // Generate seeds in both orientations for palindromes
    do {
      rec.hit.rc                   = useRc;
      savedChunkIndex[savedCount]  = chunkIndex;
      savedBucketIndex[savedCount] = bucketIndex;
      savedRec[savedCount]         = rec;
      savedCount++;
      __builtin_prefetch(&bucket[chunkIndex][bucketIndex]);
      __builtin_prefetch(&bucketAlloc[chunkIndex][bucketIndex]);
      __builtin_prefetch(&bucketCount[chunkIndex][bucketIndex]);
      __builtin_prefetch(&bucketLocks[chunkIndex][bucketIndex]);
      if (savedCount == SAVED_RECS) {
        for (j = 0; j < savedCount; j++) {
          chunkIndex  = savedChunkIndex[j];
          bucketIndex = savedBucketIndex[j];
          while (!__sync_bool_compare_and_swap(&bucketLocks[chunkIndex][bucketIndex], 0, 1))
            ;
          THREAD_BUCKET_ADD_REC(chunkIndex, bucketIndex, savedRec[j]);
#if defined(_TARGET_PPC_)
          //
          // If you run into an issue on x86 where the hash tables don't match when built multiple times
          // likely may happen when we build with more threads per core
          // uncomment out this line
          // __asm__ __volatile__ ("lfence ::: "memory");
          __asm__ __volatile__("sync" ::: "memory");
#endif
          __atomic_clear(&bucketLocks[chunkIndex][bucketIndex], __ATOMIC_SEQ_CST);
          if (fail) {
            ctx->allocFail = 0;
#ifndef LOCAL_BUILD
            WatchDogDeregister(watchdogID);
#endif
            return NULL;
          }
        }
        savedCount = 0;
      }
      numRecs++;
    } while (palindrome && !useRc++);
#ifndef LOCAL_BUILD
    WatchDogCheckin(watchdogID);
#endif
  }
  for (j = 0; j < savedCount; j++) {
    chunkIndex  = savedChunkIndex[j];
    bucketIndex = savedBucketIndex[j];
    while (!__sync_bool_compare_and_swap(&bucketLocks[chunkIndex][bucketIndex], 0, 1))
      ;
    THREAD_BUCKET_ADD_REC(chunkIndex, bucketIndex, savedRec[j]);
#if defined(_TARGET_PPC_)
    __asm__ __volatile__("sync" ::: "memory");
#endif
    __atomic_clear(&bucketLocks[chunkIndex][bucketIndex], __ATOMIC_SEQ_CST);
  }
  // Report counts
  ctx->numHashes      = numHashes;
  ctx->numPalindromes = numPalindromes;
  ctx->numRecs        = numRecs;
  ctx->allocFail      = 0;
#ifndef LOCAL_BUILD
  WatchDogDeregister(watchdogID);
#endif
  return NULL;
}

// Returns NULL on success, or error message
char* buildHashTable(
    hashTableConfig_t* config,
    uint8_t*           refSeq,
    uint64_t*          refCodeHist,
    uint8_t*           refMask,
    uint8_t*           refCode,
    uint8_t*           altMatches

)
{
  int32_t     i, j, k, n, pass, chunksComplete = 0;
  static char errMsg[256]  = "";
  uint64_t    refSeqLen    = config->hdr->refSeqLen;
  uint64_t    hashDigest   = 0;
  uint64_t    extTabDigest = 0;
  uint32_t    maxExtendIds = 0;
#if defined(_TARGET_PPC_)
  pthread_mutex_t lock;
#endif

  // Statistics
  uint64_t           bytesWritten = 0, totExtTabRecs = 0;
  uint32_t           countKmer = 0, countPal = 0, totalRecs = 0, countMultBasePos = 0, countMultBaseSeeds = 0;
  buildThreadStats_t stats     = {0};
  FILE*              statsFile = NULL;

  // Some global or table-specific stuff
  int seedLen = config->hdr->priSeedBases;

  void* crcInit = crcHash64Init(config->hdr->priCrcBits, config->hdr->priCrcPoly);
  if (!crcInit) return "Failure initializing primary hash function\n";

  // Grab config params
  int   addrBits      = config->hdr->tableAddrBits;
  int   maxThreads    = config->maxThreads;
  int   maxGB         = config->maxGB;
  int   writeCompFile = config->writeCompFile;
  int   writeHashFile = config->writeHashFile;
  char* hashFname     = writeHashFile ? config->hashFname : NULL;
  char* compFname     = writeCompFile ? config->compFname : NULL;
  char* extTabFname   = writeHashFile ? config->extTabFname : NULL;
  char* statsFname    = config->statsFname;

  // Derive hash record layout parameters
  double   squeezeRatio   = (float)config->hdr->tableSize64ths / 64;  // This floating point value is exact
  uint64_t tableBytes     = ((uint64_t)1 << addrBits) * squeezeRatio;
  int      bucketAddrBits = addrBits - HASH_BUCKET_BYTES_LOG2;
  uint64_t numBuckets     = ((uint64_t)1 << bucketAddrBits) * squeezeRatio;
  uint64_t numRecords     = ((uint64_t)1 << (bucketAddrBits + HASH_RECORDS_PER_BUCKET_LOG2)) * squeezeRatio;
  uint64_t maxPos         = refSeqLen - seedLen;

  // Compute independent chunk size (between 0.5 and 1.0 GB)
  int      tableChunks  = (addrBits < INDEPENDENT_ADDR_BITS ? 1 : 1 << addrBits - INDEPENDENT_ADDR_BITS);
  uint64_t chunkBytes   = tableBytes / tableChunks;
  uint64_t chunkBuckets = numBuckets / tableChunks;
  uint64_t chunkRecords = numRecords / tableChunks;
  // Calculate how many passes with how many threads to process these chunks
  int      maxBuildThreads = (maxThreads > maxGB ? maxGB : maxThreads);
  int      numPasses       = (tableChunks + maxBuildThreads - 1) / maxBuildThreads;
  int      numBuildThreads = (tableChunks + numPasses - 1) / numPasses;
  uint64_t passBytes       = numBuildThreads * chunkBytes;
  uint64_t passBuckets     = numBuildThreads * chunkBuckets;
  // Pick a multiple of numBuildThreads for chunks to keep in memory at once
  int passesPerHashRun = maxGB / numBuildThreads;
  if (passesPerHashRun > numPasses) passesPerHashRun = numPasses;
  int hashChunks = passesPerHashRun * numBuildThreads;
  printf("Constructing %s hash table", bytesReadable(tableBytes));
  if (numPasses == 1)
    printf(" in 1 pass\n");
  else
    printf(" in %d passes, %s each\n", numPasses, bytesReadable(passBytes));
  fflush(stdout);
  // Context record for compressed output
  writeCompHashTableCtx_t compCtx = {0};

  // Allocate empty buckets for each thread
  uint32_t*   bucketAlloc[hashChunks];
  uint32_t*   bucketCount[hashChunks];
  uint8_t*    bucketLocks[hashChunks];
  uint8_t*    chunkLitFlags[hashChunks];
  hashrec_t** bucket[hashChunks];
  hashrec_t*  physRecords[hashChunks];
  for (i = 0; i < hashChunks; i++) {
    bucketAlloc[i] = malloc(4 * chunkBuckets);
    bucketCount[i] = malloc(4 * chunkBuckets);
    bucketLocks[i] = calloc(1, chunkBuckets);
    bucket[i]      = malloc(sizeof(hashrec_t*) * chunkBuckets);
    physRecords[i] = malloc(chunkBytes);
    if (!bucketAlloc[i] || !bucketCount[i] || !bucket[i] || !physRecords[i]) {
      sprintf(errMsg, "Failed to allocate empty buckets");
      goto mainBuildHashTableError;
    }
    chunkLitFlags[i] = NULL;
    if (writeCompFile) {
      chunkLitFlags[i] = malloc(chunkRecords / 8);
      if (!chunkLitFlags[i]) {
        sprintf(errMsg, "Failed to allocate compression flags");
        goto mainBuildHashTableError;
      }
    }
  }

  // Allocate compressed format information arrays
  uint32_t       numPopRecs        = (uint32_t)ceil(refSeqLen / config->hdr->refSeedInterval);
  seedPopRec_t*  seedPopRecs       = NULL;
  uint32_t       chunkExtIndexRecs = chunkBuckets >> EXTTAB_INDEX_BUCKET_BITS;
  uint32_t       numExtIndexRecs   = chunkExtIndexRecs * tableChunks;
  extIndexRec_t* extIndexRecs      = NULL;
  if (writeCompFile) {
    seedPopRecs  = calloc(numPopRecs, sizeof(seedPopRec_t));
    extIndexRecs = calloc(numExtIndexRecs, sizeof(extIndexRec_t));
    if (!seedPopRecs || !extIndexRecs) {
      sprintf(errMsg, "Failed to allocate arrays for compression");
      goto mainBuildHashTableError;
    }
  }

  // Open the hash table output file
  FILE* hashFile   = NULL;
  FILE* extTabFile = NULL;
  if (writeHashFile) {
    hashFile = fopen(hashFname, "wb");
    if (!hashFile) {
      sprintf(errMsg, "Cannot open hash table file");
      goto mainBuildHashTableError;
    }
    extTabFile = fopen(extTabFname, "wb");
    if (!extTabFile) {
      sprintf(errMsg, "Cannot open extension table file");
      goto mainBuildHashTableError;
    }
  }

  // Open the compressed output file
  if (writeCompFile) {
    if (writeCompHashTableCtxOpen(&compCtx, compFname)) {
      sprintf(errMsg, "Cannot open compressed hash table output");
      goto mainBuildHashTableError;
    }
    if (writeCompHashTableHeader(&compCtx, config->hdr)) {
      sprintf(errMsg, "Failure writing compressed hash table (header)");
      goto mainBuildHashTableError;
    }
  }

  // To limit memory utilization, take multiple passes through the reference.  Hash all K-mers on each pass,
  // but only accumulate a subset of them in buckets, to construct a portion of the hash table.
  int abortBuildThreads = 0;
  for (pass = 0; pass < numPasses; pass++) {
    if (numPasses > 1) {
      printf("Pass %d\n", pass + 1);
      fflush(stdout);
    }
    // Don't necessarily run seed hashing on every pass.  If enough memory is allowed, we can save hash
    // records for more thread chunks than the number of build threads, and use them for multiple passes
    if (pass % passesPerHashRun == 0) {
      int      hashPass        = pass / passesPerHashRun;
      uint64_t hashPassBuckets = passBuckets * passesPerHashRun;
      uint64_t passBucketStart = hashPass * hashPassBuckets;
      uint64_t passBucketEnd   = (hashPass + 1) * hashPassBuckets;
      uint32_t countRecs       = 0;
      countKmer                = 0;
      countPal                 = 0;

      // Point virtual buckets into physical record space to start
      for (i = 0; i < hashChunks; i++) {
        for (j = 0; j < chunkBuckets; j++) {
          bucketCount[i][j] = 0;
          bucketAlloc[i][j] = HASH_RECORDS_PER_BUCKET;
          bucket[i][j]      = physRecords[i] + j * HASH_RECORDS_PER_BUCKET;
        }
      }

      // Spawn worker threads to hash reference seeds
      int numHashThreads = maxThreads, allocFail = 0;
      int remainPasses = (pass + passesPerHashRun > numPasses ? numPasses - pass : passesPerHashRun);
      printf(
          "  Spawning %d thread%s to hash reference seeds", numHashThreads, numHashThreads == 1 ? "" : "s");
      if (remainPasses == 1)
        printf(" for this pass\n");
      else
        printf(" for %d passes\n", remainPasses);
      fflush(stdout);

#if defined(_TARGET_PPC_)
      // Initialize the mutex
      //
      pthread_mutex_init(&lock, NULL);
#endif

      pthread_t hashThreads[numHashThreads];
      memset(hashThreads, 0, numHashThreads * sizeof(pthread_t));
      hashThreadCtx_t hashCtx[numHashThreads];
      for (i = 0; i < numHashThreads; i++) {
        hashCtx[i].threadId         = i;
        hashCtx[i].numThreads       = numHashThreads;
        hashCtx[i].seedLen          = seedLen;
        hashCtx[i].addrBits         = addrBits;
        hashCtx[i].binBits          = config->hdr->anchorBinBits;
        hashCtx[i].crcInit          = crcInit;
        hashCtx[i].maxPos           = maxPos;
        hashCtx[i].maxMultBaseSeeds = config->hdr->maxMultBaseSeeds;
        hashCtx[i].liftMatchSeedInt = config->hdr->liftMatchSeedInt ? config->hdr->liftMatchSeedInt : 1;
        hashCtx[i].refSeq           = refSeq;
        hashCtx[i].refMask          = refMask;
        hashCtx[i].refCode          = refCode;
        hashCtx[i].altMatches       = altMatches;
        hashCtx[i].refAltStart      = config->hdr->refAltStart;
        hashCtx[i].refSeedInterval  = config->hdr->refSeedInterval;
        hashCtx[i].squeezeRatio     = squeezeRatio;
        hashCtx[i].chunkBuckets     = chunkBuckets;
        hashCtx[i].passBucketStart  = passBucketStart;
        hashCtx[i].passBucketEnd    = passBucketEnd;
        hashCtx[i].numChunks        = hashChunks;
        hashCtx[i].bucketAlloc      = bucketAlloc;
        hashCtx[i].bucketCount      = bucketCount;
        hashCtx[i].bucketLocks      = bucketLocks;
        hashCtx[i].bucket           = bucket;
        pthread_create(&hashThreads[i], NULL, hashThread, &hashCtx[i]);
      }
      // Wait for threads to return
      for (i = 0; i < numHashThreads; i++) {
        pthread_join(hashThreads[i], NULL);
        countKmer += hashCtx[i].numHashes;
        countPal += hashCtx[i].numPalindromes;
        countRecs += hashCtx[i].numRecs;
        totalRecs += hashCtx[i].numRecs;
        countMultBasePos += hashCtx[i].numMultBasePos;
        countMultBaseSeeds += hashCtx[i].numMultBaseSeeds;
        allocFail |= hashCtx[i].allocFail;
      }
      if (allocFail) {
        sprintf(errMsg, "Thread bucket allocation failure");
        goto mainBuildHashTableError;
      }
      printf(
          "    %lu seed positions - hashed %lu k-mers, generated %lu records\n",
          maxPos + 1,
          countKmer,
          countRecs);
      fflush(stdout);
    }

    // Reduce the number of build threads if not enough chunks remain
    int origBuildThreads = numBuildThreads;
    if ((pass + 1) * numBuildThreads > tableChunks) numBuildThreads = tableChunks - pass * numBuildThreads;

    // Spawn worker threads for each chunk in this pass
    if (numBuildThreads == 1)
      printf("  Spawning 1 thread to build a %s hash table chunk\n", bytesReadable(chunkBytes));
    else
      printf(
          "  Spawning %d threads to build %s hash table chunks\n",
          numBuildThreads,
          bytesReadable(chunkBytes));
    fflush(stdout);
    pthread_t        buildThreads[numBuildThreads];
    buildThreadCtx_t buildCtx[numBuildThreads];
    int              chunkOffset = (pass % passesPerHashRun) * origBuildThreads;
    for (i = 0; i < numBuildThreads; i++) {
      int chunk            = i + pass * origBuildThreads;
      buildCtx[i].config   = config;
      buildCtx[i].refSeq   = refSeq;
      buildCtx[i].refMask  = refMask;
      buildCtx[i].refCode  = refCode;
      buildCtx[i].litFlags = chunkLitFlags[i];

      buildCtx[i].seedPopRecs  = seedPopRecs;
      buildCtx[i].extIndexRecs = extIndexRecs;
      buildCtx[i].extIndexBase = chunk * chunkExtIndexRecs;
      buildCtx[i].threadId     = i;
      buildCtx[i].chunkNum     = chunk;
      buildCtx[i].bucketStart  = (chunkOffset + i) * chunkBuckets;
      buildCtx[i].bucketAlloc  = bucketAlloc[chunkOffset + i];
      buildCtx[i].bucketCount  = bucketCount[chunkOffset + i];
      buildCtx[i].bucket       = bucket[chunkOffset + i];
      buildCtx[i].physRecords  = physRecords[chunkOffset + i];
      buildCtx[i].numBuckets   = chunkBuckets;
      buildCtx[i].abort        = &abortBuildThreads;
      pthread_create(&buildThreads[i], NULL, buildHashTableThread, &buildCtx[i]);
    }
#ifndef LOCAL_BUILD
    unsigned long watchdogID = WatchDogRegister("_write_hash_table");
#endif
    // Wait for threads to complete
    for (i = 0; i < numBuildThreads; i++) {
      char* errorMsg;
      pthread_join(buildThreads[i], (void**)&errorMsg);
      // Aggregate stats from this thread
      chunksComplete++;
      uint64_t *srcQword = (uint64_t*)&buildCtx[i].stats, *dstQword = (uint64_t*)&stats;
      for (j = 0; j < sizeof(buildThreadStats_t); j += 8) *dstQword++ += *srcQword++;
      // Track the maximum seed extension ID, for compression
      if (buildCtx[i].maxExtendIds > maxExtendIds) maxExtendIds = buildCtx[i].maxExtendIds;
      // Check status
      if (abortBuildThreads) {
        printf("    Thread %d exited\n", i);
        fflush(stdout);
        goto buildThreadCleanup;
      }
      if (errorMsg) {
        printf("    Thread %d error\n", i);
        fflush(stdout);
        strcpy(errMsg, errorMsg);
        abortBuildThreads = 1;
        goto buildThreadCleanup;
      }
      printf("    Thread %d complete", i);
      fflush(stdout);
#ifndef LOCAL_BUILD
      WatchDogCheckin(watchdogID);
#endif
      // Adjust the start fields of primary-seed INTERVAL records, which are relative to this
      // chunk's portion of the seed extension table, by adding an offset equal to the total
      // seed extension table length of prior chunks
      hashrec_t*    recp             = physRecords[chunkOffset + i];
      hashrec_t*    endp             = recp + chunkRecords;
      extend_hit_t* extendHitRecs    = buildCtx[i].extendHitRecs;
      uint32_t      numExtendHitRecs = buildCtx[i].numExtendHitRecs;
      uint64_t      extendHitBytes   = (uint64_t)numExtendHitRecs * 8;
      uint32_t      sum;
      uint32_t      low24 = totExtTabRecs & 0xFFFFFF;
      uint32_t      high8 = totExtTabRecs >> 24;
      for (; recp < endp; recp++) {
        // Only process primary (ex=0) INTERVAL records
        if (!HASH_OPC_IS_INTVL(recp->general.opcode) | recp->matchable.ex) continue;
        switch (recp->general.opcode) {
        // One of the primary INTERVAL records must be interval_st format; add the low 24
        // offset bits to the 24-bit start field, and save a carry flag on overflow
        case HASH_OPC_INTERVAL_S:
          sum                     = recp->interval_st.start + low24;
          recp->interval_st.start = sum & 0xFFFFFF;
          recp->interval_st.carry = sum >> 24;
          break;
        // Another primary INTERVAL record must be interval_sl1 or interval_sle format;
        // add the high 8 offset bits to the 8-bit start field (same position in both formats)
        case HASH_OPC_INTERVAL_SL:
        case HASH_OPC_INTERVAL_SLE:
          recp->interval_sle.start += high8;
        }
      }
      totExtTabRecs += numExtendHitRecs;
      // Pad the final chunk's extension table records to reach a 256-byte boundary
      // (the build thread padded with 256 bytes of null records)
      if (buildCtx[i].chunkNum == tableChunks - 1) {
        while (totExtTabRecs % 32) {
          totExtTabRecs++;
          numExtendHitRecs++;
          extendHitBytes += 8;
        }
      }
      // Send extend table records flagged literal=1 to the compressed output, and clear those flags
      if (writeCompFile) writeCompHashTableExtTabLiterals(&compCtx, extendHitRecs, numExtendHitRecs);
      // Write the thread's output records to disk
      recp = physRecords[chunkOffset + i];
      if (writeHashFile) {
        printf("; writing %s chunk...", bytesReadable(chunkBytes));
        fflush(stdout);
        if (fwrite(recp, HASH_RECORD_BYTES, chunkRecords, hashFile) != chunkRecords) {
          sprintf(errMsg, "Error writing to hash table file");
          abortBuildThreads = 1;
          goto buildThreadCleanup;
        }
        bytesWritten += chunkBytes;
        if (extendHitBytes) printf(" writing %llu extension table bytes...", extendHitBytes);
        fflush(stdout);
        if (fwrite(extendHitRecs, 8, numExtendHitRecs, extTabFile) != numExtendHitRecs) {
          sprintf(errMsg, "Error writing to extension table file");
          abortBuildThreads = 1;
          goto buildThreadCleanup;
        }
        bytesWritten += extendHitBytes;
        printf("  Total bytes written = %llu", bytesWritten);
      }
      printf("\n");
      fflush(stdout);

#if defined(_TARGET_PPC_)
      pthread_mutex_lock(&lock);
#endif
      for (k = 0; k < chunkRecords; k++, recp++) {
#if defined(LOCAL_BUILD) && defined(__x86_64__)
        __asm__ __volatile__(
            "crc32q\t"
            "(%1), %0"
            : "=r"(hashDigest)
            : "r"((uint64_t*)recp), "0"(hashDigest));
#elif !defined(LOCAL_BUILD)
        hashDigest   = crc32c_hw(hashDigest, (const unsigned char*)recp, 8);
#endif
      }
      for (k = 0; k < numExtendHitRecs; k++, extendHitRecs++) {
#if defined(LOCAL_BUILD) && defined(__x86_64__)
        __asm__ __volatile__(
            "crc32q\t"
            "(%1), %0"
            : "=r"(extTabDigest)
            : "r"((uint64_t*)extendHitRecs), "0"(extTabDigest));
#elif !defined(LOCAL_BUILD)
        extTabDigest = crc32c_hw(extTabDigest, (const unsigned char*)extendHitRecs, 8);
#endif
      }
#if defined(_TARGET_PPC_)
      pthread_mutex_unlock(&lock);
#endif
      if (writeCompFile) {
        if (writeCompHashTableLiterals(&compCtx, buildCtx[i].physRecords, chunkRecords, chunkLitFlags[i])) {
          sprintf(errMsg, "Failure writing compressed hash table (literals)");
          abortBuildThreads = 1;
          goto buildThreadCleanup;
        }
      }
    buildThreadCleanup:
      free(buildCtx[i].extendHitRecs);
      buildCtx[i].extendHitRecs = NULL;
    }

#ifndef LOCAL_BUILD
    WatchDogDeregister(watchdogID);
#endif

    if (abortBuildThreads) goto mainBuildHashTableError;
  }
  printf("Table occupancy %.1f%%\n", (double)(numRecords - stats.countEmptyRec) / numRecords * 100);
  if (totExtTabRecs && config->extTableAlloc)
    printf(
        "Extension table space utilized %lu / %lu records (%.1f%%)\n",
        totExtTabRecs,
        config->extTableAlloc,
        100.0 * totExtTabRecs / config->extTableAlloc);
  if (totExtTabRecs > config->extTableAlloc) {
    printf("Extension table overflow, please retry with --ht-ext-table-alloc=%lu\n", totExtTabRecs);
    sprintf(
        errMsg,
        "Extension table overflow, need more extension table space reserved, retry with --ht-ext-table-alloc=%lu",
        totExtTabRecs);
    goto mainBuildHashTableError;
  }

  // Finish compression
  if (writeCompFile) {
    printf("Compressing hash table...");
    fflush(stdout);
    // Measure how many bits are needed to represent extension table offsets within each of the bucket bins
    uint32_t val, len;
    for (i = 0; i < numExtIndexRecs; i++) {
      val = extIndexRecs[i].length;
      val -= (val != 0);
      for (len = 0; val; val = val >> 1) len++;
      extIndexRecs[i].offsetBits = len;
    }
    // Write extension table index
    if (writeCompHashTableExtIndex(&compCtx, extIndexRecs, numExtIndexRecs, totExtTabRecs)) {
      sprintf(errMsg, "Failure writing compressed hash table (extend table index)");
      goto mainBuildHashTableError;
    }
    // Write automatic hits
    if (writeCompHashTableAutoHits(&compCtx, seedPopRecs, numPopRecs, maxExtendIds, extIndexRecs)) {
      sprintf(errMsg, "Failure writing compressed hash table (automatics)");
      goto mainBuildHashTableError;
    }
    if (bytesWritten)
      printf(
          "  compressed bytes = %llu (%.1f%%)\n",
          compCtx.bitLen >> 3,
          100.0 * (compCtx.bitLen >> 3) / bytesWritten);
    fflush(stdout);
  }

mainBuildHashTableError:

  // Stats file
  statsFile = fopen(statsFname, "w");
  if (!statsFile) goto mainBuildHashTableCleanup;

  // Analyze reference code histogram
  uint64_t countBases[4] = {
      refCodeHist[BASE_A], refCodeHist[BASE_C], refCodeHist[BASE_G], refCodeHist[BASE_T]};
  uint64_t countGC            = refCodeHist[BASE_C] + refCodeHist[BASE_G];
  uint64_t countAT            = refCodeHist[BASE_A] + refCodeHist[BASE_T];
  uint64_t countUnmasked      = countGC + countAT;
  uint64_t countMasked        = config->hdr->refSeqLen - countUnmasked;
  uint64_t origLen            = config->hdr->refSeqLen - refCodeHist[BASE_PAD];
  uint64_t baseSetSizeHist[5] = {
      refCodeHist[BASE_PAD],
      refCodeHist[BASE_A] + refCodeHist[BASE_C] + refCodeHist[BASE_G] + refCodeHist[BASE_T],
      refCodeHist[BASE_M_AC] + refCodeHist[BASE_R_AG] + refCodeHist[BASE_S_CG] + refCodeHist[BASE_W_AT] +
          refCodeHist[BASE_Y_CT] + refCodeHist[BASE_K_GT],
      refCodeHist[BASE_V_ACG] + refCodeHist[BASE_H_ACT] + refCodeHist[BASE_D_AGT] + refCodeHist[BASE_B_CGT],
      refCodeHist[BASE_N]};

#define NZERO(x) ((x == 0) ? 1 : (x))
  // Find max seed extension ID utilization
  int extIdBinPctMax = 0;
  for (i = 0; i <= 100; i++)
    if (stats.extendIdBinPctHist[i]) extIdBinPctMax = i;

  // Write stats
  if (errMsg[0]) {
    fprintf(statsFile, "Build hash table terminated with error:\n");
    fprintf(statsFile, "  %s\n", errMsg);
    fprintf(statsFile, "Statistics below may be incomplete, but are presented to aid in diagnosis.\n");
    fprintf(statsFile, "\n");
  }
  fprintf(statsFile, "Reference sequence:\n");
  fprintf(statsFile, "  Original: %10llu\n", origLen);
  fprintf(statsFile, "  Encoded:  %10llu\n", config->hdr->refSeqLen);
  fprintf(
      statsFile,
      "  Masked:   %10llu  (%.1f%%)\n",
      countMasked,
      (double)countMasked / NZERO(config->hdr->refSeqLen) * 100);
  fprintf(statsFile, "  Unmasked: %10llu\n", countUnmasked);
  for (i = 0; i < 4; i++) fprintf(statsFile, "  %c bases:  %10llu\n", BASE_CHARS[i], countBases[i]);
  fprintf(statsFile, "  GC content: %4.1f%%\n", (double)countGC / NZERO(countUnmasked) * 100);
  fprintf(statsFile, "  IUPAC-IUB Codes:\n");
  fprintf(statsFile, "    0 bases (padding)     : %10llu\n", baseSetSizeHist[0]);
  fprintf(statsFile, "    1 base  (A,C,G,T)     : %10llu\n", baseSetSizeHist[1]);
  fprintf(statsFile, "    2 bases (K,M,R,S,W,Y) : %10llu\n", baseSetSizeHist[2]);
  fprintf(statsFile, "    3 bases (B,D,H,V)     : %10llu\n", baseSetSizeHist[3]);
  fprintf(statsFile, "    4 bases (N)           : %10llu\n", baseSetSizeHist[4]);
  fprintf(statsFile, "\n");
  fprintf(statsFile, "Reference K-mers: (K=%d)\n", seedLen);
  fprintf(statsFile, "  Distinct K-mers:     %10u\n", stats.countUniqKmer);
  fprintf(statsFile, "  K-mer positions:     %10u\n", countKmer);
  fprintf(statsFile, "  Palindromes:         %10u\n", countPal);
  if (config->hdr->maxMultBaseSeeds > 0) {
    fprintf(statsFile, "  Non-ACGT seed pos:   %10u\n", countMultBasePos);
    fprintf(statsFile, "    K-mers generated:  %10u\n", countMultBaseSeeds);
    fprintf(
        statsFile, "    Average K-mers:    %10.4f\n", (double)countMultBaseSeeds / NZERO(countMultBasePos));
  }
  fprintf(statsFile, "  Total K-mer records: %10u\n", totalRecs);
  fprintf(statsFile, "  Thinned out:         %10u\n", stats.countThinnedKmers);
  fprintf(statsFile, "  Populated seeds:     %10u\n", stats.countSeedHits[0] + stats.countSeedHits[1]);
  fprintf(statsFile, "  NOTE: All K-mer frequency stats are w.r.t. reference K-mer positions,\n");
  fprintf(statsFile, "        and hence a K-mer with frequency N is included N times.\n");
  fprintf(
      statsFile,
      "  Average K-mer frequency: %.2f\n",
      (double)stats.kmerFreqSumSquare / NZERO(stats.countRawKmer));
  fprintf(statsFile, "  K-mer frequency histogram:\n");
  printHistogram(statsFile, stats.kmerFreqHist, 64, 4, 1);
  fprintf(statsFile, "  Log2 K-mer frequency histogram:\n");
  printHistogram(statsFile, stats.kmerLogFreqHist, 32, 4, 0);
  fprintf(statsFile, "\n");
  fprintf(statsFile, "Alt contig K-mer positions:  %10u\n", stats.priSeedAltContig);
  fprintf(
      statsFile,
      "  Liftover K-mer matching:   %10u  (%4.1f%%)\n",
      stats.priSeedAltLiftPresent,
      (double)stats.priSeedAltLiftPresent / NZERO(stats.priSeedAltContig) * 100);
  fprintf(
      statsFile,
      "  Liftover K-mer different:  %10u  (%4.1f%%)\n",
      stats.priSeedAltLiftAbsent,
      (double)stats.priSeedAltLiftAbsent / NZERO(stats.priSeedAltContig) * 100);
  fprintf(
      statsFile,
      "  No liftover:               %10u  (%4.1f%%)\n",
      stats.priSeedAltNoLiftover,
      (double)stats.priSeedAltNoLiftover / NZERO(stats.priSeedAltContig) * 100);
  fprintf(statsFile, "\n");
  fprintf(statsFile, "Raw primary-seed liftover groups: %10u\n", stats.rawLiftoverGroupCount);
  fprintf(
      statsFile,
      "  Average liftover group size:    %10.2f\n",
      (double)stats.rawLiftoverGroupSizeSum / NZERO(stats.rawLiftoverGroupCount));
  fprintf(statsFile, "  Histogram of liftover group sizes:\n");
  printHistogram(statsFile, stats.rawLiftoverGroupSizeHist, 64, 4, 1);
  fprintf(statsFile, "  Histogram of ALT hit count with no liftover:\n");
  printHistogram(statsFile, stats.rawNoLiftoverCountHist, 64, 4, 1);
  fprintf(statsFile, "\n");
  fprintf(statsFile, "Liftover groups after possible seed extension: %10u\n", stats.popLiftoverGroupCount);
  fprintf(
      statsFile,
      "  Liftover seed matching:     %10u  (%4.1f%%)\n",
      stats.liftoverPriHitsMatching,
      (double)stats.liftoverPriHitsMatching / NZERO(stats.popLiftoverGroupCount) * 100);
  fprintf(
      statsFile,
      "  Liftover seed injected:     %10u  (%4.1f%%)\n",
      stats.liftoverPriHitsAdded,
      (double)stats.liftoverPriHitsAdded / NZERO(stats.popLiftoverGroupCount) * 100);
  fprintf(
      statsFile,
      "  No liftover position:       %10u  (%4.1f%%)\n",
      stats.liftoverDummyPriSeeds,
      (double)stats.liftoverDummyPriSeeds / NZERO(stats.popLiftoverGroupCount) * 100);
  fprintf(
      statsFile,
      "  Average liftover group size:%10.2f\n",
      (double)stats.popLiftoverGroupSizeSum / NZERO(stats.popLiftoverGroupCount));
  fprintf(statsFile, "  Histogram of liftover group sizes:\n");
  printHistogram(statsFile, stats.popLiftoverGroupSizeHist, 64, 4, 1);
  fprintf(statsFile, "  Histogram of ALT hit count with no liftover:\n");
  printHistogram(statsFile, stats.popNoLiftoverCountHist, 64, 4, 1);
  fprintf(statsFile, "\n");
  fprintf(statsFile, "Hash records:\n");
  fprintf(statsFile, "  Bytes per record:    %10u\n", HASH_RECORD_BYTES);
  fprintf(statsFile, "  Number of records:   %10llu\n", numRecords);
  fprintf(
      statsFile,
      "  Hit records:         %10u  (%4.1f%%)\n",
      stats.countHitRec,
      (double)stats.countHitRec / NZERO(numRecords) * 100);
  fprintf(
      statsFile,
      "  Extension records:   %10u  (%4.1f%%)\n",
      stats.countExtendRec,
      (double)stats.countExtendRec / NZERO(numRecords) * 100);
  fprintf(
      statsFile,
      "  Interval records:    %10u  (%4.1f%%)\n",
      stats.countIntervalRec,
      (double)stats.countIntervalRec / NZERO(numRecords) * 100);
  fprintf(
      statsFile,
      "  Chain records:       %10u  (%4.1f%%)\n",
      stats.countChainRec,
      (double)stats.countChainRec / NZERO(numRecords) * 100);
  fprintf(
      statsFile,
      "  Empty records:       %10u  (%4.1f%%)\n",
      stats.countEmptyRec,
      (double)stats.countEmptyRec / NZERO(numRecords) * 100);
  fprintf(
      statsFile, "  Raw K-mer occupancy: %5.1f%%\n", (double)stats.countRawKmer / NZERO(numRecords) * 100);
  fprintf(
      statsFile,
      "  Final occupancy:     %5.1f%%\n",
      (double)(numRecords - stats.countEmptyRec) / NZERO(numRecords) * 100);
  fprintf(statsFile, "\n");
  fprintf(statsFile, "Hash buckets:\n");
  fprintf(statsFile, "  Records per bucket: %u\n", HASH_RECORDS_PER_BUCKET);
  fprintf(statsFile, "  Number of buckets:  %u\n", numBuckets);
  fprintf(statsFile, "  Histogram of raw K-mer bucket occupancy:\n");
  printHistogram(statsFile, stats.bucketLevelRawHist, 64, 4, 1);
  fprintf(statsFile, "  Histogram of bucket occupancy after extending or rejecting high frequency seeds:\n");
  printHistogram(statsFile, stats.bucketLevelEscHist, 64, 4, 1);
  fprintf(statsFile, "  Histogram of physical bucket occupancy as mapped:\n");
  printHistogram(statsFile, stats.bucketLevelMapHist, 64, 4, 0);
  fprintf(statsFile, "\n");
  fprintf(statsFile, "Seed extensions:\n");
  fprintf(statsFile, "  Base seed length:              %u\n", seedLen);
  fprintf(
      statsFile,
      "  Average extended seed length:  %.1f\n",
      (double)stats.extendLenWtdSum / NZERO(stats.countExtendKmer));
  fprintf(
      statsFile,
      "  Average extension increment:   %.1f\n",
      (double)stats.extendIncrWtdSum / NZERO(stats.countKmerExtendIncr));
  fprintf(
      statsFile,
      "  Average extension steps:       %.2f\n",
      (double)stats.extendStepsWtdSum / NZERO(stats.countExtendKmer));
  fprintf(statsFile, "  Extension IDs utilization:     %d%%\n", extIdBinPctMax);
  fprintf(statsFile, "  Portion of reference K-mers...\n");
  fprintf(statsFile, "    All raw K-mers:            %10u  (%5.1f%%)\n", stats.countRawKmer, 100.0);
  fprintf(
      statsFile,
      "    Extended to longer seeds:  %10u  (%5.1f%%)\n",
      stats.countExtendKmer,
      (double)stats.countExtendKmer / NZERO(stats.countRawKmer) * 100);
  fprintf(
      statsFile,
      "    Remaining as primary hit:  %10u  (%5.1f%%)\n",
      stats.countSeedHits[0],
      (double)stats.countSeedHits[0] / NZERO(stats.countRawKmer) * 100);
  fprintf(
      statsFile,
      "  Space in extension table:    %10u  (%5.1f%% of unmasked K-mers)\n",
      totExtTabRecs,
      (double)totExtTabRecs / NZERO(countUnmasked) * 100);
  fprintf(statsFile, "  Average frequencies of reference K-mers...\n");
  fprintf(
      statsFile,
      "    All raw K-mers:              %8.2f\n",
      (double)stats.kmerFreqSumSquare / NZERO(stats.countRawKmer));
  fprintf(
      statsFile,
      "    Extended to longer seeds:    %8.2f\n",
      (double)stats.kmerExtendFreqSum / NZERO(stats.countExtendKmer));
  fprintf(
      statsFile,
      "    Remaining as primary hit:    %8.2f\n",
      (double)stats.kmerHitFreqSum / NZERO(stats.countSeedHits[0]));
  fprintf(
      statsFile,
      "    As extended seed hit:        %8.2f\n",
      (double)stats.secKmerHitFreqSum / NZERO(stats.countSeedHits[1]));
  fprintf(
      statsFile,
      "    As primary or extended seed: %8.2f\n",
      (double)(stats.kmerHitFreqSum + stats.secKmerHitFreqSum) /
          NZERO((stats.countSeedHits[0] + stats.countSeedHits[1])));
  fprintf(statsFile, "  Extended seed length histogram:\n");
  printHistogram(statsFile, stats.extendLengthHist, MAX_EXTENDED_LENGTH + 1, 4, 0);
  fprintf(statsFile, "  Seed extension increment histogram:\n");
  printHistogram(statsFile, stats.extendIncrHist, 64, 4, 0);
  fprintf(statsFile, "  Seed extension steps histogram:\n");
  printHistogram(statsFile, stats.extendStepsHist, 64, 4, 0);
  fprintf(statsFile, "  Pre-extended K-mer frequency histogram:\n");
  printHistogram(statsFile, stats.extendFreqHist, 64, 4, 1);
  fprintf(statsFile, "  Remaining primary hit K-mer frequency histogram:\n");
  printHistogram(statsFile, stats.hitFreqHist[0], 64, 4, 1);
  fprintf(statsFile, "  Post-extended K-mer frequency histogram:\n");
  printHistogram(statsFile, stats.hitFreqHist[1], 64, 4, 1);
  fprintf(statsFile, "\n");
  fprintf(statsFile, "Hash chaining and probing:\n");
  fprintf(statsFile, "  Number of chains: %u\n", stats.countChains);
  fprintf(statsFile, "  Chain buckets:    %u\n", stats.countChainBuckets);
  fprintf(statsFile, "  Average length beyond each bucket...\n");
  fprintf(statsFile, "    chain:  %6.4f\n", (double)stats.chainLenSum / NZERO(numBuckets));
  fprintf(statsFile, "    probe:  %6.4f\n", (double)stats.probeLenSum / NZERO(numBuckets));
  fprintf(
      statsFile, "    either: %6.4f\n", (double)(stats.chainLenSum + stats.probeLenSum) / NZERO(numBuckets));
  fprintf(statsFile, "  Histogram of bucket probe lengths replaced by chaining:\n");
  printHistogram(statsFile, stats.probeDistChainedHist, MAX_PROBES + 1, 4, 1);
  fprintf(statsFile, "  Bucket chain length histogram:\n");
  printHistogram(statsFile, stats.chainLengthHist, 64, 4, 1);
  fprintf(statsFile, "  Bucket probe length histogram:\n");
  printHistogram(statsFile, stats.probeLengthHist, 64, 4, 0);
  fprintf(statsFile, "  Chain or probe length histogram:\n");
  printHistogram(statsFile, stats.chainOrProbeLengthHist, 64, 4, 1);
  fprintf(statsFile, "\n");
  fprintf(statsFile, "Compression:          Records        Bits    Mean\n");
  fprintf(
      statsFile,
      "  auto pri hits:   %10lld  %10lld  %6.3f\n",
      compCtx.autoPriRecs,
      compCtx.autoPriBits,
      (double)compCtx.autoPriBits / NZERO(compCtx.autoPriRecs));
  fprintf(
      statsFile,
      "  auto sec hits:   %10lld  %10lld  %6.3f\n",
      compCtx.autoSecRecs,
      compCtx.autoSecBits,
      (double)compCtx.autoSecBits / NZERO(compCtx.autoSecRecs));
  fprintf(
      statsFile,
      "  auto nul hits:   %10lld  %10lld  %6.3f\n",
      compCtx.autoNulRecs,
      compCtx.autoNulBits,
      (double)compCtx.autoNulBits / NZERO(compCtx.autoNulRecs));
  fprintf(
      statsFile,
      "  special hits:    %10lld  %10lld  %6.3f\n",
      compCtx.specialRecs,
      compCtx.specialBits,
      (double)compCtx.specialBits / NZERO(compCtx.specialRecs));
  fprintf(
      statsFile,
      "  chain pointers:  %10lld  %10lld  %6.3f\n",
      compCtx.chainPtrRecs,
      compCtx.chainPtrBits,
      (double)compCtx.chainPtrBits / NZERO(compCtx.chainPtrRecs));
  fprintf(
      statsFile,
      "  chain ends:      %10lld  %10lld  %6.3f\n",
      compCtx.chainEndRecs,
      compCtx.chainEndBits,
      (double)compCtx.chainEndBits / NZERO(compCtx.chainEndRecs));
  fprintf(
      statsFile,
      "  literals:        %10lld  %10lld  %6.3f\n",
      compCtx.literalRecs,
      compCtx.literalBits,
      (double)compCtx.literalBits / NZERO(compCtx.literalRecs));
  fprintf(
      statsFile,
      "  ext literals:    %10lld  %10lld  %6.3f\n",
      compCtx.extLitRecs,
      compCtx.extLitBits,
      (double)compCtx.extLitBits / NZERO(compCtx.extLitRecs));
  fprintf(
      statsFile,
      "  TOTAL:           %10lld  %10lld  %6.3f\n",
      compCtx.totalRecs,
      compCtx.totalBits,
      (double)compCtx.totalBits / NZERO(compCtx.totalRecs));
  fprintf(statsFile, "  Misc bits:   %lld\n", compCtx.miscBits);
  fprintf(statsFile, "  Final bits:  %lld\n", compCtx.bitLen);
  fprintf(statsFile, "  Final bytes: %lld\n", compCtx.bitLen / 8);
  fprintf(statsFile, "\n");
  if (stats.cyclesOverhead) {
    fprintf(statsFile, "Build thread cycle counts:\n");
    fprintf(statsFile, "  cyclesOverhead:        %lld\n", stats.cyclesOverhead);
    fprintf(statsFile, "  cyclesBucketOverhead:  %lld\n", stats.cyclesBucketOverhead);
    fprintf(statsFile, "  cyclesExtendPrep:      %lld\n", stats.cyclesExtendPrep);
    fprintf(statsFile, "  cyclesExtendSort:      %lld\n", stats.cyclesExtendSort);
    fprintf(statsFile, "  cyclesLiftover:        %lld\n", stats.cyclesLiftover);
    fprintf(statsFile, "  cyclesPriSort:         %lld\n", stats.cyclesPriSort);
    fprintf(statsFile, "  cyclesExtendFreq:      %lld\n", stats.cyclesExtendFreq);
    fprintf(statsFile, "  cyclesExtendDynProg:   %lld\n", stats.cyclesExtendDynProg);
    fprintf(statsFile, "  cyclesExtendIntervals: %lld\n", stats.cyclesExtendIntervals);
    fprintf(statsFile, "  cyclesExtendConstruct: %lld\n", stats.cyclesExtendConstruct);
    fprintf(statsFile, "  cyclesBucketSort:      %lld\n", stats.cyclesBucketSort);
    fprintf(statsFile, "  cyclesBucketOrganize:  %lld\n", stats.cyclesBucketOrganize);
    fprintf(statsFile, "  cyclesBucketChain:     %lld\n", stats.cyclesBucketChain);
    fprintf(statsFile, "  cyclesBucketWrite:     %lld\n", stats.cyclesBucketWrite);
    fprintf(statsFile, "  cyclesBucketCompress:  %lld\n", stats.cyclesBucketCompress);
    fprintf(statsFile, "\n");
  }
  fclose(statsFile);
  printf("Wrote statistics to '%s'\n", statsFname);

mainBuildHashTableCleanup:

  config->hdr->extTabRecs   = totExtTabRecs;
  config->hdr->extTabDigest = extTabDigest;
  config->hdr->hashDigest   = hashDigest;
  free(crcInit);
  if (writeHashFile) {
    fclose(hashFile);
    fclose(extTabFile);
  }
  writeCompHashTableClose(&compCtx);
  for (i = 0; i < hashChunks; i++) {
    for (n = 0; n < chunkBuckets; n++)
      if (bucketAlloc[i][n] > HASH_RECORDS_PER_BUCKET) free(bucket[i][n]);
    free(bucketAlloc[i]);
    free(bucketCount[i]);
    free(bucketLocks[i]);
    free(chunkLitFlags[i]);
    free(bucket[i]);
    free(physRecords[i]);
  }
  free(seedPopRecs);
  free(extIndexRecs);

#if defined(_TARGET_PPC_)
  // Destroy the mutex
  //
  pthread_mutex_destroy(&lock);
#endif
  if (errMsg[0])
    return errMsg;
  else
    return NULL;
}
