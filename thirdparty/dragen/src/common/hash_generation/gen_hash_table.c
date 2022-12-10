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
//

#define __STDC_FORMAT_MACROS  // for PRIu64
#include <assert.h>
#include <ctype.h>
#include <inttypes.h>
#include <libgen.h>
#include <math.h>
#include <pthread.h>
#include <semaphore.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>

#include "crc_hash.h"
#include "gen_hash_table.h"
#include "hash_cfg_file.h"
#include "hash_table.h"
#include "host_version_hashtable.h"
#include "liftover.h"
#include "methylation_hash_table.h"
#ifndef LOCAL_BUILD
#include "alt_contig_tracker.h"
#include "crc32_hw.h"
#include "reference_names.h"
#include "syslogger.h"
#include "watchdog.h"
#ifdef _TARGET_PPC_
#include "crc32_powerpc.h"
#endif

/* ============ end of #include directives ============ */
#endif

// The following disable these specific compiler warnings, which are numerous in
// this file
#pragma GCC diagnostic ignored "-Wtype-limits"
#pragma GCC diagnostic ignored "-Wformat"
#pragma GCC diagnostic ignored "-Wchar-subscripts"
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"

#ifndef LOCAL_BUILD
extern uint32_t crc32c_hw(uint32_t crc, const void* buf, size_t len);
#endif

#define MAX_LINE (REF_SEQ_END_PAD_BASES + 4096)
#define MAX_LIFT_LINE (1024 * 1024)
#define SEQ_REPORT_INTERVAL 500000000

#ifndef LOCAL_BUILD
static unsigned long watchdogID = ULONG_MAX;
#endif

static char ERR_MSG[16384] = "";

static char defHashFname[]      = "hash_table.bin";
static char defCompFname[]      = "hash_table.cmp";
static char defExtTabFname[]    = "extend_table.bin";
static char defRefFname[]       = "reference.bin";
static char defRefIdxFname[]    = "ref_index.bin";
static char defRepMaskFname[]   = "repeat_mask.bin";
static char defStatsFname[]     = "hash_table_stats.txt";
static char defConfigFname[]    = "hash_table.cfg";
static char defConfigBinFname[] = "hash_table.cfg.bin";
static char defStrFname[]       = "str_table.bin";
static char defPopSnpsFname[]   = "ref_pop_snps.bin";

// When building hashtables for bisulfite sequenced methylation analysis, need to
// append strings to each contig's sequence name.
static const char seqnameSuffixDefault[] = "";
const char        seqnameSuffixCtoT[]    = "_CT_converted";
const char        seqnameSuffixGtoA[]    = "_GA_converted";

// Short Tandem Repeat (STR) parameters
#define STR_MAX_PERIOD 8
#define STR_MAX_REPEATLEN 20
// Initial log2 decimation factors, selected to keep human exome bins over 1K samples,
// when possible, and with short repeats sampled extra because their absolute indel
// rates are low.  Further decimation will occur as aligned reads are processed,
// but initial decimation reduces STRs output during hash table construction
// from 1.5 billion to 6.5 million for human genome.
static const int strDecimationLog2[STR_MAX_PERIOD + 1][STR_MAX_REPEATLEN + 1] = {
    // RepeatLen 0  1  2  3  4  5  6  7  8  9 ...
    {0},                               // Period 0
    {0, 10, 10, 9, 8, 7, 5, 3, 1, 0},  // Period 1
    {0, 0, 9, 6, 3, 0},                // Period 2
    {0, 0, 8, 4, 1, 0},                // Period 3
    {0, 0, 6, 0},                      // Period 4
    {0, 0, 5, 0},                      // Period 5
    {0, 0, 4, 0},                      // Period 6
    {0, 0, 1},                         // Period 7
    {0}};                              // Period 8
// Base character lookup ([Aa],[Cc],[Gg],[Tt]) => (1,2,3,4)
static const char STR_BASE[256] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 2, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 2, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

#define APPLY_MIN(val, min) ((val) = ((val) < (min) ? (min) : (val)))
#define APPLY_MAX(val, max) ((val) = ((val) > (max) ? (max) : (val)))
#define APPLY_MINMAX(val, min, max) ((val) = ((val) < (min) ? (min) : (val) > (max) ? (max) : (val)))

#define DEFAULT_GB 8  // Default maximum 1GB table chunks

// Expected chromosome 1-22,X,Y lengths for GRCh37-38
#define GRCH_VERSION_LO 37
#define GRCH_VERSION_HI 38
#define GRCH_VERSIONS (GRCH_VERSION_HI - GRCH_VERSION_LO + 1)
static const uint32_t grchChromLens[GRCH_VERSIONS][24] = {
    {249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022,
     141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753,
     81195210,  78077248,  59128983,  63025520,  48129895,  51304566,  155270560, 59373566},
    {248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636,
     138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345,
     83257441,  80373285,  58617616,  64444167,  46709983,  50818468,  156040895, 57227415}};

//-------------------------------------------------------------------------------swhitmore
// setDefaultHashParams - Sets the default parameters shared between the host software
// and build_hash_table.
//
void setDefaultHashParams(hashTableConfig_t* defConfig, const char* destDir, HashTableType hashTableType)
{
  defConfig->hdr = (hashTableHeader_t*)malloc(sizeof(hashTableHeader_t));
  memset(defConfig->hdr, 0, sizeof(hashTableHeader_t));
  defConfig->hdr->hashTableVersion = HT_VERSION;

  size_t dirLen = 0;
  if (destDir) {
    dirLen = strlen(destDir) + 1;  // Plus one for ending '/'
  }

  // Determine the type of hash table being generated. If not the normal one, then
  // automatically append reference subdirectory names.
  char* dir = 0;
  if (hashTableType == HT_TYPE_NORMAL) {
    if (destDir) dir = strdup(destDir);
  } else {
    const char* SUBDIR_METHYL_C_TO_T   = "CT_converted";
    const char* SUBDIR_METHYL_G_TO_A   = "GA_converted";
    const char* SUBDIR_METHYL_COMBINED = "methyl_converted";
    const char* SUBDIR_ANCHORED        = "anchored_rna";
    const char* subDir                 = 0;

    switch (hashTableType) {
    case HT_TYPE_METHYL_C_TO_T:
      subDir = SUBDIR_METHYL_C_TO_T;
      break;
    case HT_TYPE_METHYL_G_TO_A:
      subDir = SUBDIR_METHYL_G_TO_A;
      break;
    case HT_TYPE_METHYL_COMBINED:
      subDir = SUBDIR_METHYL_COMBINED;
      break;
    case HT_TYPE_ANCHORED:
      subDir = SUBDIR_ANCHORED;
      break;
    default:
      assert(0 && "Unrecognized hash table type in building hashtable");
    }

    const size_t fullDirLen = dirLen + strlen(subDir) + 2;  // +1 for /, +1 for null terminator
    dir                     = (char*)malloc(fullDirLen);
    if (destDir) {
      snprintf(dir, fullDirLen, "%s/%s", destDir, subDir);
    } else {
      snprintf(dir, fullDirLen, "%s", subDir);
    }
    dirLen = fullDirLen;

    // Create the directory if it does not exist
    struct stat st = {0};
    if (stat(dir, &st) == -1) {
      mkdir(dir, S_IRWXU | S_IRWXG | S_IROTH);
    }
  }

  // Assemble full paths to the various files:
  size_t configLen    = strlen(defConfigFname) + 1;  // Plus one for end of string
  size_t configBinLen = strlen(defConfigBinFname) + 1;
  size_t hashLen      = strlen(defHashFname) + 1;
  size_t compLen      = strlen(defCompFname) + 1;
  size_t extTabLen    = strlen(defExtTabFname) + 1;
  size_t refLen       = strlen(defRefFname) + 1;
  size_t refIdxLen    = strlen(defRefIdxFname) + 1;
  size_t repMaskLen   = strlen(defRepMaskFname) + 1;
  size_t statsLen     = strlen(defStatsFname) + 1;
  size_t strLen       = strlen(defStrFname) + 1;

  defConfig->configFname    = (char*)malloc((dir ? (dirLen + configLen) : configLen));
  defConfig->configBinFname = (char*)malloc((dir ? (dirLen + configBinLen) : configBinLen));
  defConfig->hashFname      = (char*)malloc((dir ? (dirLen + hashLen) : hashLen));
  defConfig->compFname      = (char*)malloc((dir ? (dirLen + compLen) : compLen));
  defConfig->extTabFname    = (char*)malloc((dir ? (dirLen + extTabLen) : extTabLen));
  defConfig->refOutput      = (char*)malloc((dir ? (dirLen + refLen) : refLen));
  defConfig->refIdxFname    = (char*)malloc((dir ? (dirLen + refIdxLen) : refIdxLen));
  defConfig->repMaskFname   = (char*)malloc((dir ? (dirLen + repMaskLen) : repMaskLen));
  defConfig->statsFname     = (char*)malloc((dir ? (dirLen + statsLen) : statsLen));
  defConfig->strFname       = (char*)malloc((dir ? (dirLen + strLen) : strLen));

  if (dir) {
    snprintf(defConfig->configFname, dirLen + configLen, "%s/%s", dir, defConfigFname);
    snprintf(defConfig->configBinFname, dirLen + configBinLen, "%s/%s", dir, defConfigBinFname);
    snprintf(defConfig->hashFname, dirLen + hashLen, "%s/%s", dir, defHashFname);
    snprintf(defConfig->compFname, dirLen + compLen, "%s/%s", dir, defCompFname);
    snprintf(defConfig->extTabFname, dirLen + extTabLen, "%s/%s", dir, defExtTabFname);
    snprintf(defConfig->refOutput, dirLen + refLen, "%s/%s", dir, defRefFname);
    snprintf(defConfig->refIdxFname, dirLen + refIdxLen, "%s/%s", dir, defRefIdxFname);
    snprintf(defConfig->repMaskFname, dirLen + repMaskLen, "%s/%s", dir, defRepMaskFname);
    snprintf(defConfig->statsFname, dirLen + statsLen, "%s/%s", dir, defStatsFname);
    snprintf(defConfig->strFname, dirLen + strLen, "%s/%s", dir, defStrFname);

  } else {
    strcpy(defConfig->configFname, defConfigFname);
    strcpy(defConfig->configBinFname, defConfigBinFname);
    strcpy(defConfig->hashFname, defHashFname);
    strcpy(defConfig->compFname, defCompFname);
    strcpy(defConfig->extTabFname, defExtTabFname);
    strcpy(defConfig->refOutput, defRefFname);
    strcpy(defConfig->refIdxFname, defRefIdxFname);
    strcpy(defConfig->repMaskFname, defRepMaskFname);
    strcpy(defConfig->statsFname, defStatsFname);
    strcpy(defConfig->strFname, defStrFname);

  }

  if (dir) {
    free(dir);
  }

  defConfig->hostVersion = (char*)getHostVersion();
}

//-------------------------------------------------------------------------------swhitmore
// freeHashParams - Frees allocated buffers inside a hash tables config record.
//
void freeHashParams(hashTableConfig_t* config)
{
  if (config->usedReadBuf) {
    free(config->readBuf);
  } else {
    size_t i;
    for (i = 0; i < config->hdr->numRefSeqs; ++i) {
      free(config->refSeq[i]);
      free(config->seqName[i]);
    }
    free(config->hdr);
    free(config->configFname);
    free(config->configBinFname);
    free(config->hashFname);
    free(config->compFname);
    free(config->extTabFname);
    free(config->refOutput);
    free(config->refIdxFname);
    free(config->repMaskFname);
    free(config->statsFname);
    free(config->cmdLine);
    free(config->strFname);
    free(config->popSnpsOutput);
  }
  free(config->refSeq);
  free(config->seqName);
}

//-------------------------------------------------------------------------miker/swhitmore
// decodeSizeArg - Parse hash table/memory limit size strings
//
int decodeSizeArg(const char* s, uint64_t* bytes)
{
  if (!s) {
    *bytes = 0;
    return 1;
  }
  char*  tail;
  double fbytes = strtof(s, &tail);
  switch (toupper(*tail)) {
  case '\0':
    break;
  case 'B':
    break;
  case 'K':
    fbytes *= (1 << 10);
    break;
  case 'M':
    fbytes *= (1 << 20);
    break;
  case 'G':
    fbytes *= (1 << 30);
    break;
  default:
    snprintf(
        ERR_MSG,
        sizeof(ERR_MSG),
        "Invalid size units on argument '%s'.  May use B (or blank), KB, MB, GB.\n",
        s);
    return 0;
  }
  *bytes = fbytes;
  if (*bytes == 0) {
    return 1;
  }

  if (*bytes != fbytes) {
    snprintf(ERR_MSG, sizeof(ERR_MSG), "Invalid size units on argument '%s' - fractional bytes.\n", s);
    return 0;
  }
  return 1;
}

//-------------------------------------------------------------------------------swhitmore
// setCmdLine - Update config with command line used to build hash table.
//
void setCmdLine(hashTableConfig_t* config, int argc, char* argv[])
{
  size_t cmdLineLen = 0;
  int    i;
  // Calculate length of string, which will include NULL termination charater.
  for (i = 0; i < argc; ++i) {
    cmdLineLen += strlen(argv[i]) + 1;
  }
  config->cmdLine = (char*)malloc(cmdLineLen);
  memset(config->cmdLine, 0, cmdLineLen);
  for (i = 0; i < argc; ++i) {
    strcat(config->cmdLine, argv[i]);
    if (i < (argc - 1)) {
      strcat(config->cmdLine, " ");
    } else {
      strcat(config->cmdLine, "\0");
    }
  }
}

//-------------------------------------------------------------------------------swhitmore
// printInternalParams - print all parameters not already printed to stdout (for testing).
//
void printInternalParams(hashTableConfig_t* config)
{
  printf("\nInternal Settings:\n");
  printf("  Maximum threads                    : %d\n", config->maxThreads);
  printf("  Maximum GB table chunks            : %d\n", config->maxGB);
  printf("  Write uncompressed hash table      : %d\n", config->writeHashFile);
  printf("  Write compressed hash table        : %d\n", config->writeCompFile);
  printf("  Size of hash table                 : %s\n", config->sizeStr);
  printf("  Memory limit                       : %s\n", config->memSizeStr);
  printf("  SJ reserved space                  : %s\n", config->sjSizeStr);
  printf("  Methylated DNA conversion          : %u\n", config->methylatedConv);
  printf("  Primary polynomial index           : %d\n", config->priPolyIndex);
  printf("  Extended polynomial index          : %d\n", config->secPolyIndex);
  printf("  Name of FASTA input file           : %s\n", config->refInput);
  printf("  Name of ALT contig liftover file   : %s\n", config->altLiftover);
  printf("  Hash table configuration file (txt): %s\n", config->configFname);
  printf("  Hash table configuration file (bin): %s\n", config->configBinFname);
  printf("  Hash table output file             : %s\n", config->hashFname);
  printf("  Extension table output file        : %s\n", config->extTabFname);
  printf("  Reference output file              : %s\n", config->refOutput);
  printf("  Stats output file                  : %s\n", config->statsFname);
  if (config->popSnpsInput && config->popSnpsOutput) {
    printf("  Population SNPs output file        : %s\n", config->popSnpsOutput);
  }
  if (config->maskBed) {
    printf("  Mask bed file                      : %s\n", config->maskBed);
  }
  printf("  Host Version                       : %s\n", config->hostVersion);
  printf("  Override hash table check          : %d\n", config->overrideCheck);
  printf("  Test only                          : %d\n", config->testOnly);
  printf("\n");
}

void freeRefRec(refRec_t* p)
{
  free(p->name);
  free(p->seq);
  free(p->str);
  free(p->lifted);
  free(p->liftMatch);
  p->name = p->seq = NULL;
  p->str           = NULL;
  p->lifted        = NULL;
  p->liftMatch     = NULL;
}

typedef struct {
  refRec_t*          refRecs;
  int32_t            numRefSeqs;
  int32_t*           nextRefSeq;
  pthread_mutex_t*   mut;
  uint8_t*           seedHashCounts;
  uint8_t*           seedHashLocks;
  uint64_t           validSeeds;
  uint64_t           extendedSeeds;
  hashTableHeader_t* cfghdr;
} strScanCtx_t;

// Thread function to scan reference sequences for short tandem repeats (STRs)
void* strScanThreadMultBase(void* ctxPtr)
{
  strScanCtx_t* ctx = (strScanCtx_t*)ctxPtr;
  refRec_t*     r   = NULL;
  int32_t       per, len, refSeq;
  // Iniitialize decimation masks
  uint32_t masks[STR_MAX_PERIOD + 1][STR_MAX_REPEATLEN + 1] = {{0}};
  for (per = 0; per <= STR_MAX_PERIOD; per++)
    for (len = 0; len <= STR_MAX_REPEATLEN; len++) masks[per][len] = (1 << strDecimationLog2[per][len]) - 1;
  while (1) {
    // Grab the next reference sequence, if any
    pthread_mutex_lock(ctx->mut);
    r = *ctx->nextRefSeq < ctx->numRefSeqs ? &ctx->refRecs[refSeq = (*ctx->nextRefSeq)++] : NULL;
    pthread_mutex_unlock(ctx->mut);
    if (!r) return NULL;
    int32_t        pos = 0, prevBeg = -1, beg, end;
    int32_t        bestLen = 0, bestPer = 0, bestBeg = 0, bestEnd = 0, bestRun;
    unsigned char* seq       = (unsigned char*)r->seq;
    uint8_t*       lifted    = r->lifted;
    uint8_t*       liftMatch = r->liftMatch;
    uint32_t       cnt, count[STR_MAX_PERIOD + 1][STR_MAX_REPEATLEN + 1] = {{0}};
    uint32_t       begTrim = r->begTrim;
    uint32_t       trimLen = r->trimLen;
    htStrRec_t     str     = {0};
    str.refId              = r->seqNum;
    int num                = 0;
    int isAlt              = r->isAlt;
    int isLifted, liftMatches, incBy, extIncBy;
    // Pre-hash seed k-mers to count extended seeds and size the seed extension table
    uint8_t* seedHashCounts = ctx->seedHashCounts;
    uint8_t* seedHashLocks  = ctx->seedHashLocks;
    uint8_t* lockPtr;
    uint8_t* countPtr;
    uint32_t seedLength       = ctx->cfghdr->priSeedBases;
    uint32_t anchorBinBits    = ctx->cfghdr->anchorBinBits ? ctx->cfghdr->anchorBinBits : 31;
    double   refSeedInterval  = ctx->cfghdr->refSeedInterval;
    uint32_t liftMatchSeedInt = ctx->cfghdr->liftMatchSeedInt ? ctx->cfghdr->liftMatchSeedInt : 1;
    int32_t  bp = 0, nextPos = 0, nextIndex = 1;
    uint64_t base, comp, fwBaseSeed = 0, rcBaseSeed = 0, seed, seedMask = (1ULL << (2 * seedLength)) - 1,
                         hash;
    uint64_t      anchorBin;
    uint64_t      bitMask          = (1ULL << seedLength) - 1;
    uint32_t      anchorRefSeq     = ctx->cfghdr->anchorBinBits ? refSeq : 0;
    uint64_t*     seedPtr          = &seed;
    int           minFreqToExtend  = ctx->cfghdr->minFreqToExtend;
    int           anchorShift      = seedLength * 2;
    int           compShift        = (seedLength - 1) * 2;
    char*         charToCode       = r->charToCode;
    uint32_t      maxMultBaseSeeds = ctx->cfghdr->maxMultBaseSeeds;
    unsigned char c;
    int           numb, zero, code, numZeros = seedLength, i, j, palindrome;
    uint8_t       multNum[32], multIdx[32], multNumPos;
    uint8_t       numHist[32]   = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    uint8_t       zerosHist[32] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                             1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    uint8_t       codeHist[32]  = {0};
    uint64_t      numKmers      = 1, fwSeed, rcSeed, multFwSeeds[32][4], multRcSeeds[32][4];
    for (bp = 0; bp < trimLen; bp++) {
      // Convert IUB character to 4-bit code (mask of matching bases),
      // and count of matching bases (minimum of 1, separate flag for zero)
      c    = seq[bp];
      code = charToCode[c];
      numb = baseCodeNumBasesMin1[code];
      zero = (code == 0);
      // Keep a history of count and zero flag
      codeHist[bp & 31]  = code;
      numHist[bp & 31]   = numb;
      zerosHist[bp & 31] = zero;
      // Update number of zeros and product of counts within the seed length
      numZeros   = numZeros + zero - zerosHist[(bp - seedLength) & 31];
      numKmers   = numKmers * numb / numHist[(bp - seedLength) & 31];
      anchorBin  = bp >> anchorBinBits;
      base       = baseCodeBase[code];
      comp       = base ^ 3;
      fwBaseSeed = ((fwBaseSeed << 2) & seedMask) | base;
      rcBaseSeed = (rcBaseSeed >> 2) | (comp << compShift);
      // Positions where seeds should get populated
      // (This is onlt approximate / statistical, because here we restart at
      // position zero in each chromosome, rather than proper global coordinates.)
      if (bp == nextPos) {
        nextPos = floor(nextIndex++ * refSeedInterval);
        // Not if overlapping an invalid base, or too many multi-base codes
        if (!numZeros && numKmers <= maxMultBaseSeeds) {
          // For alt contig positions where the whole seed matches the liftover sequence, skip most
          // positions, as we only populate one per liftMatchSeedInt seed indexes
          liftMatches =
              isAlt && liftMatch && ((*((uint64_t*)&liftMatch[bp >> 3]) >> (bp & 7)) & bitMask) == bitMask;
          if (liftMatches && (nextIndex % liftMatchSeedInt)) continue;
          // Decide when to count an extra seed to populate, because nonmatching liftover hits
          // may get injected along with alt contig seeds.  Always count double for seeds with
          // liftover at nonmatching positions.  For seeds matching their liftover, count once
          // until the tally is high enough for seed extension, then count twice because some
          // of the extended seeds will stop matching the liftover position.
          isLifted = isAlt && (!lifted || ((*((uint64_t*)&lifted[bp >> 3]) >> (bp & 7)) & bitMask));
          incBy    = 1 + (isLifted && !liftMatches);
          extIncBy = 1 + isLifted;
          // Prepare iteration over choices for multi-base codes
          multNumPos = 0;
          if (numKmers > 1) {
            for (i = 0; i < seedLength; i++) {
              if ((numb = numHist[(bp - i) & 31]) > 1) {
                code                = codeHist[(bp - i) & 31];
                multNum[multNumPos] = numb;
                multIdx[multNumPos] = 0;
                // Erase this position from the base seeds
                fwBaseSeed &= ~((uint64_t)3 << (i << 1));
                rcBaseSeed &= ~((uint64_t)3 << ((seedLength - 1 - i) << 1));
                for (j = 0; j < numb; j++) {
                  multFwSeeds[multNumPos][j] = (uint64_t)baseCodeBaseList[code][j] << (i << 1);
                  multRcSeeds[multNumPos][j] = ((uint64_t)baseCodeBaseList[code][j] ^ 3)
                                               << ((seedLength - 1 - i) << 1);
                }
                multNumPos++;
              }
            }
          }
          // Loop through possibly multiple combinations of choices for multi-base codes
          while (1) {
            ctx->validSeeds++;
            // Grab the baseline k-mer, and its reverse complement
            fwSeed = fwBaseSeed;
            rcSeed = rcBaseSeed;
            // Inject current choices for multi-base positions
            for (i = 0; i < multNumPos; i++) {
              fwSeed |= multFwSeeds[i][multIdx[i]];
              rcSeed |= multRcSeeds[i][multIdx[i]];
            }
            seed       = (fwSeed < rcSeed ? fwSeed : rcSeed) | (anchorBin << anchorShift);
            palindrome = (fwSeed == rcSeed);
            // CRC32c hash the seed
#if defined(LOCAL_BUILD) && defined(__x86_64__)
            __asm__ __volatile__(
                "crc32q\t"
                "(%1), %0"
                : "=r"(hash)
                : "r"((uint64_t*)seedPtr), "0"(0));
#elif !defined(LOCAL_BUILD)
            hash = crc32c_hw(0, (const unsigned char*)seedPtr, 8);
#endif
            // Add reference sequence ID in anchored mode, to separate difference sequences
            hash += anchorRefSeq;
            // Obtain a lock corresponding to a 20-bit segment of the hash
            lockPtr  = &seedHashLocks[hash & 0xFFFFF];
            countPtr = &seedHashCounts[hash & 0xFFFFFFFF];
            while (!__sync_bool_compare_and_swap(lockPtr, 0, 1))
              ;
            // Count palindromes twice
            for (; palindrome >= 0; palindrome--) {
              // Increment count for this hash
              if (*countPtr < minFreqToExtend) {
                // When count reaches the extension threshold, count them all as extended
                if ((*countPtr += incBy) >= minFreqToExtend) {
                  ctx->extendedSeeds += *countPtr + extIncBy - incBy;
                }
              }
              // And count further seeds with this hash as extended
              else {
                ctx->extendedSeeds += extIncBy;
              }
            }
            // Release the lock
#if defined(_TARGET_PPC_)
            __asm__ __volatile__("sync" ::: "memory");
#endif
            *lockPtr = 0;
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
      }
    }

    // skip decoy and population alt contigs
    if (r->isDecoy || r->isPopAlt) continue;

    // Allocate STR records in proportion to sequence length, assuming decimation
    int alloc = (trimLen >> (strDecimationLog2[1][1] - 2)) + 256;
    if (!(r->str = calloc(alloc, sizeof(htStrRec_t)))) continue;
    // Scan the reference sequence
    while (pos < trimLen) {
      bestLen = 0;
      bestEnd = pos + 1;
      // Find the best STR, with greatest repeatLen, breaking ties to lower period
      for (per = 1; per <= STR_MAX_PERIOD; per++) {
        // Skip unless [ACGT] through first whole period
        if (!STR_BASE[seq[pos + per - 1]]) break;
        // Leftward extension of repeat with this period
        for (beg = pos - 1; beg > prevBeg; beg--) {
          if (r->charToBase[seq[beg]] != r->charToBase[seq[beg + per]]) break;
        }
        beg++;
        // Rightward extension of repeat with this period
        for (end = pos + per; end < trimLen; end++) {
          if (r->charToBase[seq[end]] != r->charToBase[seq[end - per]]) break;
        }
        // Track maximum number of whole periods
        len = (end - beg) / per;
        if (len > bestLen) {
          bestLen = len;
          bestPer = per;
          bestBeg = beg;
          bestEnd = end;
        }
      }
      // Process valid STR
      if (bestLen) {
        // Cap repeatLength
        if (bestLen > STR_MAX_REPEATLEN) bestLen = STR_MAX_REPEATLEN;
        bestRun = bestEnd - bestBeg;
        // Test incrementing count per bin against decimation mask
        // Adjust upward by reference ID to randomize downsampling positions, important for many short contigs
        cnt = str.refId + count[bestPer][bestLen]++;
        if (!(cnt & masks[bestPer][bestLen])) {
          // Enlarge allocation 25% if needed
          if (num >= alloc) {
            alloc += alloc >> 2;
            if (!(r->str = realloc(r->str, alloc * sizeof(htStrRec_t)))) break;
          }
          // Save STR record
          str.mask      = cnt;
          str.begPos    = bestBeg + begTrim;
          str.length    = bestRun > 0xFFFF ? 0xFFFF : bestRun;
          str.period    = bestPer;
          str.repeatLen = bestLen;
          r->str[num++] = str;
        }
      }
      pos = bestEnd;
    }
    // Final STR record count
    r->str    = realloc(r->str, (num ? num : 1) * sizeof(htStrRec_t));
    r->strNum = num;
  }
}

// Thread function to scan reference sequences for short tandem repeats (STRs)
// TODO If (maxMultBaseSeeds > 0), call the multi base version. Remove this version and make the multi base
// the only version of this function for the next HT version.
void* strScanThread(void* ctxPtr)
{
  strScanCtx_t* ctx = (strScanCtx_t*)ctxPtr;

  uint32_t maxMultBaseSeeds = ctx->cfghdr->maxMultBaseSeeds;
  if (maxMultBaseSeeds > 0) {
    return strScanThreadMultBase(ctxPtr);
  }

  refRec_t* r = NULL;
  int32_t   per, len, refSeq;
  // Iniitialize decimation masks
  uint32_t masks[STR_MAX_PERIOD + 1][STR_MAX_REPEATLEN + 1] = {{0}};
  for (per = 0; per <= STR_MAX_PERIOD; per++)
    for (len = 0; len <= STR_MAX_REPEATLEN; len++) masks[per][len] = (1 << strDecimationLog2[per][len]) - 1;
  while (1) {
    // Grab the next reference sequence, if any
    pthread_mutex_lock(ctx->mut);
    r = *ctx->nextRefSeq < ctx->numRefSeqs ? &ctx->refRecs[refSeq = (*ctx->nextRefSeq)++] : NULL;
    pthread_mutex_unlock(ctx->mut);
    if (!r) return NULL;
    int32_t        pos = 0, prevBeg = -1, beg, end;
    int32_t        bestLen = 0, bestPer = 0, bestBeg = 0, bestEnd = 0, bestRun;
    unsigned char* seq = (unsigned char*)r->seq;
    uint32_t       cnt, count[STR_MAX_PERIOD + 1][STR_MAX_REPEATLEN + 1] = {{0}};
    uint32_t       begTrim = r->begTrim;
    uint32_t       trimLen = r->trimLen;
    htStrRec_t     str     = {0};
    str.refId              = r->seqNum;
    int num                = 0;
    // Pre-hash seed k-mers to count extended seeds and size the seed extension table
    uint8_t*  seedHashCounts = ctx->seedHashCounts;
    uint8_t*  seedHashLocks  = ctx->seedHashLocks;
    uint8_t*  lockPtr;
    uint8_t*  countPtr;
    uint32_t  seedLength      = ctx->cfghdr->priSeedBases;
    uint32_t  anchorBinBits   = ctx->cfghdr->anchorBinBits ? ctx->cfghdr->anchorBinBits : 31;
    double    refSeedInterval = ctx->cfghdr->refSeedInterval;
    int32_t   bp = 0, nextPos = 0, nextIndex = 1;
    uint64_t  x, base, comp, seedFw = 0, seedRc = 0, seed, seedMask = (1ULL << (2 * seedLength)) - 1, hash;
    uint64_t  anchorBin;
    uint32_t  anchorRefSeq = ctx->cfghdr->anchorBinBits ? refSeq : 0;
    uint32_t  invl, seedIn = 0xFFFFFFFF, invlMask = (1 << seedLength) - 1;
    uint64_t* seedPtr         = &seed;
    int       minFreqToExtend = ctx->cfghdr->minFreqToExtend;
    int       anchorShift     = seedLength * 2;
    int       compShift       = (seedLength - 1) * 2;
    int*      charToBase      = r->charToBase;
    for (bp = 0; bp < trimLen; bp++) {
      anchorBin = bp >> anchorBinBits;
      x         = charToBase[seq[bp]];
      base      = x & 3;
      comp      = base ^ 3;
      invl      = x >> 2;
      seedFw    = ((seedFw << 2) & seedMask) | base;
      seedRc    = (seedRc >> 2) | (comp << compShift);
      seedIn    = ((seedIn << 1) & invlMask) | invl;
      seed      = (seedFw < seedRc ? seedFw : seedRc) | (anchorBin << anchorShift);
      if (bp == nextPos) {
        nextPos = floor(nextIndex++ * refSeedInterval);
        // Not for invalid (non-ACGT) seeds
        if (!seedIn) {
          ctx->validSeeds++;
          // CRC32c hash the seed
#if defined(LOCAL_BUILD) && defined(__x86_64__)
          __asm__ __volatile__(
              "crc32q\t"
              "(%1), %0"
              : "=r"(hash)
              : "r"((uint64_t*)seedPtr), "0"(0));
#elif !defined(LOCAL_BUILD)
          hash = crc32c_hw(0, (const unsigned char*)seedPtr, 8);
#endif
          // Add reference sequence ID in anchored mode, to separate difference sequences
          hash += anchorRefSeq;
          // Obtain a lock corresponding to a 20-bit segment of the hash
          lockPtr  = &seedHashLocks[hash & 0xFFFFF];
          countPtr = &seedHashCounts[hash & 0xFFFFFFFF];
          while (!__sync_bool_compare_and_swap(lockPtr, 0, 1))
            ;
          // Increment count for this hash
          if (*countPtr < minFreqToExtend) {
            // When count reaches the extension threshold, count them all as extended
            if (++*countPtr == minFreqToExtend) {
              ctx->extendedSeeds += minFreqToExtend;
            }
          }
          // And count further seeds with this hash as extended
          else {
            ctx->extendedSeeds++;
          }
          // Release the lock
#if defined(_TARGET_PPC_)
          __asm__ __volatile__("sync" ::: "memory");
#endif
          *lockPtr = 0;
        }
      }
    }
    // skip decoy and population alt contigs
    if (r->isDecoy || r->isPopAlt) continue;

    // Allocate STR records in proportion to sequence length, assuming decimation
    int alloc = (trimLen >> (strDecimationLog2[1][1] - 2)) + 256;
    if (!(r->str = calloc(alloc, sizeof(htStrRec_t)))) continue;
    // Scan the reference sequence
    while (pos < trimLen) {
      bestLen = 0;
      bestEnd = pos + 1;
      // Find the best STR, with greatest repeatLen, breaking ties to lower period
      for (per = 1; per <= STR_MAX_PERIOD; per++) {
        // Skip final partial period
        if ((pos + per - 1) >= trimLen) break;
        // Skip unless [ACGT] through first whole period
        if (!STR_BASE[seq[pos + per - 1]]) break;
        // Leftward extension of repeat with this period
        for (beg = pos - 1; beg > prevBeg; beg--) {
          if (charToBase[seq[beg]] != charToBase[seq[beg + per]]) break;
        }
        beg++;
        // Rightward extension of repeat with this period
        for (end = pos + per; end < trimLen; end++) {
          if (charToBase[seq[end]] != charToBase[seq[end - per]]) break;
        }
        // Track maximum number of whole periods
        len = (end - beg) / per;
        if (len > bestLen) {
          bestLen = len;
          bestPer = per;
          bestBeg = beg;
          bestEnd = end;
        }
      }
      // Process valid STR
      if (bestLen) {
        // Cap repeatLength
        if (bestLen > STR_MAX_REPEATLEN) bestLen = STR_MAX_REPEATLEN;
        bestRun = bestEnd - bestBeg;
        // Test incrementing count per bin against decimation mask
        // Adjust upward by reference ID to randomize downsampling positions, important for many short contigs
        cnt = str.refId + count[bestPer][bestLen]++;
        if (!(cnt & masks[bestPer][bestLen])) {
          // Enlarge allocation 25% if needed
          if (num >= alloc) {
            alloc += alloc >> 2;
            if (!(r->str = realloc(r->str, alloc * sizeof(htStrRec_t)))) break;
          }
          // Save STR record
          str.mask      = cnt;
          str.begPos    = bestBeg + begTrim;
          str.length    = bestRun > 0xFFFF ? 0xFFFF : bestRun;
          str.period    = bestPer;
          str.repeatLen = bestLen;
          r->str[num++] = str;
        }
      }
      pos = bestEnd;
    }
    // Final STR record count
    r->str    = realloc(r->str, (num ? num : 1) * sizeof(htStrRec_t));
    r->strNum = num;
  }
}

int refRecCompareName(const void* a, const void* b)
{
  return strcmp((*(refRec_t**)a)->name, (*(refRec_t**)b)->name);
}

//-------------------------------------------------------------------------------
// Returns sorted array of #numRefSeqs# pointers to ths #refRecs#s. Sorted by name.
// The pointers are allocated with malloc and need to be freed by the caller.
refRec_t** SortRefRecs(refRec_t* refRecs, int32_t numRefSeqs)
{
  refRec_t** sortedRefRecs = malloc(numRefSeqs * sizeof(refRec_t*));
  if (!sortedRefRecs) {
    snprintf(ERR_MSG, sizeof(ERR_MSG), "Cannot allocate sortedSeqNums\n");
    return NULL;
  }
  int i = 0;
  for (i = 0; i < numRefSeqs; i++) {
    sortedRefRecs[i] = &refRecs[i];
  }
  qsort(sortedRefRecs, numRefSeqs, sizeof(refRec_t*), refRecCompareName);
  return sortedRefRecs;
}

void printRefRec(refRec_t* r)
{
  printf("name:      %s\n", r->name);
  printf("seqLen:    %lld\n", r->seqLen);
  printf("begTrim:   %lld\n", r->begTrim);
  printf("endTrim:   %lld\n", r->endTrim);
  printf("trimLen:   %lld\n", r->trimLen);
  printf("endPad:    %lld\n", r->endPad);
  printf("blockLen:  %lld\n", r->blockLen);
  printf("strNum:    %lld\n", r->strNum);
  printf("seqNum:    %lld\n", r->seqNum);
  printf("isPri:     %d\n", r->isPri);
  printf("isAlt:     %d\n", r->isAlt);
  printf("isDecoy:   %d\n", r->isDecoy);
  printf("isPopAlt:  %d\n\n", r->isPopAlt);
  fflush(stdout);
}

//-------------------------------------------------------------------------------
// Compute the number of padding bits for a refRec_t. #trimLen# is the length of the contig after
// trimming undesired bases.
int ComputePadding(uint32_t trimLen, int isPopAlt)
{
  int padLen           = REF_SEQ_MIN_PAD_BASES;
  int refSeqAlignBases = REF_SEQ_ALIGN_BASES;
  if (isPopAlt) {
    // use reduced padding between pop alt contigs
    padLen           = REF_SEQ_MIN_PAD_BASES_POP_ALT;
    refSeqAlignBases = REF_SEQ_ALIGN_BASES_POP_ALT;
  }
  padLen += (refSeqAlignBases - ((trimLen + padLen) % refSeqAlignBases)) % refSeqAlignBases;
  return padLen;
}

//-------------------------------------------------------------------------------
// Compute the number of padding bits for the final refRec_t so that the reference is
// aligned to a 1024 base boundary.
int ComputeEndPadding(uint32_t trimLen, uint64_t totRefLen)
{
  int padLen = REF_SEQ_MIN_PAD_BASES;
  padLen +=
      (REF_SEQ_ALIGN_BASES - ((totRefLen + trimLen + padLen) % REF_SEQ_ALIGN_BASES)) % REF_SEQ_ALIGN_BASES;
  return padLen;
}

//-------------------------------------------------------------------------------
// Mask all intervals in the reference sequences stored in #refRecs# that are in
// the ranged specified in #maskBed# with N's.Updates the trimmed regions and
// padding of each refRec_t. Update #totRefLen# to take into account the new trim
// and padding, and subtracts from #refBasesNotN# the number of non-N bases that
// got masked. Returns the digest of the bed file.
uint32_t MaskBedRegions(
    const char* maskBedPath,
    refRec_t*   refRecs,
    int32_t     numRefSeqs,
    uint64_t*   totRefLen,
    uint64_t*   refBasesNotN)
{
  FILE*    maskBedFile   = NULL;
  uint32_t maskBedDigest = 0;
  if (!(maskBedFile = fopen(maskBedPath, "r"))) {
    snprintf(ERR_MSG, sizeof(ERR_MSG), "Cannot open mask bed file %s\n", maskBedPath);
    return maskBedDigest;
  }

  printf("\nReading BED file %s and masking regions...\n", maskBedPath);

  // Sort the refRecs
  refRec_t** sortedRefRecs = SortRefRecs(refRecs, numRefSeqs);
  if (!sortedRefRecs) {
    fclose(maskBedFile);
    return maskBedDigest;
  }

  const char nCode             = ENCODE_BASE['N'];
  char       bedLine[MAX_LINE] = "";
  int        i;

  // Pointers to key and match for use in the binary search
  refRec_t   r, *searchKey = &r;
  refRec_t** match = NULL;

  while (fgets(bedLine, sizeof(bedLine), maskBedFile)) {
#ifndef LOCAL_BUILD
    maskBedDigest = crc32c_hw(maskBedDigest, bedLine, strlen(bedLine));
#endif
    char* contigName     = strtok(bedLine, " \t");
    char* startPosString = strtok(NULL, " \t");
    char* endPosString   = strtok(NULL, " \t\n");

    if (!(contigName && startPosString && endPosString)) {
      snprintf(ERR_MSG, sizeof(ERR_MSG), "Invalid line in mask bed file.");
      break;
    }

    // Find the refRec_t with this contig name
    searchKey->name = contigName;
    const refRec_t* curRefRec =
        (match =
             (refRec_t**)bsearch(&searchKey, sortedRefRecs, numRefSeqs, sizeof(refRec_t*), refRecCompareName))
            ? *match
            : NULL;
    if (!match) {
      snprintf(ERR_MSG, sizeof(ERR_MSG), "Invalid contig name in mask bed file: %s\n", contigName);
      break;
    }

    uint32_t startPos = strtoul(startPosString, NULL, 10);  // 0-based
    uint32_t endPos   = strtoul(endPosString, NULL, 10);    // 1-based

    if (endPos > curRefRec->seqLen) {
      snprintf(
          ERR_MSG,
          sizeof(ERR_MSG),
          "Position in mask bed file (%u) out of range (%u) in contig %s",
          endPos,
          curRefRec->seqLen,
          curRefRec->name);
      break;
    }

    if (endPos < startPos) continue;
    size_t len = endPos - startPos;

    // Adjust the position for the trimmed N's from the front
    if (startPos < curRefRec->begTrim) {
      uint32_t lengthAdjustment = curRefRec->begTrim - startPos;
      APPLY_MAX(lengthAdjustment, len);
      startPos = 0;
      len -= lengthAdjustment;
    } else {
      startPos -= curRefRec->begTrim;
    }

    // And for the trimmed N's from the back
    APPLY_MAX(startPos, curRefRec->trimLen);
    APPLY_MAX(len, curRefRec->trimLen - startPos);

    for (i = 0; i < len; i++) {
      // Only consider it a masked base if it wasnt N before the masking
      if (ENCODE_BASE[curRefRec->seq[startPos + i]] < nCode) {
        --(*refBasesNotN);
      }
      curRefRec->seq[startPos + i] = 'N';
    }

#ifndef LOCAL_BUILD
    WatchDogCheckin(watchdogID);
#endif
  }
  fclose(maskBedFile);
  free(sortedRefRecs);

  // If the starts or ends of a contig were masked, we need to re-trim it.
  // Reset the total reference length after trimming. As we iterate through the sequences and re-trim
  // them, we tally up the new total reference length.
  *totRefLen = REF_SEQ_END_PAD_BASES;
  for (i = 0; i < numRefSeqs; i++) {
    const uint32_t oldLen = refRecs[i].trimLen;
    char*          oldSeq = refRecs[i].seq;

    uint32_t numLeadingNs = 0;
    while ((numLeadingNs < oldLen) && (ENCODE_BASE[oldSeq[numLeadingNs]] == nCode)) ++numLeadingNs;

    uint32_t numTrailingNs = 0;
    while ((numTrailingNs + numLeadingNs < oldLen) &&
           (ENCODE_BASE[oldSeq[oldLen - numTrailingNs - 1]] == nCode))
      ++numTrailingNs;

    // Updated the lengths of the trimming and padding
    if (numLeadingNs == oldLen) {
      refRecs[i].begTrim = 0;
      refRecs[i].endTrim = refRecs[i].seqLen;
    } else {
      refRecs[i].begTrim += (numLeadingNs / REF_SEQ_TRIM_GRAN) * REF_SEQ_TRIM_GRAN;
      refRecs[i].endTrim += numTrailingNs;
    }

    refRecs[i].trimLen = refRecs[i].seqLen - refRecs[i].endTrim - refRecs[i].begTrim;
    refRecs[i].endPad  = (i + 1 == numRefSeqs && refRecs[i].isPopAlt)
                            ? ComputeEndPadding(refRecs[i].trimLen, *totRefLen)
                            : ComputePadding(refRecs[i].trimLen, refRecs[i].isPopAlt);
    refRecs[i].blockLen = refRecs[i].trimLen + refRecs[i].endPad;

    (*totRefLen) += refRecs[i].blockLen;

    // Avoid copying the sequence if it has not changed.
    if (refRecs[i].trimLen != oldLen) {
      // Trim the actual sequence
      refRecs[i].seq = malloc(refRecs[i].trimLen);
      memcpy(refRecs[i].seq, oldSeq + refRecs[i].begTrim, refRecs[i].trimLen);
      free(oldSeq);
    }
  }
  return maskBedDigest;
}

#ifndef LOCAL_BUILD
//-------------------------------------------------------------------------------
// Validate that there is some input for dealing with hg19/38 alt contigs if they
// are present. If neither a mask bed nor liftover sam are provided, will attempt
// to autodetect bed file.
// On error, an exeption is thrown from the c++ code.
void AltContigValidate(hashTableConfig_t* config, const refRec_t* refRecs, int32_t numRefSeqs)
{
  int haveAltAwareInput = (config->maskBed != NULL);



  // If no alt aware input was given and autodetect is enabled, try to find a mask bed
  // file for the alt contigs
  if (!haveAltAwareInput && config->autoDetectValidate) {
    const char* autoDetectedAltMask = AutoDetectAltMask(config->autoDetectDir);
    if (autoDetectedAltMask) {
      config->maskBed   = strdup(autoDetectedAltMask);
      haveAltAwareInput = 1;
    }
  }

  if (config->altContigValidate && !haveAltAwareInput) {
    // Not building hash table with liftover - setup to detect and report error for
    // hg19/hg38 references
    AltContigTrackerInit(config->refInput);

    // Loop over the reference sequences and feed them to the alt contig tracker
    // Throws an exception if hg19/hg38 (without liftover) is being built
    int i = 0;
    for (i = 0; i < numRefSeqs; i++) AltContigTrackerTally(refRecs[i].name);
  }
}
#endif

//-------------------------------------------------------------------------swhitmore/miker
// generateHashTable - Returns NULL on success, or error message
//
char* generateHashTable(hashTableConfig_t* config, int argc, char* argv[])
{
  uint64_t refSeqLen = 0, refIdxLen = 0, refSeqAlloc, refMaskAlloc, refCodeAlloc, nextReport, popSnpsLen = 0;
  uint64_t altMatchesAlloc = 0;
  uint8_t *refSeq = NULL, *refMask = NULL, *refCode = NULL, *altMatches = NULL;
  char     c, line[MAX_LINE] = "", seqTag[MAX_LINE] = "";
  FILE *   outRefFile = NULL, *inpRefFile = NULL, *refIdxFile = NULL, *repMaskFile = NULL;
  int i, n, base, code, seqByte = 0, maskByte = 0, refByte = 0, repMaskByte = 0, repMaskBits = 0, done = 0,
                        addrBits;
  uint64_t refCodeHist[16] = {0};
  uint32_t refDigest = 0, refIndexDigest = 0, hashDigest = 0, extTabDigest = 0, popSnpsDigest = 0;
  char     nullPaddingBlock[REF_SEQ_END_PAD_BASES] = {0};

  // Get the maximum number of threads in system if it is greater than 48 use 48 if it is less than 32 use 32
  uint32_t maxThreadsInSystem = ((sysconf(_SC_NPROCESSORS_ONLN) > 48) ? 48 : sysconf(_SC_NPROCESSORS_ONLN));
  maxThreadsInSystem          = ((maxThreadsInSystem < 32) ? 32 : maxThreadsInSystem);
  uint32_t minThreadsInSystem = sysconf(_SC_NPROCESSORS_ONLN);

  hashTableHeader_t* cfghdr = config->hdr;

  APPLY_MIN(cfghdr->refSeedInterval, 1.0);
  APPLY_MAX(cfghdr->refSeedInterval, 255.9375);
  cfghdr->refSeedInterval = round(cfghdr->refSeedInterval * 16) / 16;
  APPLY_MIN(cfghdr->maxSeedFreq, 1);
  APPLY_MAX(cfghdr->maxSeedFreq, MAX_SEED_HIT_FREQ);
  if (cfghdr->priMaxSeedFreq < 1 || cfghdr->priMaxSeedFreq > cfghdr->maxSeedFreq)
    cfghdr->priMaxSeedFreq = cfghdr->maxSeedFreq;
  APPLY_MINMAX(cfghdr->targetSeedFreq, 1, cfghdr->maxSeedFreq);
  APPLY_MIN(cfghdr->priSeedBases, 1);
  APPLY_MINMAX(cfghdr->maxExtIncrement, 2, MAX_SEED_EXTENSION_INCR);
  APPLY_MIN(cfghdr->thinningFreqCap, 1);
  APPLY_MINMAX(cfghdr->thinningPeriod, 1, THINNING_MAX_PERIOD);
  APPLY_MIN(cfghdr->seedLenCost, 0.0);
  APPLY_MIN(cfghdr->seedFreqCost, 0.0);
  APPLY_MIN(cfghdr->extensionCost, 0.0);
  APPLY_MIN(cfghdr->extStepCost, 0.0);
  APPLY_MIN(cfghdr->extRecCost, 0.0);
  APPLY_MIN(cfghdr->extRandHitFreq, 0);
  if (config->maxThreads)
    APPLY_MINMAX(config->maxThreads, 1, maxThreadsInSystem);
  else
    APPLY_MINMAX(config->maxThreads, minThreadsInSystem, maxThreadsInSystem);
  APPLY_MINMAX(cfghdr->minRepairProb, 0.0, 1.0);
  APPLY_MIN(cfghdr->anchorBinBits, 0);
  APPLY_MINMAX(config->writeHashFile, 0, 1);
  APPLY_MINMAX(config->writeCompFile, 0, 1);

  // By default, if the user increases maxThreads, increase maxGB to match
  if (config->maxGB <= 0) {
    config->maxGB = 32;
    if (config->maxThreads > config->maxGB) {
      config->maxGB = config->maxThreads;
    }
  }

  struct stat s;
  int         ret      = stat(config->refInput, &s);
  int64_t     inputLen = ret < 0 ? -1 : s.st_size;
  if (!(inpRefFile = fopen(config->refInput, "r"))) {
    snprintf(ERR_MSG, sizeof(ERR_MSG), "Cannot open reference input file %s\n", config->refInput);
    return ERR_MSG;
  }

#ifndef LOCAL_BUILD
  watchdogID = WatchDogRegister("_gen_hash_table");

  // Create an object for quick lookup by reference sequence name. This will make it easier for us to check
  // for duplicate names.
  ReferenceNamesInit();
#endif

  gzFile inpRefGzFile = NULL;

  // Read input FASTA (reference, decoy, pop alt contig) into sequence records
  printf("\nReading reference input file %s", config->refInput);
  if (inputLen >= 0) printf(": %lld bytes", inputLen);
  printf("...\n");
  fflush(stdout);
  int32_t   numRefSeqs = -1, origRefSeqs = -1;
  uint32_t  decoySeqs = 0, decoySeqLen = 0;
  uint32_t  popAltSeqs = 0, popAltSeqLen = 0;
  uint64_t  curSeqLen = 0, totSeqLen = 0, totRefLen = REF_SEQ_END_PAD_BASES, curBasesNotN = 0;
  uint64_t  refBasesNotN = 0;
  int64_t   adjBasesNotN = 0;
  uint32_t  inpSeqAlloc = 0, refRecsAlloc = 0;
  int64_t   firstNotN = -1, lastNotN = -1;
  char*     inpSeq                    = NULL;
  int       grchCounts[GRCH_VERSIONS] = {0};
  int       grchVersion               = 0;
  int       startDecoys               = 0;
  int       endDecoys                 = 0;
  int       inDecoys                  = 0;
  int       startPopAlts              = 0;
  int       inPopAlts                 = 0;
  refRec_t* refRecs                   = NULL;
  done                                = 0;
  nextReport                          = SEQ_REPORT_INTERVAL;
  while (!done) {
    char* ret = NULL;
    if (inpRefFile) {
      ret = fgets(line, sizeof(line), inpRefFile);
    } else if (inpRefGzFile) {
      ret = gzgets(inpRefGzFile, line, sizeof(line));
    }
    if (!ret) {
      if (inpRefFile) {
        fclose(inpRefFile);
      }
      if (inpRefGzFile) {
        gzclose(inpRefGzFile);
      }
#if 0
      printf(
          "  %d sequence%s, %lld bases (%lld after trimming/padding)\n",
          numRefSeqs,
          numRefSeqs == 1 ? "" : "s",
          totSeqLen,
          totRefLen);
      fflush(stdout);
#endif

      if (!startDecoys && !endDecoys && config->decoyFname) {
        if (!(inpRefFile = fopen(config->decoyFname, "r"))) {
          snprintf(ERR_MSG, sizeof(ERR_MSG), "Cannot open automatic decoys file %s\n", config->decoyFname);
          return ERR_MSG;
        }
        // printf("\nReading decoys input file %s...\n", config->decoyFname);
        // Decide GRCh version
        for (n = 0; n < GRCH_VERSIONS; n++)
          if (grchCounts[n] >= 22) grchVersion = GRCH_VERSION_LO + n;
        startDecoys = 1;
        continue;
      } else if (!startPopAlts && config->popAltContigsFname) {
        // Check if file ends in .gz
        if (strcmp(config->popAltContigsFname + strlen(config->popAltContigsFname) - 3, ".gz") == 0) {
          if (!(inpRefGzFile = gzopen(config->popAltContigsFname, "rb"))) {
            snprintf(
                ERR_MSG,
                sizeof(ERR_MSG),
                "Cannot open population alternate contigs file %s\n",
                config->popAltContigsFname);
            return ERR_MSG;
          }
          inpRefFile = NULL;
        } else {
          if (!(inpRefFile = fopen(config->popAltContigsFname, "r"))) {
            snprintf(
                ERR_MSG,
                sizeof(ERR_MSG),
                "Cannot open population alternate contigs file %s\n",
                config->popAltContigsFname);
            return ERR_MSG;
          }
          inpRefGzFile = NULL;
        }
        printf("\nReading population alternate contigs file %s...\n", config->popAltContigsFname);
        startDecoys  = 0;
        endDecoys    = 1;
        startPopAlts = 1;
        continue;
      } else {
        done = 1;
      }
    }
    if (done || line[0] == '>') {
      int suppress = 0;
      if (numRefSeqs >= 0) {
        if (numRefSeqs >= refRecsAlloc) {
          refRecsAlloc = (refRecsAlloc ? refRecsAlloc * 2 : 1024);
          refRecs      = realloc(refRecs, refRecsAlloc * sizeof(*refRecs));
          if (!refRecs) {
            snprintf(ERR_MSG, sizeof(ERR_MSG), "Cannot allocate %lu reference records\n", refRecsAlloc);
            return ERR_MSG;
          }
        }
        if (curSeqLen > 0xFFFFFFFF) {
          snprintf(
              ERR_MSG,
              sizeof(ERR_MSG),
              "Reference sequence #%d exceeds 2^32-1 bases (4.29Gbp): %s\n",
              numRefSeqs,
              seqTag);
          return ERR_MSG;
        }
        refRec_t* r = &refRecs[numRefSeqs];
        r->seqLen   = curSeqLen;
        r->begTrim  = firstNotN < 0 ? 0 : (firstNotN / REF_SEQ_TRIM_GRAN) * REF_SEQ_TRIM_GRAN;
        r->endTrim  = curSeqLen - lastNotN - 1;
        r->trimLen  = r->seqLen - r->begTrim - r->endTrim;
        r->endPad   = (done && inPopAlts) ? ComputeEndPadding(r->trimLen, totRefLen)
                                        : ComputePadding(r->trimLen, inPopAlts);
        r->blockLen = r->trimLen + r->endPad;
        r->isPri    = 0;
        r->isAlt    = 0;
        r->isDecoy  = inDecoys;
        r->isPopAlt = inPopAlts;
        // Detect GRCh37,.. by testing sequence lengths
        if (!inDecoys)
          for (n = 0; n < GRCH_VERSIONS; n++)
            for (i = 0; i < 24; i++) grchCounts[n] += (curSeqLen == grchChromLens[n][i]);
        // Detect standard Edico decoy sequences from GRCh*
        if (inDecoys && !strncmp(seqTag, "Edico_decoy_GRCh", 16) && isdigit(seqTag[16]) &&
            isdigit(seqTag[17]) && seqTag[18] == '_') {
          n = (seqTag[16] - '0') * 10 + (seqTag[17] - '0');
          if (n >= GRCH_VERSION_LO && n <= GRCH_VERSION_HI) {
            char* decoyName = seqTag + 19;
            // Suppress this decoy sequence if not using its source GRCh version, or its name is already in
            // the input FASTA
            suppress = (n != grchVersion);
#ifdef LOCAL_BUILD
            for (i = 0; i < origRefSeqs; i++)
              if (!strcmp(decoyName, refRecs[i].name)) suppress = 1;
#else
            if (ReferenceNamesContains(decoyName)) {
              suppress = 1;
            }
#endif
          }
        }
        // Silently supress pop alt duplicates so that we dont need to add them to the autodetector.
#ifndef LOCAL_BUILD
        if (inPopAlts && ReferenceNamesContains(seqTag)) suppress = 1;
#endif
        if (!suppress) {
          r->name      = malloc(strlen(seqTag) + 1);
          r->seq       = malloc(r->trimLen);
          r->seqNum    = numRefSeqs;
          r->str       = NULL;
          r->lifted    = NULL;
          r->liftMatch = NULL;
          r->strNum    = 0;
          if (!r->name || !r->seq) {
            snprintf(ERR_MSG, sizeof(ERR_MSG), "Cannot allocate %lu byte reference sequence\n", r->trimLen);
            return ERR_MSG;
          }
          strcpy(r->name, seqTag);
          memcpy(r->seq, &inpSeq[r->begTrim], r->trimLen);
          totSeqLen += curSeqLen * !inDecoys;
          totRefLen += r->blockLen;
          decoySeqs += inDecoys;
          decoySeqLen += inDecoys * curSeqLen;
          popAltSeqs += inPopAlts;
          popAltSeqLen += inPopAlts * curSeqLen;
          refBasesNotN += curBasesNotN;

#ifndef LOCAL_BUILD
          // This will throw an exception if the reference sequence name has already been inserted (if the
          // same name is used multiple times in the input fasta).
          if (!inDecoys && !inPopAlts) {
            ReferenceNamesInsert(r->name, r->seqNum, r->seqLen);
          }
#endif
        }
        curBasesNotN = 0;
        curSeqLen    = 0;
        firstNotN    = -1;
        lastNotN     = -1;
      }
      numRefSeqs += !suppress;
      origRefSeqs += !inDecoys & !inPopAlts;
      inDecoys  = startDecoys;
      inPopAlts = startPopAlts;
      if (!done) {
        strcpy(seqTag, line + 1);
        for (n = strlen(seqTag) - 1; n >= 0; n--)
          if (isspace(seqTag[n])) seqTag[n] = 0;
      }
      if (totSeqLen >= nextReport) {
        printf(
            "  %d sequence%s, %lld bases (%lld after trimming/padding)\n",
            origRefSeqs,
            origRefSeqs == 1 ? "" : "s",
            totSeqLen,
            totRefLen);
        fflush(stdout);
        nextReport += SEQ_REPORT_INTERVAL;
      }
      continue;
    } else if (line[0] == ';' || numRefSeqs < 0)
      continue;
    if (curSeqLen + sizeof(line) > inpSeqAlloc) {
      inpSeqAlloc = (inpSeqAlloc ? inpSeqAlloc * 2 : 64 * 1024 * 1024);
      inpSeq      = realloc(inpSeq, inpSeqAlloc);
      if (!inpSeq) {
        snprintf(ERR_MSG, sizeof(ERR_MSG), "Cannot allocate %llu bytes for liftover text\n", inpSeqAlloc);
        return ERR_MSG;
      }
    }
    char* p = line;
    int   notN;
    while ((c = *p++)) {
      code = ENCODE_BASE[c];
      notN = (code > 0 && code < 0xF);
      curBasesNotN += notN;
      if (notN) {
        lastNotN = curSeqLen;
        if (firstNotN < 0) firstNotN = curSeqLen;
      }
      inpSeq[curSeqLen] = c;
      curSeqLen += (code >= 0);
    }
#ifndef LOCAL_BUILD
    WatchDogCheckin(watchdogID);
#endif
  }

  // Done extracting reference sequences from input fastas.
  if (origRefSeqs + decoySeqs + popAltSeqs != numRefSeqs) {
    snprintf(ERR_MSG, sizeof(ERR_MSG), "Invalid number of sequences\n");
    return ERR_MSG;
  }

  free(inpSeq);
  inpSeq = NULL;

  // If there are known alt contigs present and the user wants this validation,
  // make sure there is some alt-aware method available (masking or liftover)
  // And if autodetect is enabled, autodetect a mask bed if it is needed.
#ifndef LOCAL_BUILD
  AltContigValidate(config, refRecs, numRefSeqs);
#endif

  // Sort pointers to reference sequences by name for faster lookup
  refRec_t** sortedRefRecs = SortRefRecs(refRecs, numRefSeqs);
  if (!sortedRefRecs) {
    return ERR_MSG;
  }

  // TODO: make this function take the sorted refRecs instead of re-sorting.
  // If a BED file for masking is provided, change bases in the intervals to N's
  config->maskBedDigest = 0;
  if (config->maskBed) {
    ERR_MSG[0]            = '\0';
    config->maskBedDigest = MaskBedRegions(config->maskBed, refRecs, numRefSeqs, &totRefLen, &refBasesNotN);
    if (strlen(ERR_MSG)) return ERR_MSG;
  }

  // Pointer for getting error messages from helper functions.
  char* returnedErrMsg = NULL;

  uint32_t liftoverDigest = 0;


  // After contigs and liftover files have been read and liftover relations established, if we are building a
  // double methyl converted hash table, add and modify the different conversion types to refRecs, and modify
  // global size parameters
  if (config->methylatedConv == HT_TYPE_METHYL_COMBINED) {
    printf(
        "\nBefore methylation conversions: %d sequence%s, %lld bases (%lld after trimming/padding)\n",
        numRefSeqs,
        numRefSeqs == 1 ? "" : "s",
        totSeqLen,
        totRefLen);

    refRecsAlloc *= 2;
    refRecs = realloc(refRecs, refRecsAlloc * sizeof(*refRecs));
    if (!refRecs) {
      snprintf(ERR_MSG, sizeof(ERR_MSG), "Cannot allocate %lu reference records\n", refRecsAlloc);
      return ERR_MSG;
    }

    printf("Applying methylation conversions to contigs...\n");
    fflush(stdout);
    UpdateContigsForMethylation(refRecs, numRefSeqs, &totRefLen);



    totSeqLen *= 2;
    decoySeqs *= 2;
    decoySeqLen *= 2;
    popAltSeqs *= 2;
    popAltSeqLen *= 2;
    refBasesNotN *= 2;
    numRefSeqs *= 2;
    origRefSeqs *= 2;
  }

  totRefLen += REF_SEQ_END_PAD_BASES;
  printf(
      "\nTotal: %d sequence%s, %lld bases (%lld after trimming/padding)\n",
      numRefSeqs,
      numRefSeqs == 1 ? "" : "s",
      totSeqLen,
      totRefLen);

#ifndef LOCAL_BUILD
  char buf[1024];
  snprintf(
      buf,
      sizeof(buf),
      "Hash generation: included %u extra decoy sequences, %u bases\n",
      decoySeqs,
      decoySeqLen);
  SysLoggerLogInfo(buf);
  snprintf(
      buf,
      sizeof(buf),
      "Hash generation: included %u extra population alt sequences, %u bases\n",
      popAltSeqs,
      popAltSeqLen);
  SysLoggerLogInfo(buf);
#endif
  if (totRefLen >> MAPPER_REF_POS_BITS != 0) {
    snprintf(
        ERR_MSG,
        sizeof(ERR_MSG),
        "Padded reference sequence too long (%lld); maximum %lld supported\n",
        totRefLen,
        1ull << MAPPER_REF_POS_BITS);
    return ERR_MSG;
  }
  if (origRefSeqs > MAX_REF_SEQS) {
    snprintf(
        ERR_MSG,
        sizeof(ERR_MSG),
        "Too many reference sequences (%d); only %d supported\n",
        origRefSeqs,
        MAX_REF_SEQS);
    return ERR_MSG;
  }

  // Each refRec gets a lookup table to translate from the characters read from the FASTQ file,
  // into the accepted character values.  This is a quick way to interpret upper- or
  // lower-case letters the same (e.g. 'a' is equivalent to 'A').
  //
  // For the case where the user is doing processing of methylated DNA, they may want
  // to force all G's to be interpreted as A's, or all C's as T's.  This lookup table
  // is a great place to force this conversion.
  // Convert multi-base codes too, replacing the source base with the destination
  // base in each matching subset.  But leave 'N' alone because it has special
  // status in alignment scoring.

  // Default lookup table, no base conversions.
  char BASE_ENCODING[256];
  if (cfghdr->maxMultBaseSeeds > 0) {
    // Negative values not wanted here
    for (i = 0; i < 256; i++) BASE_ENCODING[i] = ENCODE_BASE[i] < 0 ? 0 : ENCODE_BASE[i];
  } else {
    memcpy(BASE_ENCODING, ENCODE_BASE, 256);
  }

  // Conversion tables used for methylation.
  char BASE_ENCODING_G_TO_A[256];
  char BASE_ENCODING_C_TO_T[256];
  memcpy(BASE_ENCODING_G_TO_A, BASE_ENCODING, 256);
  memcpy(BASE_ENCODING_C_TO_T, BASE_ENCODING, 256);

  if (cfghdr->maxMultBaseSeeds > 0) {
    BASE_ENCODING_G_TO_A['g'] = BASE_ENCODING_G_TO_A['G'] = BASE_ENCODING_G_TO_A['A'];
    BASE_ENCODING_G_TO_A['r'] = BASE_ENCODING_G_TO_A['R'] = BASE_ENCODING_G_TO_A['A'];
    BASE_ENCODING_G_TO_A['s'] = BASE_ENCODING_G_TO_A['S'] = BASE_ENCODING_G_TO_A['M'];
    BASE_ENCODING_G_TO_A['v'] = BASE_ENCODING_G_TO_A['V'] = BASE_ENCODING_G_TO_A['M'];
    BASE_ENCODING_G_TO_A['k'] = BASE_ENCODING_G_TO_A['K'] = BASE_ENCODING_G_TO_A['W'];
    BASE_ENCODING_G_TO_A['d'] = BASE_ENCODING_G_TO_A['D'] = BASE_ENCODING_G_TO_A['W'];
    BASE_ENCODING_G_TO_A['b'] = BASE_ENCODING_G_TO_A['B'] = BASE_ENCODING_G_TO_A['H'];

    BASE_ENCODING_C_TO_T['c'] = BASE_ENCODING_C_TO_T['C'] = BASE_ENCODING_C_TO_T['T'];
    BASE_ENCODING_C_TO_T['m'] = BASE_ENCODING_C_TO_T['M'] = BASE_ENCODING_C_TO_T['W'];
    BASE_ENCODING_C_TO_T['s'] = BASE_ENCODING_C_TO_T['S'] = BASE_ENCODING_C_TO_T['K'];
    BASE_ENCODING_C_TO_T['v'] = BASE_ENCODING_C_TO_T['V'] = BASE_ENCODING_C_TO_T['D'];
    BASE_ENCODING_C_TO_T['y'] = BASE_ENCODING_C_TO_T['Y'] = BASE_ENCODING_C_TO_T['T'];
    BASE_ENCODING_C_TO_T['h'] = BASE_ENCODING_C_TO_T['H'] = BASE_ENCODING_C_TO_T['W'];
    BASE_ENCODING_C_TO_T['b'] = BASE_ENCODING_C_TO_T['B'] = BASE_ENCODING_C_TO_T['K'];
  } else {
    BASE_ENCODING_G_TO_A['g'] = BASE_ENCODING_G_TO_A['a'];
    BASE_ENCODING_G_TO_A['G'] = BASE_ENCODING_G_TO_A['A'];

    BASE_ENCODING_C_TO_T['c'] = BASE_ENCODING_C_TO_T['t'];
    BASE_ENCODING_C_TO_T['C'] = BASE_ENCODING_C_TO_T['T'];
  }

  // Construct a code-to-base lookup table
  char codeToBase[16];
  memset(codeToBase, 4, 16);
  codeToBase[BASE_A] = 0;
  codeToBase[BASE_C] = 1;
  codeToBase[BASE_G] = 2;
  codeToBase[BASE_T] = 3;

  // This is only used when (maxMultBaseSeeds == 0)
  // And a direct char-to-base table
  int charToBase[256];
  int charToBaseGToA[256];
  int charToBaseCToT[256];
  for (i = 0; i < 256; i++) {
    // For those cases where BASE_ENCODING returns -1, index codeToBase by BASE_PAD
    // (i.e. 0), which returns the default value of 4.
    charToBase[i]     = codeToBase[(BASE_ENCODING[i] == -1 ? BASE_PAD : BASE_ENCODING[i])];
    charToBaseGToA[i] = codeToBase[(BASE_ENCODING_G_TO_A[i] == -1 ? BASE_PAD : BASE_ENCODING_G_TO_A[i])];
    charToBaseCToT[i] = codeToBase[(BASE_ENCODING_C_TO_T[i] == -1 ? BASE_PAD : BASE_ENCODING_C_TO_T[i])];
  }

  // Set each refRec with its lookup table and set a sequence name suffix
  // in the case of a single methylation table. For the double methylation table,
  // the contig names are already updated in the refRecs.
  const char* seqnameSuffix = seqnameSuffixDefault;
  if (config->methylatedConv == HT_TYPE_METHYL_COMBINED) {
    for (i = 0; i < numRefSeqs; i++) {
      refRecs[i].charToBase = !(i % 2) ? charToBaseCToT : charToBaseGToA;
      refRecs[i].charToCode = !(i % 2) ? BASE_ENCODING_C_TO_T : BASE_ENCODING_G_TO_A;
    }
  } else if (config->methylatedConv == HT_TYPE_METHYL_G_TO_A) {
    for (i = 0; i < numRefSeqs; i++) {
      refRecs[i].charToBase = charToBaseGToA;
      refRecs[i].charToCode = BASE_ENCODING_G_TO_A;
    }
    seqnameSuffix = seqnameSuffixGtoA;
  } else if (config->methylatedConv == HT_TYPE_METHYL_C_TO_T) {
    for (i = 0; i < numRefSeqs; i++) {
      refRecs[i].charToBase = charToBaseCToT;
      refRecs[i].charToCode = BASE_ENCODING_C_TO_T;
    }
    seqnameSuffix = seqnameSuffixCtoT;
  } else {
    for (i = 0; i < numRefSeqs; i++) {
      refRecs[i].charToBase = charToBase;
      refRecs[i].charToCode = BASE_ENCODING;
    }
  }

  uint64_t altMatchSeeds = 0;
  // Next hash table version, always do this
  if (cfghdr->maxMultBaseSeeds > 0) {

    // Deduct 4/5 of the alt contig positions matching liftover positions from the non-N base count,
    // because we will populate those seeds with only 20% normal density
    cfghdr->liftMatchSeedInt = 5;
    adjBasesNotN = refBasesNotN - (altMatchSeeds * (cfghdr->liftMatchSeedInt - 1) / cfghdr->liftMatchSeedInt);
    APPLY_MIN(adjBasesNotN, 0);
  }
  // Next hash table version, eliminate this
  else {
    cfghdr->liftMatchSeedInt = 0;
    adjBasesNotN             = refBasesNotN;
  }

  // Begin size calculations
  uint64_t bytes;
  if (!decodeSizeArg(config->sizeStr, &bytes)) {
    snprintf(ERR_MSG, sizeof(ERR_MSG), "Could not decode %s\n", config->sizeStr);
    return ERR_MSG;
  }

  // Default space to reserve for annotated SJ table
  uint64_t sjReserve = totRefLen / REF_BASES_PER_SJ_RESERVE_BYTE;
  // User override
  if (config->sjSizeStr && (isdigit(config->sjSizeStr[0]) || config->sjSizeStr[0] == '.')) {
    if (!decodeSizeArg(config->sjSizeStr, &sjReserve)) {
      snprintf(ERR_MSG, sizeof(ERR_MSG), "Could not decode SJ reserve string %s\n", config->sjSizeStr);
      return ERR_MSG;
    }
  }

  // Convert hash table size string to bytes, address bits, and 64ths
  if (bytes) {
    addrBits = 0;
    int sixtyfourths;
    while (bytes && !(bytes & 1)) {
      bytes >>= 1;
      addrBits++;
    }
    if (bytes < 1 || bytes > 63) {
      snprintf(
          ERR_MSG,
          sizeof(ERR_MSG),
          "Invalid table size '%s'.  Must be a power of two times an integer under 64.\n",
          config->sizeStr);
      return ERR_MSG;
    }
    while (bytes <= 32) {
      bytes <<= 1;
      addrBits--;
    }
    addrBits += 6;
    sixtyfourths           = bytes;
    cfghdr->tableAddrBits  = addrBits;
    cfghdr->tableSize64ths = sixtyfourths;
  }

  printf("\n");

  // Things to fit in memory limit besides hash table and extension table
  uint64_t refSize      = totRefLen / 2;
  uint64_t idxSize      = totRefLen / REF_BASES_PER_IDX_RESERVE_BYTE;
  uint64_t totOtherSize = sjReserve + refSize + idxSize;

  // Decode memory size string
  uint64_t cfgMemSize;
  if (!decodeSizeArg(config->memSizeStr, &cfgMemSize)) {
    snprintf(ERR_MSG, sizeof(ERR_MSG), "Could not decode memory size string %s\n", config->memSizeStr);
    return ERR_MSG;
  }
  uint64_t memSize = cfgMemSize;
  uint64_t defaultSize;
  if (!decodeSizeArg(DEFAULT_MEM_SIZE_STR, &defaultSize)) {
    snprintf(ERR_MSG, sizeof(ERR_MSG), "Could not decode size string %s\n", DEFAULT_MEM_SIZE_STR);
    return ERR_MSG;
  }
  // If not specified, come up with a generous estimate with a max of 32GB (DEFAULT)
  if (!memSize) {
    int64_t generousSize = adjBasesNotN * 3 * HASH_RECORD_BYTES / cfghdr->refSeedInterval;
    if (generousSize < 4 * MAX_WRAP_BYTES) generousSize = 4 * MAX_WRAP_BYTES;
    memSize = generousSize + totOtherSize;
    if (memSize > defaultSize) memSize = defaultSize;
  }
  uint64_t maxSize;
  if (!decodeSizeArg(MAX_MEM_SIZE_STR, &maxSize)) {
    snprintf(ERR_MSG, sizeof(ERR_MSG), "Could not decode size string %s\n", MAX_MEM_SIZE_STR);
    return ERR_MSG;
  }
  if (memSize > maxSize) memSize = maxSize;

  // Initial assumed hash table capacity
  int64_t tableBytes = bytes ? bytes : memSize - totOtherSize;

  // Adjust seed interval for large references, if not explicitly specified
  double rawOcc   = (double)adjBasesNotN / (tableBytes / HASH_RECORD_BYTES);
  double reduceBy = rawOcc <= THRESH_OCCUPANCY ? 1.0 : rawOcc / TARGET_OCCUPANCY;
  double overlong = (double)totRefLen / MAX_SEED_INDEXES;
  if (overlong > reduceBy) reduceBy = overlong;
  if (reduceBy > 1.0 && cfghdr->refSeedInterval == 1.0) {
    reduceBy = ceil(reduceBy * 16) / 16;
    APPLY_MAX(reduceBy, 255.9375);
    cfghdr->refSeedInterval = reduceBy;
    printf("Increasing seed interval to %g for long reference\n\n", cfghdr->refSeedInterval);
  }
  if (totRefLen / cfghdr->refSeedInterval > MAX_SEED_INDEXES) {
    overlong = ceil(overlong * 16) / 16;
    snprintf(
        ERR_MSG,
        sizeof(ERR_MSG),
        "These reference sequences require minimum ht-ref-seed-interval %g\n",
        overlong);
    return ERR_MSG;
  }
  // When the selected reference seed interval is exactly 1+3/16=1.1875, populating every 5th seed index
  // where alt contigs match their liftover positions is problematic.  5*1.1875 = 5.9375 is close
  // to an even integer, and since the mapper looks up every 2nd seed position by default, it could
  // miss populated seeds for a long stretch.  Populate every 6th seed index instead in that case.
  // Other refSeedInterval values up to 1.5 are okay, and above that assumedly we'd have to query
  // seeds more densely at run time.
  if (cfghdr->liftMatchSeedInt == 5 && cfghdr->refSeedInterval == 1.1875) {
    cfghdr->liftMatchSeedInt = 6;
  }

  // Now that we have refSeedInterval finalized, interrupt size calculations to load reference sequences
  // and find out how large the seed extension table needs to be...

  // Estimate what primary seed frequency the cost function will typically want to extend,
  // and save that value as a firm policy never to extend any lower-frequency seeds,
  // to help ensure the seed extension table will not overflow what we allocate
  cfghdr->minFreqToExtend =
      cfghdr->seedFreqCost > 0
          ? cfghdr->targetSeedFreq +
                round((cfghdr->extStepCost + 6 * cfghdr->seedLenCost) / cfghdr->seedFreqCost)
          : 999;
  if (cfghdr->minFreqToExtend > cfghdr->priMaxSeedFreq) cfghdr->minFreqToExtend = cfghdr->priMaxSeedFreq + 1;

  // Spawn threads to scan reference sequences for short tandem repeats (STRs)
  int             allowSeedExt  = cfghdr->maxSeedBases > cfghdr->priSeedBases;
  int             useStrThreads = config->strFname || allowSeedExt;
  int             numStrThreads = config->maxThreads <= 1 ? 1 : config->maxThreads - 1;
  int32_t         nextRefSeq    = 0;
  pthread_t       strThreads[numStrThreads];
  strScanCtx_t    strCtx[numStrThreads];
  pthread_mutex_t strMut;
  uint8_t*        seedHashCounts = NULL;
  uint8_t*        seedHashLocks  = NULL;
  if (useStrThreads) {
    // Allocate seed hash counters
    seedHashCounts = calloc(1ULL << 32, 1);
    seedHashLocks  = calloc(1ULL << 20, 1);
    if (!seedHashCounts || !seedHashLocks) {
      snprintf(ERR_MSG, sizeof(ERR_MSG), "Cannot allocate 4GB for seed hash counters\n");
      return ERR_MSG;
    }
    memset(strCtx, 0, numStrThreads * sizeof(strScanCtx_t));
    printf("Spawning %d threads build STR table...\n", numStrThreads);
    fflush(stdout);
    pthread_mutex_init(&strMut, 0);
    pthread_mutex_lock(&strMut);
    for (i = 0; i < numStrThreads; i++) {
      strCtx[i].mut            = &strMut;
      strCtx[i].refRecs        = refRecs;
      strCtx[i].numRefSeqs     = numRefSeqs;
      strCtx[i].nextRefSeq     = &nextRefSeq;
      strCtx[i].seedHashCounts = seedHashCounts;
      strCtx[i].seedHashLocks  = seedHashLocks;
      strCtx[i].cfghdr         = cfghdr;
      pthread_create(&strThreads[i], NULL, strScanThread, &strCtx[i]);
    }
    pthread_mutex_unlock(&strMut);
  }

  // Allocate reference sequence space according to input length
  refSeqAlloc  = totRefLen / 4 + 64;
  refMaskAlloc = refSeqAlloc / 2;
  refCodeAlloc = cfghdr->maxMultBaseSeeds > 1 ? refSeqAlloc * 2 : 0;
  if (!(refSeq = malloc(refSeqAlloc))) {
    snprintf(ERR_MSG, sizeof(ERR_MSG), "Cannot allocate %llu bytes for reference sequence\n", refSeqAlloc);
    return ERR_MSG;
  }
  if (!(refMask = malloc(refMaskAlloc))) {
    snprintf(ERR_MSG, sizeof(ERR_MSG), "Cannot allocate %llu bytes for reference mask\n", refMaskAlloc);
    return ERR_MSG;
  }
  if (refCodeAlloc && !(refCode = malloc(refCodeAlloc))) {
    snprintf(
        ERR_MSG,
        sizeof(ERR_MSG),
        "Cannot allocate %llu bytes for reference multi-base codes\n",
        refCodeAlloc);
    return ERR_MSG;
  }

  // Open reference sequence output file
  if (config->refOutput && !(outRefFile = fopen(config->refOutput, "wb"))) {
    snprintf(ERR_MSG, sizeof(ERR_MSG), "Cannot open reference output file %s\n", config->refOutput);
    return ERR_MSG;
  }
  // Open reference index output file
  if (config->refIdxFname && !(refIdxFile = fopen(config->refIdxFname, "wb"))) {
    snprintf(ERR_MSG, sizeof(ERR_MSG), "Cannot open reference index file %s\n", config->refIdxFname);
    return ERR_MSG;
  }
  // Open the reference repeat mask output file
  if (config->repMaskFname && !(repMaskFile = fopen(config->repMaskFname, "wb"))) {
    snprintf(ERR_MSG, sizeof(ERR_MSG), "Cannot open reference repeat mask file %s\n", config->repMaskFname);
    return ERR_MSG;
  }

  // Process reference sequences
  printf("Encoding binary reference sequence...\n");
  fflush(stdout);
  char     outRefBuf[1024];
  uint8_t  repMaskBuf[1024];
  uint32_t refIdxBuf[256];
  int      outRefCount  = 0;
  int      repMaskCount = 0;

  config->refSeq  = (hashTableSeq_t**)malloc(numRefSeqs * sizeof(hashTableSeq_t*));
  config->seqName = (char**)malloc(numRefSeqs * sizeof(char*));

  // Process primary assembly sequences, then alternate contigs
  int     alt, ref, realSeq, shift;
  int32_t seqNum = 0, numSeq = 0;
  //uint32_t seqLen;
  uint32_t cpyLen, blkLen, srcPos;
  uint8_t *srcPtr, *dstPtr;
  uint64_t seqBeg = 0, seqEnd, seqEff, blkTot = 0;
  uint64_t refAltStart;
  char*    seqSrc;
  char*    curCharToCode;
  nextReport = SEQ_REPORT_INTERVAL;
  totSeqLen  = 0;
  for (alt = 0; alt <= 1; alt++) {
    // Between primary and alternate, capture the alternate section's start position
    if (alt == 1) {
      config->hdr->refAltSeed  = ceil(seqBeg / config->hdr->refSeedInterval);
      config->hdr->refAltStart = refAltStart = seqBeg;
      // Allocate a buffer for liftover match flags, only for the alt sequences remaining
      altMatchesAlloc = refMaskAlloc - refAltStart / 8;
      if (!(altMatches = calloc(altMatchesAlloc, 1))) {
        snprintf(
            ERR_MSG,
            sizeof(ERR_MSG),
            "Cannot allocate %llu bytes for alt contig match flags\n",
            altMatchesAlloc);
        return ERR_MSG;
      }
    }
    // Extra 'ref' iteration before alt=0 and after alt=1
    for (ref = alt - 1; ref < numRefSeqs + alt; ref++) {
      // Real reference sequences
      realSeq = (ref >= 0 && ref < numRefSeqs);
      done    = (ref == numRefSeqs);
      if (realSeq) {
        if (refRecs[ref].isAlt != alt) continue;
        seqSrc = refRecs[ref].seq;
        //seqLen = refRecs[ref].seqLen;
        blkLen               = refRecs[ref].blockLen;
        cpyLen               = refRecs[ref].trimLen;
        seqEnd               = seqBeg + blkLen;
        seqEff               = seqBeg - refRecs[ref].begTrim;
        curCharToCode        = refRecs[ref].charToCode;
        const size_t namelen = strlen(refRecs[ref].name) + strlen(seqnameSuffix) + 1;
        config->seqName[ref] = (char*)malloc(namelen);
        snprintf(config->seqName[ref], namelen, "%s%s", refRecs[ref].name, seqnameSuffix);
        config->refSeq[ref]           = (hashTableSeq_t*)malloc(sizeof(hashTableSeq_t));
        config->refSeq[ref]->seqStart = seqBeg;
        config->refSeq[ref]->seqLen   = refRecs[ref].seqLen;
        config->refSeq[ref]->begTrim  = refRecs[ref].begTrim;
        config->refSeq[ref]->endTrim  = refRecs[ref].endTrim;
        refRecs[ref].seqStart         = seqBeg;  // for modifying reference buffers with pop snps
        if (ref < origRefSeqs) {
          seqNum = ref;
          totSeqLen += refRecs[ref].seqLen;
        } else {
          seqNum = -1;
          if (popAltSeqs) totSeqLen += refRecs[ref].seqLen;
        }
        // In alt contigs, copy the mask of positions matching the liftover sequence
        if (alt && refRecs[ref].liftMatch) {
          srcPtr = refRecs[ref].liftMatch;
          dstPtr = &altMatches[(refSeqLen - refAltStart) >> 3];
          shift  = (refSeqLen - refAltStart) & 7;
          // Copy 7 bytes at a time, shifting to arbitrary bit alignment
          for (srcPos = 0; srcPos < cpyLen; srcPos += 7 * 8, srcPtr += 7, dstPtr += 7) {
            *(uint64_t*)dstPtr |= *(uint64_t*)srcPtr << shift;
          }
        }
      }
      // Begin/end padding
      else {
        seqNum = ref;
        seqSrc = nullPaddingBlock;
        //seqLen = 0;
        blkLen = REF_SEQ_END_PAD_BASES;
        cpyLen = 0;
        seqEnd = seqBeg;
        seqEff = seqBeg;
        // Default conversion tables for converting the null padding.
        curCharToCode = BASE_ENCODING;
      }

      // Fill a buffer with copies of a 16-byte sequence descriptor word
      refIdxBuf[0] = seqNum;        // 32 bits: Sequence ID
      refIdxBuf[1] = seqBeg / 256;  // 32 bits: Sequence start coordinate / 256
      refIdxBuf[2] = seqEnd / 256;  // 32 bits: Sequence end coordinate / 256
      refIdxBuf[3] = seqEff / 256;  // 32 bits: Trim-adjusted effective start coordinate / 256
      for (i = 4; i < sizeof(refIdxBuf) / 4; i += 4) memcpy(refIdxBuf + i, refIdxBuf, 16);
      // Write reference index section for this sequence - one 16-byte word for each 1024 bases with padding
      blkTot += blkLen;
      uint64_t idxTot   = blkTot / 1024 * 16;
      int32_t  idxBytes = idxTot - refIdxLen;
      refIdxLen += idxBytes;
      // Extend after end to 256-byte boundary
      if (done) {
        int32_t idxPad =
            (DRAM_FILE_ALIGN_BYTES - (refIdxLen % DRAM_FILE_ALIGN_BYTES)) % DRAM_FILE_ALIGN_BYTES;
        idxBytes += idxPad;
        refIdxLen += idxPad;
      }
      // Write copies of index word to the reference index file
      int32_t idxLeft = idxBytes;
      while (idxLeft > 0) {
        int32_t idxChunk = (idxLeft > sizeof(refIdxBuf) ? sizeof(refIdxBuf) : idxLeft);
        fwrite(refIdxBuf, 1, idxChunk, refIdxFile);
        idxLeft -= idxChunk;
      }

#ifndef LOCAL_BUILD
      // Update ref_index.bin digest
      for (idxLeft = idxBytes; idxLeft > 0; idxLeft -= 16) {
        // This may look wrong - we hash the same buffer repeatedly - but we are correctly
        // hashing the same data we just wrote, so many copies of the first 16 bytes in refIdxBuf
        uint64_t* qp = (uint64_t*)refIdxBuf;

        refIndexDigest = crc32c_hw(refIndexDigest, (unsigned char*)qp, 8);
        qp++;
        refIndexDigest = crc32c_hw(refIndexDigest, (unsigned char*)qp, 8);
      }
#endif

      char* p      = seqSrc;
      char* cpyEnd = seqSrc + cpyLen;
      char* blkEnd = seqSrc + blkLen;
      for (; p < blkEnd; p++) {
        if (p < cpyEnd) {
          if (*p >= 'a') {
            // Lowercase base, repeat region detected, set the bit
            repMaskByte |= 1 << (7 - repMaskBits);
          }
        }
        repMaskBits++;
        c    = p < cpyEnd ? *p : '~';
        code = curCharToCode[c];
        base = codeToBase[code];
        seqByte |= ((base & 3) << ((refSeqLen & 3) << 1));
        maskByte |= ((base >> 2) << (refSeqLen & 7));
        refByte |= code << ((refSeqLen & 1) << 2);
        refCodeHist[code]++;
        if ((refSeqLen & 3) == 3) {
          refSeq[refSeqLen >> 2] = seqByte;
          seqByte                = 0;
        }
        if ((refSeqLen & 7) == 7) {
          refMask[refSeqLen >> 3] = maskByte;
          maskByte                = 0;
        }
        if ((refSeqLen & 1) == 1) {
          if (refCode != NULL) {
            refCode[refSeqLen >> 1] = refByte;
          }
          if (outRefFile) {
            outRefBuf[outRefCount++] = refByte;
#ifndef LOCAL_BUILD
            // Update reference.bin digest
            refDigest = crc32c_hw(refDigest, &refByte, 1);
#endif
            if (outRefCount == sizeof(outRefBuf)) {
              fwrite(outRefBuf, 1, outRefCount, outRefFile);
              outRefCount = 0;
            }
          }
          refByte = 0;
        }
        // For repeat mask file, encoded as 1 bit per base !
        if (repMaskFile) {
          if (repMaskBits == 8) {
            repMaskBuf[repMaskCount++] = repMaskByte;
            repMaskByte                = 0;
            repMaskBits                = 0;
          }
          if (repMaskCount == sizeof(repMaskBuf)) {
            fwrite(repMaskBuf, 1, repMaskCount, repMaskFile);
            repMaskCount = 0;
          }
        }
        refSeqLen++;
        if (refSeqLen >> 2 > refSeqAlloc) {
          snprintf(ERR_MSG, sizeof(ERR_MSG), "Reference sequence longer than allocated buffer\n");
          printf("refSeqLen: %llu\n", refSeqLen);
          printf("refSeqAlloc: %llu\n", refSeqAlloc);
          return ERR_MSG;
        }
      }
      seqBeg = refSeqLen;
      if (realSeq) {
        if (!config->strFname) freeRefRec(&refRecs[ref]);
        numSeq += (ref < origRefSeqs);
      }
      if (totSeqLen >= nextReport || done) {
        printf(
            "  %d sequence%s, %lld bases (%lld after trimming/padding)\n",
            numSeq,
            numSeq == 1 ? "" : "s",
            totSeqLen,
            refSeqLen);
        fflush(stdout);
        nextReport += SEQ_REPORT_INTERVAL;
      }
    }
#ifndef LOCAL_BUILD
    WatchDogCheckin(watchdogID);
#endif
  }

  // Wait for STR thread completion
  uint64_t extendedSeeds    = 0;
  uint64_t validSeeds       = 0;
  uint64_t nonExtendedSeeds = ceil(adjBasesNotN / cfghdr->refSeedInterval);
  if (useStrThreads) {
    for (i = 0; i < numStrThreads; i++) {
      pthread_join(strThreads[i], NULL);
      extendedSeeds += strCtx[i].extendedSeeds;
      validSeeds += strCtx[i].validSeeds;
    }
    if (allowSeedExt) {
      printf("Estimated extended seeds: %llu / %llu\n", extendedSeeds, validSeeds);
      nonExtendedSeeds = validSeeds - extendedSeeds;
    } else {
      extendedSeeds    = 0;
      nonExtendedSeeds = validSeeds;
    }
    free(seedHashCounts);
    free(seedHashLocks);
    seedHashCounts = NULL;
    seedHashLocks  = NULL;
    pthread_mutex_destroy(&strMut);
  }
  // Write STR table
  if (config->strFname) {
    FILE*   strFile = fopen(config->strFname, "wb");
    int32_t strTot  = 0;
    if (!strFile) {
      snprintf(ERR_MSG, sizeof(ERR_MSG), "Cannot open STR file %s\n", config->strFname);
      return ERR_MSG;
    }
    for (ref = 0; ref < origRefSeqs; ref++) {
      refRec_t* r = &refRecs[ref];
      if (!r->str) {
        snprintf(ERR_MSG, sizeof(ERR_MSG), "STR alloc failure\n");
        return ERR_MSG;
      }
      if (r->strNum && fwrite(r->str, sizeof(htStrRec_t), r->strNum, strFile) != r->strNum) {
        snprintf(ERR_MSG, sizeof(ERR_MSG), "Write failure to STR file %s\n", config->strFname);
        return ERR_MSG;
      }
      strTot += r->strNum;
    }
    fclose(strFile);
    printf(
        "Wrote STR table to '%s': %d records, %lld bytes\n",
        config->strFname,
        strTot,
        (int64_t)strTot * sizeof(htStrRec_t));
    fflush(stdout);
  }

  cfghdr->refSeqLen  = totRefLen;
  cfghdr->refLenRaw  = totSeqLen;
  cfghdr->refLenNotN = refBasesNotN;

#ifndef LOCAL_BUILD
  WatchDogDeregister(watchdogID);
#endif

  // Pad to keep memory accesses safe
  memset(&refSeq[refSeqLen >> 2], 0, refSeqAlloc - (refSeqLen >> 2));
  memset(&refMask[refSeqLen >> 3], 0, refMaskAlloc - (refSeqLen >> 3));
  if (refCode != NULL) {
    memset(&refCode[refSeqLen >> 1], 0, refCodeAlloc - (refSeqLen >> 1));
  }
  if (outRefFile) {
    if (outRefCount) fwrite(outRefBuf, 1, outRefCount, outRefFile);
    fclose(outRefFile);
    printf("Wrote binary reference to '%s': %lu bytes\n", config->refOutput, refSeqLen >> 1);
    fflush(stdout);
  }
  if (config->refIdxFname) {
    fclose(refIdxFile);
    printf("Wrote reference index to '%s': %lu bytes\n", config->refIdxFname, refIdxLen);
    fflush(stdout);
  }
  if (config->repMaskFname) {
    if (repMaskBits) {
      repMaskBuf[repMaskCount++] = repMaskByte;
    }
    if (repMaskCount) fwrite(repMaskBuf, 1, repMaskCount, repMaskFile);
    fclose(repMaskFile);
    printf("Wrote repeat mask bitmap to '%s': %lu bytes\n", config->repMaskFname, refSeqLen >> 3);
    fflush(stdout);
  }

  char outSnpsBuf[1024];
  int  outSnpsCount = 0;



  free(sortedRefRecs);

  // Size extension table according to the estimated number of extended seeds found by pre-hashing,
  // but override with configured VALUE if specified, and always zero if no seed extension allowed
  int64_t extTabRecs, extTabReserve;
  if (!allowSeedExt)
    extTabRecs = 0;
  else if (config->extTableAlloc)
    extTabRecs = config->extTableAlloc;
  else if (bytes) {
    extTabRecs = (memSize - bytes - totOtherSize) / 8;
    extTabRecs -= (extTabRecs % 1024);
  } else {
    // The estimate should be pretty tight, even a little high, but add a little for extra liftover hits
    extTabRecs = extendedSeeds * 1.01 + 8192;
  }
  if (extTabRecs % 1024) extTabRecs += (1024 - (extTabRecs % 1024));
  extTabReserve = extTabRecs * sizeof(extend_hit_t);

  // If not specified, automatic hash table size based on input length and memory limit
  if (!bytes) {
    int64_t maxTableSize = memSize - totOtherSize - extTabReserve;
    // Enforce 32GB cap on hash table proper, so max of 32 construction chunks
    if (maxTableSize >> 35) maxTableSize = 1ULL << 35;
    if (cfgMemSize || memSize == defaultSize) {
      printf("Memory limit:       %s\n", bytesReadable(memSize));
      printf("Reference reserved: %s\n", bytesReadable(refSize));
      printf("Seq index reserved: %s\n", bytesReadable(idxSize));
      printf("SJ table reserved:  %s\n", bytesReadable(sjReserve));
      printf("Ext table reserved: %s\n", bytesReadable(extTabReserve));
      if (maxTableSize < 0) {
        printf("Remaining space:    -%s\n", bytesReadable(-maxTableSize));
        fflush(stdout);
        snprintf(ERR_MSG, sizeof(ERR_MSG), "Negative room available for hash table within memory limit\n");
        return ERR_MSG;
      }
      printf("Remaining space:    %s\n", bytesReadable(maxTableSize));
    }
    // Pick the largest number of address bits where a minimum 33/64 squeeze factor will fit
    int sixtyfourths;
    for (addrBits = 0; (1ull << (addrBits + 1)) * 33 / 64 <= maxTableSize; addrBits++)
      ;
    // Pick the largest squeeze factor 33/64 .. 64/64 that will fit
    for (sixtyfourths = 64; (1ull << addrBits) * sixtyfourths / 64 > maxTableSize; sixtyfourths--)
      ;
    if (addrBits < MAX_WRAP_BYTES_LOG2) {
      addrBits     = MAX_WRAP_BYTES_LOG2;
      sixtyfourths = 64;
    }
    cfghdr->tableAddrBits  = addrBits;
    cfghdr->tableSize64ths = sixtyfourths;
  }

  // Recompute and report planned table size
  tableBytes             = ((uint64_t)1 << cfghdr->tableAddrBits) * cfghdr->tableSize64ths / 64;
  cfghdr->hashTableBytes = tableBytes;
  printf("Hash table size:    %s\n", bytesReadable(tableBytes));
  if (tableBytes < MAX_WRAP_BYTES) {
    snprintf(
        ERR_MSG,
        sizeof(ERR_MSG),
        "Hash table size %llu smaller than minimum supported %d KB\n",
        tableBytes,
        (1 << (MAX_WRAP_BYTES_LOG2 - 10)));
    return ERR_MSG;
  }
  // Adjust extension table allocation to take up any slack
  extTabReserve = memSize - tableBytes - totOtherSize;
  extTabReserve -= extTabReserve % 256;
  extTabRecs = extTabReserve / sizeof(extend_hit_t);
  if (extTabRecs < extendedSeeds && !config->extTableAlloc) {
    snprintf(
        ERR_MSG,
        sizeof(ERR_MSG),
        "Space for %lld extension table recs less than estimated %lld needed\n",
        extTabRecs,
        extendedSeeds);
    return ERR_MSG;
  }
  config->extTableAlloc = extTabRecs;
  printf("Ext table space:    %s\n\n", bytesReadable(extTabReserve));

  if (tableBytes / HASH_RECORD_BYTES < nonExtendedSeeds * 0.75 && !config->overrideCheck) {
    snprintf(ERR_MSG, sizeof(ERR_MSG), "The requested hash table size looks much too small!\n");
    return ERR_MSG;
  }

  int      tableChunks  = (addrBits < INDEPENDENT_ADDR_BITS ? 1 : 1 << (addrBits - INDEPENDENT_ADDR_BITS));
  uint64_t chunkBytes   = tableBytes / tableChunks;
  int64_t  popBufBytes  = config->writeCompFile ? totRefLen * 5 : 0;
  int      chunksInMem  = tableChunks > config->maxGB ? config->maxGB : tableChunks;
  int      chunksActive = chunksInMem > config->maxThreads ? config->maxThreads : chunksInMem;
  int      estGB =
      (int)((1.6 * chunkBytes * chunksInMem + 1.0 * chunkBytes * chunksActive + popBufBytes + totRefLen * 3 / 8) / (1 << 30)) +
      1;
  printf("Max worker threads: %d\n", config->maxThreads);
  printf("Thread chunk size:  %s\n", bytesReadable(chunkBytes));
  printf("Max chunks in mem:  %d\n", config->maxGB);
  printf("Est peak memory:    %d GB\n", estGB);
#ifdef LOCAL_BUILD
  printf("Write uncompressed: %s\n", config->writeHashFile ? "true" : "false");
  printf("Write compressed:   %s\n", config->writeCompFile ? "true" : "false");
#endif

  // Compute CRC lengths
  cfghdr->priCrcBits = cfghdr->tableAddrBits + HASH_BYTE_ADDR_START;
  cfghdr->secCrcBits = cfghdr->priCrcBits > MAX_SEC_CRC_BITS ? MAX_SEC_CRC_BITS : cfghdr->priCrcBits;

  // Cap seed lengths so hash bits not covered by the address fit in the hash records
  unsigned int maxPriSeed = cfghdr->priCrcBits / 2;
  if (cfghdr->priSeedBases > maxPriSeed) {
    printf(
        "\nBase seed length of %u bases too long for %s hash table.\n"
        "  The limit increases by 1 base for every 2 table address bits.\n"
        "  Using %u bases, the maximum seed length for this table size.\n",
        cfghdr->priSeedBases,
        bytesReadable(tableBytes),
        maxPriSeed);
    cfghdr->priSeedBases = maxPriSeed;
  }
  int initialMaxSeedBases = cfghdr->maxSeedBases;
  APPLY_MIN(cfghdr->maxSeedBases, cfghdr->priSeedBases);
  int reqSeedExt = cfghdr->maxSeedBases - cfghdr->priSeedBases;
  reqSeedExt -= (reqSeedExt & 1);
  if (reqSeedExt > MAX_NET_SEED_EXTENSION) {
    if (initialMaxSeedBases != 999)
      printf(
          "\nMaximum extended seed of %u bases is too long for current settings.\n"
          "  Limit is %d longer than the primary seed.  Using %u bases.\n",
          cfghdr->maxSeedBases,
          MAX_NET_SEED_EXTENSION,
          cfghdr->priSeedBases + MAX_NET_SEED_EXTENSION);
    reqSeedExt = MAX_NET_SEED_EXTENSION;
  }
  cfghdr->maxSeedBases = cfghdr->priSeedBases + reqSeedExt;
  APPLY_MINMAX(cfghdr->maxSeedFreqLen, cfghdr->priSeedBases, cfghdr->maxSeedBases);

  // Enforce minimum anchored seed bin size
  if (cfghdr->anchorBinBits && cfghdr->anchorBinBits < MIN_ANCHOR_BIN_BITS) {
    printf("\nSetting anchor-bin-bits to minimum %d\n", MIN_ANCHOR_BIN_BITS);
    cfghdr->anchorBinBits = MIN_ANCHOR_BIN_BITS;
  }
  // Check if the seed looks too short to get any mapping uniqueness
  if ((1ULL << (2 * cfghdr->priSeedBases)) < adjBasesNotN &&
      (!cfghdr->anchorBinBits || (cfghdr->priSeedBases << 1) < cfghdr->anchorBinBits) &&
      !config->overrideCheck) {
    if (cfghdr->anchorBinBits)
      snprintf(
          ERR_MSG,
          sizeof(ERR_MSG),
          "The requested seed length looks too short for unique mapping in this anchor bin size!\n");
    else
      snprintf(
          ERR_MSG,
          sizeof(ERR_MSG),
          "The requested seed length looks too short for unique mapping in this reference!\n");
    return ERR_MSG;
  }
  // Verify seed length and anchor bin size are compatible
  if (cfghdr->anchorBinBits && (cfghdr->priSeedBases << 1) > KEY_ANCHOR_OFFSET + cfghdr->anchorBinBits) {
    snprintf(
        ERR_MSG,
        sizeof(ERR_MSG),
        "For anchor-bin-bits=%d, maximum seed length is %d\n",
        cfghdr->anchorBinBits,
        (KEY_ANCHOR_OFFSET + cfghdr->anchorBinBits) / 2);
    return ERR_MSG;
  }

  // Fetch appropriate polynomials
  memcpy(cfghdr->priCrcPoly, CRC_POLYS[cfghdr->priCrcBits][config->priPolyIndex], 8);
  memcpy(cfghdr->secCrcPoly, CRC_POLYS[cfghdr->secCrcBits][config->secPolyIndex], 8);

  printf("\nSettings:\n");
  printf("  Initial seed bases                 : %u\n", cfghdr->priSeedBases);
  printf("  Maximum seed bases                 : %u\n", cfghdr->maxSeedBases);
  printf("  Maximum seed extension increment   : %u\n", cfghdr->maxExtIncrement);
  printf("  Reference seed interval            : %g\n", cfghdr->refSeedInterval);
  printf("  Hash table address bits            : %u\n", cfghdr->tableAddrBits);
  printf("  Hash table size 64ths              : %u\n", cfghdr->tableSize64ths);
  printf("  Maximum seed frequency             : %u\n", cfghdr->maxSeedFreq);
  printf("  Maximum primary seed frequency     : %u\n", cfghdr->priMaxSeedFreq);
  printf("  Ext seed length for max frequency  : %u\n", cfghdr->maxSeedFreqLen);
  printf("  Target seed frequency              : %g\n", cfghdr->targetSeedFreq);
  if (cfghdr->maxMultBaseSeeds > 0) {
    printf("  Maximum seeds for multi-base codes : %u\n", cfghdr->maxMultBaseSeeds);
  }
  if (cfghdr->liftMatchSeedInt > 0) {
    printf("  Liftover-matching seeds interval   : %u\n", cfghdr->liftMatchSeedInt);
  }
  printf("  Thinning frequency cap             : %g\n", cfghdr->thinningFreqCap);
  printf("  Max thinning factor                : %u\n", cfghdr->thinningPeriod);
  printf("  Primary hash CRC length            : %u\n", cfghdr->priCrcBits);
  printf("  Secondary hash CRC length          : %u\n", cfghdr->secCrcBits);
  printf("  Cost coefficient of seed length    : %g\n", cfghdr->seedLenCost);
  printf("  Cost coefficient of seed frequency : %g\n", cfghdr->seedFreqCost);
  printf("  Cost penalty for seed extension    : %g\n", cfghdr->extensionCost);
  printf("  Cost penalty for extension steps   : %g\n", cfghdr->extStepCost);
  printf("  Cost penalty for extension records : %g\n", cfghdr->extRecCost);
  printf("  Seed extension repair strategy     : %u\n", cfghdr->repairStrategy);
  printf("  Min probability for seed repair    : %g\n", cfghdr->minRepairProb);
  printf("  Log2 bin for anchored seed search  : %u\n", cfghdr->anchorBinBits);
  printf("  Random sample with HIFREQ seed     : %u\n", cfghdr->hiFreqRandHit);
  printf("  Extend freq for random sample      : %u\n", cfghdr->extRandHitFreq);
  printf("  Primary hash CRC polynomial        : %llX\n", ((long long int*)cfghdr->priCrcPoly)[0]);
  printf("  Secondary hash CRC polynomial      : %llX\n", ((long long int*)cfghdr->secCrcPoly)[0]);
  printf("  Base conversion for methylation    : ");
  if (config->methylatedConv == HT_TYPE_METHYL_G_TO_A) {
    printf("G to A\n");
  } else if (config->methylatedConv == HT_TYPE_METHYL_C_TO_T) {
    printf("C to T\n");
  } else if (config->methylatedConv == HT_TYPE_METHYL_COMBINED) {
    printf("G to A and C to T (combined)\n");
  } else {
    printf("disabled\n");
  }
  printf("\n");

  if (config->showIntParams) {
    printInternalParams(config);
  }

  if (config->testOnly) return 0;

  setCmdLine(config, argc, argv);

  // Quit if not writing compressed or uncompressed hash table
  if (!(config->writeHashFile || config->writeCompFile)) {
    printf("Skipping actual hash table building\n");
    if (config->hashFname) printf("  '%s' not modified\n", config->hashFname);
    if (config->compFname) printf("  '%s' not modified\n", config->compFname);
    if (config->extTabFname) printf("  '%s' not modified\n", config->extTabFname);
    if (config->configFname) printf("  '%s' not modified\n", config->configFname);
    if (config->configBinFname) printf("  '%s' not modified\n", config->configBinFname);
    if (config->statsFname) printf("  '%s' not modified\n", config->statsFname);
    goto generateHashTableCleanup;
  }


  // Some memory we can free before actually building the hash table
  if (refRecs)
    for (i = 0; i < numRefSeqs; i++) {
      // In combined methyl, pairs of refRecs share the same seq, so need to avoid freeing twice.
      if (config->methylatedConv == HT_TYPE_METHYL_COMBINED && (i % 2)) refRecs[i].seq = NULL;
      freeRefRec(&refRecs[i]);
    }
  free(refRecs);
  refRecs = NULL;

  free(inpSeq);
  inpSeq = NULL;

  char* errStr = buildHashTable(
      config,
      refSeq,
      refCodeHist,
      refMask,
      refCode,
      altMatches

  );
  if (errStr) {
    return errStr;
  }

  hashDigest   = config->hdr->hashDigest;
  extTabDigest = config->hdr->extTabDigest;
  printf("reference.bin digest:    0x%08lX\n", refDigest);
  printf("ref_index.bin digest:    0x%08lX\n", refIndexDigest);
  printf("hash_table.bin digest:   0x%08lX\n", hashDigest);
  printf("extend_table.bin digest: 0x%08lX\n", extTabDigest);
  if (liftoverDigest) {
    printf("ALT liftover digest:     0x%08lX\n", liftoverDigest);
  }
  if (config->maskBedDigest) {
    printf("Mask bed digest:         0x%08lX\n", config->maskBedDigest);
  }
  if (popSnpsDigest) {
    printf("ref_pop_snps.bin digest: 0x%08lX\n", popSnpsDigest);
  }
  config->hdr->digest = refDigest + refIndexDigest + hashDigest + extTabDigest + popSnpsDigest;
  // Don't allow zero digest, because then SW would think the reference
  // is already loaded right after reset
  if (config->hdr->digest == 0) config->hdr->digest = 1;
  printf("Table digest: 0x%08lX\n", config->hdr->digest);

  config->hdr->numRefSeqs     = origRefSeqs;
  config->hdr->digestType     = DIGEST_CRC32C;
  config->hdr->refDigest      = refDigest;
  config->hdr->refIndexDigest = refIndexDigest;
  config->hdr->hashDigest     = hashDigest;
  config->hdr->liftoverDigest = liftoverDigest;
  config->hdr->popSnpsDigest  = popSnpsDigest;
  config->usedReadBuf         = 0;
  writeHashCfgBin(config->configBinFname, config);
  writeHashCfgTxt(config->configFname, config);

  printf("Wrote configuration to '%s'\n", config->configFname);

generateHashTableCleanup:

  free(inpSeq);
  if (refRecs)
    for (i = 0; i < numRefSeqs; i++) {
      // In combined methyl, pairs of refRecs share the same seq, so need to avoid freeing twice.
      if (config->methylatedConv == HT_TYPE_METHYL_COMBINED && (i % 2)) refRecs[i].seq = NULL;
      freeRefRec(&refRecs[i]);
    }
  free(refRecs);
  free(refSeq);
  free(refMask);
  free(refCode);
  free(altMatches);
  return 0;
}
