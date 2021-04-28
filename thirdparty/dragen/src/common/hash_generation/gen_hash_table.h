// Copyright 2013-2018 Edico Genome Corporation. All rights reserved.
//
// This file contains confidential and proprietary information of the Edico Genome
// Corporation and is protected under the U.S. and international copyright and other
// intellectual property laws.
//
// $Id$
// $Author$:
// $Change$
// $DateTime$
//
//

#ifndef _gen_hash_table_h_
#define _gen_hash_table_h_

#ifdef __cplusplus
extern "C" {
#endif

#include <inttypes.h>

// STR descriptor record used in hashtable builder: 16 bytes
typedef struct {
  uint32_t mask;   // Mask for decimation; currently just a separate counter for each (refId,period,repeatLen)
  uint32_t refId;  // Reference sequence index
  uint32_t begPos;  // Start position in reference sequence, 0-based
  uint16_t length;  // Number of bases in STR, at least period*repeatLen, but less than period*(repeatLen+1)
  uint8_t  period;  // Length of repeating unit
  uint8_t  repeatLen;  // Number of full units repeated, including the first
} htStrRec_t;

// Local record for reference sequence info
typedef struct {
  char*       name;
  char*       seq;
  htStrRec_t* str;
  uint8_t*    lifted;
  uint8_t*    liftMatch;
  uint64_t    seqStart;
  uint32_t    seqLen;
  uint32_t    begTrim;
  uint32_t    endTrim;
  uint32_t    trimLen;
  uint32_t    endPad;
  uint32_t    blockLen;
  uint32_t    strNum;
  int32_t     seqNum;
  int         isPri;
  int         isAlt;
  int         isDecoy;
  int         isPopAlt;
  int*        charToBase;  // Conversion arrays that specify how to interpret bases from seq.
  char*       charToCode;
} refRec_t;

// -------------------------------------------------------------------------------
// Comparator for refRec_t pointers. Calls strcmp on the name field of the
// refRec_t's that #a# and #b# point to.
int refRecCompareName(const void* a, const void* b);

// IMPORTANT: the version should only be incremented when the hash table is not
// backwords compatible.  This is used to tell user's when they need to regenerate
// hash tables and by QA to load the correct hash table in their automation framework.
#define HT_VERSION 8
static const int HASH_HEADER_SIZE = 512;

// Padding for reference.bin (also used by host software to generate virtual "global"
// position)
#define REF_SEQ_ALIGN_BASES 1024
#define REF_SEQ_MIN_PAD_BASES 256
#define REF_SEQ_END_PAD_BASES 163840

// Padding for reference.bin for pop-alt contigs
#define REF_SEQ_ALIGN_BASES_POP_ALT 1
#define REF_SEQ_MIN_PAD_BASES_POP_ALT 64

// Hash table types for generation.
typedef enum {
  HT_TYPE_NORMAL = 0,
  HT_TYPE_METHYL_G_TO_A,
  HT_TYPE_METHYL_C_TO_T,
  HT_TYPE_METHYL_COMBINED,
  HT_TYPE_ANCHORED,
  HT_TYPE_NUM_MAX,
} HashTableType;

// Algorithm used to generate digest
typedef enum {
  DIGEST_CRC32  = 0,
  DIGEST_CRC32C = 1,
} DigestType;

#define DEFAULT_MEM_SIZE_STR "32GB"
#define MAX_MEM_SIZE_STR "64GB"
#define REF_BASES_PER_IDX_RESERVE_BYTE 64
#define DRAM_FILE_ALIGN_BYTES 256

typedef struct {
  uint32_t hashTableVersion;  // Version number of this hash table
  uint64_t hashTableBytes;    // #bytes in references hash table
  uint32_t priSeedBases;      // Initial seed length to store in hash table
  uint32_t maxSeedBases;      // Max extended seed len to store in secondary hash table
  uint32_t maxExtIncrement;   // Maximum bases to extend a seed by in one step
  double   refSeedInterval;   // Number of positions per reference seed
  uint32_t tableAddrBits;     // Ceiling log2 of the hash table size in bytes
  uint32_t tableSize64ths;    // (33-64) Hash table is (2^tableAddrBits)*(tableSize64ths/64) bytes
  uint32_t maxSeedFreq;       // Max allowed freq for a seed match after extension attempts
  uint32_t priMaxSeedFreq;    // Maximum frequency for a primary seed match (0 => use maxSeedFreq)
  uint32_t maxSeedFreqLen;    // Ramp from priMaxSeedFreq reaches maxSeedFreq at this seed length
  double   targetSeedFreq;    // Target seed frequency for seed extension
  double   thinningFreqCap;   // Soft seed frequency cap for thinning
  uint32_t thinningPeriod;    // Maximum decimation factor for seed thinning
  uint32_t priCrcBits;        // Length of CRC polynomial for primary seeds
  uint32_t secCrcBits;        // Length of CRC polynomial for extended seeds
  double   seedLenCost;       // Cost coefficient of extended seed length
  double   seedFreqCost;      // Cost coefficient of extended seed frequency
  double   extensionCost;     // Cost penalty to extend a seed by any number of bases
  double   extStepCost;       // Cost penalty to incrementally extend a seed another step
  uint32_t repairStrategy;    // Seed extension repair: 0=none, 1=best, 2=rand
  double   minRepairProb;     // Minimum probability of success for 'best' seed repair
  uint32_t anchorBinBits;     // Bits defining reference bins for anchored seed search, 0=none
  uint32_t hiFreqRandHit;     // Include a random hit with each HIFREQ record
  uint32_t extRandHitFreq;    // Minimum EXTEND frequency to include a random hit
  uint8_t  priCrcPoly[8];     // CRC polynomial for primary seeds
  uint8_t  secCrcPoly[8];     // CRC polynomial for extended seeds
  uint64_t refSeqLen;         // Number of bases in reference (including padding)
  uint64_t refLenRaw;         // Raw sequence length (no padding)
  uint64_t refLenNotN;        // Raw sequence length minus the the N's
  uint32_t digest;            // Digest of reference, ref_index and hashtables
  uint32_t numRefSeqs;        // Number of reference sequences
  uint32_t digestType;        // Method for computing digest: 0=CRC32 (old) 1=CRC32C (new)
  uint32_t refDigest;         // Digest of reference.bin
  uint32_t refIndexDigest;    // Digest of ref_index.bin
  uint32_t hashDigest;        // Digest of hash_table.bin
  uint32_t liftoverDigest;    // Digest of ALT liftover file (if used)
  uint32_t refAltSeed;        // Threshold seed index beginning alternate contigs
  uint64_t refAltStart;       // Threshold flat reference position beginning alternate contigs
  uint32_t extTabRecs;        // Number of 8-byte records in extend_table.bin
  uint32_t extTabDigest;      // Digest of extend_table.bin
  double   extRecCost;        // Cost penalty for each EXTEND or INTERVAL record
  uint32_t minFreqToExtend;   // Minimum seed frequency eligible for seed extension
  uint32_t maxMultBaseSeeds;  // Maximum number of seeds to populate overlapping multi-base codes
  uint32_t popSnpsDigest;     // Digest of ref_pop_snps.bin
  uint32_t liftMatchSeedInt;  // Interval between seed indexes where alt contigs match liftover
  uint8_t  padding[264];      // Reserved for future use
} __attribute__((packed)) hashTableHeader_t;

typedef struct {
  uint64_t seqStart;  // Sequence start offset in reference.bin (untrimmed portion)
  uint32_t begTrim;   // Portion of seqLen leading 'N's trimmed from sequence in reference.bin
  uint32_t endTrim;   // Portion of seqLen trailing 'N's trimmed from sequence in reference.bin
  uint32_t seqLen;    // Reference sequence len (original length, including trimmed portions)
} __attribute__((packed)) hashTableSeq_t;

// For reading version 4 hash_table.cfg.bin files and reporting an error
typedef struct {
  uint64_t seqStart;  // Sequence start offset in reference.bin (untrimmed portion)
  uint32_t seqLen;    // Reference sequence len
} __attribute__((packed)) hashTableSeqv4_t;

typedef struct {
  hashTableHeader_t* hdr;
  hashTableSeq_t**   refSeq;
  char**             seqName;

  int         maxThreads;      // Maximum worker threads
  int         maxGB;           // Maxiumum ~1GB thread table chunks in memory at once
  int         writeHashFile;   // Boolean - write uncompressed hash_table.bin
  int         writeCompFile;   // Boolean - write compressed hash_table.cmp
  const char* sizeStr;         // Size of hash table, units B|KB|MB|GB
  const char* memSizeStr;      // Memory limit (hash table + reference) units B|KB|MB|GB
  const char* sjSizeStr;       // Space to reserve for RNA annotated SJs, NULL/empty for automatic
  uint32_t    methylatedConv;  // For methlated DNA processing, convert G to A or C to T
  uint32_t    extTableAlloc;   // 8-byte records to reserve in extend_table.bin, 0=automatic
  int         priPolyIndex;    // Index of CRC polynomial for hashing primary seeds
  int         secPolyIndex;    // Index of CRC polynomial for hashing extended seeds

  char*    refInput;            // Name of FASTA input file
  char*    altLiftover;         // Name of SAM format liftover of alternate contigs in refInput
  char*    configFname;         // Name of hash table configuration file (txt)
  char*    configBinFname;      // Name of hash table configuration file (bin)
  char*    hashFname;           // Name of hash table output file
  char*    compFname;           // Name of compressed hash table output file
  char*    extTabFname;         // Name of extension table output file
  char*    refOutput;           // Name of reference output file
  char*    refIdxFname;         // Name of reference index output file
  char*    repMaskFname;        // Name of repeat mask bitmap output file
  char*    strFname;            // Name of short tandem repeats file
  char*    statsFname;          // Name of hash table stats output file
  char*    decoyFname;          // Name of additional decoy FASTA file
  char*    maskBed;             // Name of BED file for masking
  uint32_t maskBedDigest;       // Digest of the BED file for masking
  char*    hostVersion;         // Software version string
  char*    cmdLine;             // Command line used to generate hash table
  int      overrideCheck;       // Override hash tables size check
  int      testOnly;            // Testing - show user parameters, but don't do anything
  int      showIntParams;       // Testing - show internal parameters
  uint8_t* readBuf;             // Buffer used when reading binary config file
  int      usedReadBuf;         // 1 if readBuf used, 0 otherwise
  int      altContigValidate;   // If false, skips hg38/hg19 alt-contigs validation (dragen only)
  int      autoDetectValidate;  // If false, skips auto-detect of input fasta for liftover file (dragen only)

  const char* autoDetectDir;  // Which directory the auto-detect should scan for autodetect files. Used for
                              // testing (dragen only).
  char* popAltContigsFname;   // Name of population based alternate contigs FASTA input file
  char* popAltLiftoverFname;  // Name of population based SAM format liftover of alternate contigs
  char* popSnpsInput;         // Name of population based SNPs VCF input file
  char* popSnpsOutput;        // Name of population based SNPs binary output file
} hashTableConfig_t;

void setDefaultHashParams(hashTableConfig_t* defConfig, const char* dir, HashTableType hashTableType);

char* generateHashTable(hashTableConfig_t* config, int argc, char* argv[]);

void freeHashParams(hashTableConfig_t* config);

#ifndef LOCAL_BUILD
void AltContigValidate(hashTableConfig_t* config, const refRec_t* refRecs, int32_t numRefSeqs);
#endif

#ifdef __cplusplus
}
#endif

#endif
