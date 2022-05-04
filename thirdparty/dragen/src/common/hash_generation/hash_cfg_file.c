// Copyright 2013-2017 Edico Genome Corporation. All rights reserved.
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

#define __STDC_FORMAT_MACROS  // for PRIu64

#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <inttypes.h>
#include <libgen.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "gen_hash_table.h"

// Disable the following compiler warnings
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#pragma GCC diagnostic ignored "-Wformat"

#define BUF_SIZE 1 << 20
#define MAX_LINE 2048

char    errorMsg[2048];
uint8_t buf[BUF_SIZE];
size_t  bytesWritten = 0;

//-------------------------------------------------------------------------------swhitmore
// getHashCfgError - Return error message.
//
const char* getHashCfgError()
{
  return &errorMsg[0];
}

//-------------------------------------------------------------------------------swhitmore
// flushBufToDisk - Save buffer contents to disk.
//
void flushBufToDisk(FILE* fp)
{
  fwrite(&buf, bytesWritten, 1, fp);
  bytesWritten = 0;
}

//-------------------------------------------------------------------------------swhitmore
// fillBuf - Copy the given data to the buffer and flush the buffer to disk if it is full.
//
void fillBuf(const void* data, size_t len, FILE* fp)
{
  if (bytesWritten + len >= BUF_SIZE) {
    flushBufToDisk(fp);
  }

  memcpy(buf + bytesWritten, data, len);
  bytesWritten += len;
}

//-------------------------------------------------------------------------------swhitmore
//
void writeHashCfgString(FILE* fp, const char* str)
{
  fillBuf(str, strlen(str) + 1, fp);
}

//-------------------------------------------------------------------------------swhitmore
// writeHashCfgBin - write binary hash table config file.
//
void writeHashCfgBin(const char* filename, const hashTableConfig_t* config)
{
  assert(sizeof(hashTableHeader_t) == HASH_HEADER_SIZE);

  // Set umask so everyone can read the file
  mode_t origmask = umask(002);

  FILE* fp     = fopen(filename, "wb");
  bytesWritten = 0;
  fillBuf(config->hdr, sizeof(hashTableHeader_t), fp);

  size_t i;
  for (i = 0; i < config->hdr->numRefSeqs; ++i) {
    fillBuf(config->refSeq[i], sizeof(hashTableSeq_t), fp);
  }
  for (i = 0; i < config->hdr->numRefSeqs; ++i) {
    writeHashCfgString(fp, config->seqName[i]);
  }
  writeHashCfgString(fp, config->hostVersion);
  writeHashCfgString(fp, config->cmdLine);
  writeHashCfgString(fp, config->refInput);
  writeHashCfgString(fp, config->refOutput);
  writeHashCfgString(fp, config->refIdxFname);
  writeHashCfgString(fp, config->hashFname);
  writeHashCfgString(fp, config->extTabFname);

  if (config->altLiftover) {
    writeHashCfgString(fp, config->altLiftover);
  } else {
    writeHashCfgString(fp, "");
  }

  // Everything above here will be assumed to exist when reading hash_table.cfg.bin.
  // The data below may not exist if the file was generated with older code, but
  // will all exist if and only if there is still more content in the file after reading
  // the last string above.


  if (config->maskBed) {
    writeHashCfgString(fp, config->maskBed);
  } else {
    writeHashCfgString(fp, "");
  }

  flushBufToDisk(fp);
  fclose(fp);
  umask(origmask);
}

//-------------------------------------------------------------------------------swhitmore
// getNumRefSeqs - Return the number of reference sequences in the hash_table.cfg
// file.
//
uint32_t getNumRefSeqs(const char* filename)
{
  FILE* fp             = fopen(filename, "r");
  char  line[MAX_LINE] = "";
  char* match;

  while (fgets(line, MAX_LINE, fp)) {
    if ((match = strstr(line, "reference_sequences"))) {
      // Before =
      char* val = strtok(match, " = ");
      if (!val) break;
      // After =
      val = strtok(NULL, " = ");
      if (!val) break;
      return strtoul(val, 0, 0);
    }
  }

  fprintf(stderr, "Could not get the number of reference sequences from %s\n", filename);
  return 0;
}

//-------------------------------------------------------------------------------swhitmore
// Populate config with host software version, commmand line, and hash table version.
//
void processComment(const char* line, const size_t lineNum, hashTableConfig_t* config)
{
  char* match;

  if (lineNum == 1) {
    // Extract host software version
    match               = strstr(line, "version") + strlen("version") + 1;
    size_t len          = strlen(match) - 1;
    config->hostVersion = (char*)malloc(len);
    strncpy(config->hostVersion, match, len - 1);
    config->hostVersion[len - 1] = '\0';
  } else if (lineNum == 2) {
    // Extract command line
    match           = strstr(line, "Command line: ") + strlen("Command line: ");
    size_t len      = strlen(match);
    config->cmdLine = (char*)malloc(len);
    strncpy(config->cmdLine, match, len - 1);
    config->cmdLine[len - 1] = '\0';
  } else if (lineNum == 3) {
    // Extract hash table version
    match                         = strstr(line, "Hash table version ") + strlen("Hash table version ");
    config->hdr->hashTableVersion = strtoul(match, 0, 0);
  }
}

//-------------------------------------------------------------------------------seant
// isOldHashTable - Helper function to just check for really old hash tables.
// If the hash table is older than version 3, then there is no version string in the
// config file, and we are unable to detect it further downstream. Instead, we check if
// the version string is missing earlier on and then return an error. This mimics a
// portion of readHashCfgTxt and processComment.
//
int isOldHashTable(const char* filename)
{
  char   line[MAX_LINE] = "";
  size_t lineNum        = 0;
  FILE*  fp;

  if (!(fp = fopen(filename, "r"))) {
    sprintf(errorMsg, "Cannot open config file %s\n", filename);
    return -1;
  }

  // Only need to parse the line where we expect to find the hash table version string
  while (fgets(line, MAX_LINE, fp) && lineNum < 4) {
    ++lineNum;
    if (lineNum == 3) {
      // Extract hash table version
      if (NULL == strstr(line, "Hash table version ")) {
        // Hash tables older than version 3 did not even have a version number
        sprintf(
            errorMsg,
            "Invalid hash table version detected. Check that the hash table in use is at least version %d\n",
            HT_VERSION);
        return -1;
      }
    }
  }

  return 0;
}

//-------------------------------------------------------------------------------swhitmore
// readHashCfgTxt - Reads the given hash_table.cfg text file and populates
// hashTableConfig_t.   Not optimized for performance since host software does
// not read text file unless it is converting it to binary.
//
int readHashCfgTxt(const char* filename, hashTableConfig_t* config)
{
  char     line[MAX_LINE] = "";
  char*    val;
  uint32_t refSeqIdx = 0;
  size_t   i, j;
  size_t   lineNum     = 0;
  int      foundRefSeq = 0;
  // TODO: remove int extTabIdx = 0;

  FILE* fp;
  if (!(fp = fopen(filename, "r"))) {
    sprintf(errorMsg, "Cannot open config file %s\n", filename);
    return 0;
  }

  memset(config, 0, sizeof(hashTableConfig_t));
  config->hdr = (hashTableHeader_t*)malloc(sizeof(hashTableHeader_t));
  memset(config->hdr, 0, sizeof(hashTableHeader_t));

  // Older config files do not define the number of reference sequences so we need
  // do get this number from the file
  config->hdr->numRefSeqs = getNumRefSeqs(filename);
  if (!config->hdr->numRefSeqs) {
    sprintf(
        errorMsg,
        "Hash table %s appears to be version 3 or earlier - version %d required",
        filename,
        HT_VERSION);
    fclose(fp);
    return 0;
  }

  config->refSeq  = (hashTableSeq_t**)malloc(config->hdr->numRefSeqs * sizeof(hashTableSeq_t*));
  config->seqName = (char**)malloc(config->hdr->numRefSeqs * sizeof(char*));

  while (fgets(line, MAX_LINE, fp)) {
    ++lineNum;
    if (!isalpha(line[0])) {
      // extract hostVersion, command line, and hash table version from comments
      processComment(line, lineNum, config);
      continue;
    }
    for (i = 0; !isspace(line[i]); ++i)
      ;
    line[i] = 0;
    for (i++; isspace(line[i]); ++i)
      ;
    if (line[i] != '=') {
      fprintf(stderr, "Invalid format in config file %s on line %zu\n", filename, lineNum);
      fclose(fp);
      return 0;
    }
    for (i++; isspace(line[i]); i++)
      ;
    if (!line[i]) {
      fprintf(stderr, "Invalid format in config file %s on line %zu\n", filename, lineNum);
      exit(1);
    }
    if (line[i] == '\'') {
      ++i;
      for (j = i; line[j] && line[j] != '\''; ++j)
        ;
      if (!line[j]) {
        fprintf(stderr, "Invalid format in config file %s on line %zu\n", filename, lineNum);
        exit(1);
      }
      line[j] = 0;
    } else {
      for (j = i; !isspace(line[j]); j++)
        ;
      line[j] = 0;
    }
    val = &line[i];

    size_t len = strlen(val) + 1;
    if (!strcmp(line, "reference_sequences")) {
      // Make sure the value in the file matches our expected value
      assert(config->hdr->numRefSeqs == strtoul(val, 0, 0));
    } else if (!strncmp(line, "reference_sequence", strlen("reference_sequence"))) {
      foundRefSeq                = 1;
      config->refSeq[refSeqIdx]  = (hashTableSeq_t*)malloc(sizeof(hashTableSeq_t));
      config->seqName[refSeqIdx] = (char*)malloc(strlen(val) + 1);
      strcpy(config->seqName[refSeqIdx], val);
    } else if (!strncmp(line, "reference_start", strlen("reference_start"))) {
      config->refSeq[refSeqIdx]->seqStart = strtoull(val, 0, 0);
    } else if (!strncmp(line, "reference_beg_trim", strlen("reference_beg_trim"))) {
      config->refSeq[refSeqIdx]->begTrim = strtoul(val, 0, 0);
    } else if (!strncmp(line, "reference_end_trim", strlen("reference_end_trim"))) {
      config->refSeq[refSeqIdx]->endTrim = strtoul(val, 0, 0);
    } else if (!strcmp(line, "reference_len_raw")) {
      config->hdr->refLenRaw = strtoull(val, 0, 0);
    } else if (!strcmp(line, "reference_len_not_n")) {
      config->hdr->refLenNotN = strtoull(val, 0, 0);
    } else if (!strcmp(line, "reference_alt_seed")) {
      config->hdr->refAltSeed = strtoul(val, 0, 0);
    } else if (!strcmp(line, "reference_alt_start")) {
      config->hdr->refAltStart = strtoull(val, 0, 0);
    } else if (!strncmp(line, "reference_len", strlen("reference_len"))) {
      if (foundRefSeq) {
        // Reference length of sequence
        config->refSeq[refSeqIdx]->seqLen = strtoul(val, 0, 0);
        foundRefSeq                       = 0;
        ++refSeqIdx;
      } else {
        // Reference length of entire genome
        config->hdr->refSeqLen = strtoull(val, 0, 0);
      }
    } else if (!strcmp(line, "reference_source")) {
      config->refInput = (char*)malloc(len);
      strcpy(config->refInput, val);
    } else if (!strcmp(line, "reference_name")) {
      config->refOutput = (char*)malloc(len);
      strcpy(config->refOutput, val);
    } else if (!strcmp(line, "reference_index")) {
      config->refIdxFname = (char*)malloc(len);
      strcpy(config->refIdxFname, val);
    } else if (!strcmp(line, "hash_table")) {
      config->hashFname = (char*)malloc(len);
      strcpy(config->hashFname, val);
    } else if (!strcmp(line, "extend_table")) {
      config->extTabFname = (char*)malloc(len);
      strcpy(config->extTabFname, val);
    } else if (!strcmp(line, "alt_liftover")) {
      config->altLiftover = (char*)malloc(len);
      strcpy(config->altLiftover, val);
    } else if (!strcmp(line, "pop_alt_contigs")) {
      config->popAltContigsFname = (char*)malloc(len);
      strcpy(config->popAltContigsFname, val);
    } else if (!strcmp(line, "pop_alt_liftover")) {
      config->popAltLiftoverFname = (char*)malloc(len);
      strcpy(config->popAltLiftoverFname, val);
    } else if (!strcmp(line, "pop_snps_source")) {
      config->popSnpsInput = (char*)malloc(len);
      strcpy(config->popSnpsInput, val);
    } else if (!strcmp(line, "pop_snps")) {
      config->popSnpsOutput = (char*)malloc(len);
      strcpy(config->popSnpsOutput, val);
    } else if (!strcmp(line, "mask_bed")) {
      config->maskBed = (char*)malloc(len);
      strcpy(config->maskBed, val);
    } else if (!strcmp(line, "hash_table_bytes")) {
      config->hdr->hashTableBytes = strtoull(val, 0, 0);
    } else if (!strcmp(line, "extend_table_records")) {
      config->hdr->extTabRecs = strtoul(val, 0, 0);
    } else if (!strcmp(line, "pri_seed_bases")) {
      config->hdr->priSeedBases = strtoul(val, 0, 0);
    } else if (!strcmp(line, "max_seed_bases")) {
      config->hdr->maxSeedBases = strtoul(val, 0, 0);
    } else if (!strcmp(line, "max_ext_increment")) {
      config->hdr->maxExtIncrement = strtoul(val, 0, 0);
    } else if (!strcmp(line, "ref_seed_interval")) {
      config->hdr->refSeedInterval = strtod(val, 0);
    } else if (!strcmp(line, "table_addr_bits")) {
      config->hdr->tableAddrBits = strtoul(val, 0, 0);
    } else if (!strcmp(line, "table_size_64ths")) {
      config->hdr->tableSize64ths = strtoul(val, 0, 0);
    } else if (!strcmp(line, "max_seed_freq")) {
      config->hdr->maxSeedFreq = strtoul(val, 0, 0);
    } else if (!strcmp(line, "pri_max_seed_freq")) {
      config->hdr->priMaxSeedFreq = strtoul(val, 0, 0);
    } else if (!strcmp(line, "max_seed_freq_len")) {
      config->hdr->maxSeedFreqLen = strtoul(val, 0, 0);
    } else if (!strcmp(line, "target_seed_freq")) {
      config->hdr->targetSeedFreq = strtod(val, 0);
    } else if (!strcmp(line, "max_mult_base_seeds")) {
      config->hdr->maxMultBaseSeeds = strtoul(val, 0, 0);
    } else if (!strcmp(line, "lift_match_seed_int")) {
      config->hdr->liftMatchSeedInt = strtoul(val, 0, 0);
    } else if (!strcmp(line, "min_freq_to_extend")) {
      config->hdr->minFreqToExtend = strtoul(val, 0, 0);
    } else if (!strcmp(line, "thinning_freq_cap")) {
      config->hdr->thinningFreqCap = strtod(val, 0);
    } else if (!strcmp(line, "max_thinning_factor")) {
      config->hdr->thinningPeriod = strtoul(val, 0, 0);
    } else if (!strcmp(line, "pri_crc_bits")) {
      config->hdr->priCrcBits = strtoul(val, 0, 0);
    } else if (!strcmp(line, "sec_crc_bits")) {
      config->hdr->secCrcBits = strtoul(val, 0, 0);
    } else if (!strcmp(line, "seed_len_cost")) {
      config->hdr->seedLenCost = strtod(val, 0);
    } else if (!strcmp(line, "seed_freq_cost")) {
      config->hdr->seedFreqCost = strtod(val, 0);
    } else if (!strcmp(line, "extension_cost")) {
      config->hdr->extensionCost = strtod(val, 0);
    } else if (!strcmp(line, "ext_step_cost")) {
      config->hdr->extStepCost = strtod(val, 0);
    } else if (!strcmp(line, "ext_rec_cost")) {
      config->hdr->extRecCost = strtod(val, 0);
    } else if (!strcmp(line, "repair_strategy")) {
      config->hdr->repairStrategy = strtoul(val, 0, 0);
    } else if (!strcmp(line, "min_repair_prob")) {
      config->hdr->minRepairProb = strtod(val, 0);
    } else if (!strcmp(line, "anchor_bin_bits")) {
      config->hdr->anchorBinBits = strtoul(val, 0, 0);
    } else if (!strcmp(line, "hi_freq_rand_hit")) {
      config->hdr->hiFreqRandHit = strtoul(val, 0, 0);
    } else if (!strcmp(line, "ext_rand_hit_freq")) {
      config->hdr->extRandHitFreq = strtoul(val, 0, 0);
    } else if (!strcmp(line, "pri_crc_poly")) {
      *((uint64_t*)config->hdr->priCrcPoly) = strtoull(val, 0, 16);
    } else if (!strcmp(line, "sec_crc_poly")) {
      *((uint64_t*)config->hdr->secCrcPoly) = strtoull(val, 0, 16);
    } else if (!strcmp(line, "digest_type")) {
      config->hdr->digestType = strtoul(val, 0, 0);
    } else if (!strcmp(line, "ref_digest")) {
      config->hdr->refDigest = strtoul(val, 0, 0);
    } else if (!strcmp(line, "ref_index_digest")) {
      config->hdr->refIndexDigest = strtoul(val, 0, 0);
    } else if (!strcmp(line, "hash_digest")) {
      config->hdr->hashDigest = strtoul(val, 0, 0);
    } else if (!strcmp(line, "liftover_digest")) {
      config->hdr->liftoverDigest = strtoul(val, 0, 0);
    } else if (!strcmp(line, "extend_table_digest")) {
      config->hdr->extTabDigest = strtoul(val, 0, 0);
    } else if (!strcmp(line, "pop_snps_digest")) {
      config->hdr->popSnpsDigest = strtoul(val, 0, 0);
    } else if (!strcmp(line, "digest")) {
      config->hdr->digest = strtoul(val, 0, 0);
    }
  }
  config->usedReadBuf = 0;
  fclose(fp);
  return 1;
}

//-------------------------------------------------------------------------------swhitmore
// regenHashTableError - Add information to errorMsg on how to regenerate
// the hash table.
//
void regenHashTableError(const char* binFile, const char* liftoverFile)
{
  // binfile is an earlier version of hash_table.cfg that cannot be read, but
  // we can try to read the .txt file
  char* ext = strstr(binFile, ".bin");
  if (!ext) {
    return;
  }

  *ext = '\0';
  hashTableConfig_t config;
  if (!readHashCfgTxt(binFile, &config)) {
    return;
  }

  strcat(errorMsg, "Please regenerate hash table with the following command:\n");

  // Replace --output-dir /staging/foo with --output-dir <NEW-OUTPUT-DIR>
  char* outputStart = strstr(config.cmdLine, "--output-dir");

  if (outputStart) {
    // This is a dragen command
    char* pathStart = strstr(outputStart, " ");
    if (pathStart) {
      while (*pathStart == ' ') {
        // Strip off spaces after --output-dir
        *pathStart = '\0';
        pathStart += 1;
      }
      pathStart += 1;
      char* restOfArgs = strstr(pathStart, " ");

      sprintf(errorMsg + strlen(errorMsg), "  %s <NEW-OUTPUT-DIR> ", config.cmdLine);
      if (liftoverFile) {
        sprintf(errorMsg + strlen(errorMsg), "--ht-alt-liftover %s ", liftoverFile);
      }
      if (restOfArgs) {
        strncat(errorMsg + strlen(errorMsg), restOfArgs, 2048 - strlen(errorMsg) - 1);
      } else {
        sprintf(errorMsg + strlen(errorMsg), "\n");
      }
      return;
    }
  }
  sprintf(errorMsg + strlen(errorMsg), "  Old: %s\n", config.cmdLine);
  sprintf(errorMsg + strlen(errorMsg), "  New: dragen --build-hash-table ");
  sprintf(errorMsg + strlen(errorMsg), "--output-directory <NEW-OUTPUT-DIR> ");
  sprintf(errorMsg + strlen(errorMsg), "--ht-reference <PATH-TO-FASTA> ");
  if (liftoverFile) {
    sprintf(errorMsg + strlen(errorMsg), "--ht-alt-liftover %s ", liftoverFile);
  }
  sprintf(errorMsg + strlen(errorMsg), "<args>\n");
}

//-------------------------------------------------------------------------------swhitmore
// getBuildHashTableCmdline - Return message with instructions for rebuilding the hash
// table.
//
const char* getBuildHashTableCmdline(const char* binFile, const char* liftoverFile)
{
  regenHashTableError(binFile, liftoverFile);
  return getHashCfgError();
}

//-------------------------------------------------------------------------------swhitmore
// readHashCfgBin - Reads the given hash_table.cfg.bin binary file and populates
// hashTableConfig_t.  Returns 1 on success and 0 on error.
//
int readHashCfgBin(const char* filename, hashTableConfig_t* config)
{
  if (access(filename, R_OK) == -1) {
    sprintf(errorMsg, "Could not open %s -- %s\n", filename, strerror(errno));
    return 0;
  }
  FILE* fp = fopen(filename, "rb");
  fseek(fp, 0L, SEEK_END);
  long fileSize = ftell(fp);
  fseek(fp, 0L, SEEK_SET);

  if (fileSize <= HASH_HEADER_SIZE) {
    sprintf(errorMsg, "%s is invalid - filesize = %lu\n", filename, fileSize);
    fclose(fp);
    return 0;
  }

  config->readBuf = (uint8_t*)malloc(fileSize);
  if (fread(config->readBuf, fileSize, 1, fp) != fileSize)
  {
    sprintf(errorMsg, "failed to read %lu bytes from %s: %s\n", fileSize, filename, (feof(fp) ? "EOF" : strerror(errno)));
    fclose(fp);
    return 0;
  }

  config->hdr = (hashTableHeader_t*)config->readBuf;
  // size_t bufOffset = sizeof(hashTableHeader_t);
  long   bufOffset = sizeof(hashTableHeader_t);
  size_t i;

  int hashTableVersion = config->hdr->hashTableVersion;
  if (hashTableVersion == HT_VERSION) {
    config->refSeq = (hashTableSeq_t**)malloc(config->hdr->numRefSeqs * sizeof(hashTableSeq_t*));
    for (i = 0; i < config->hdr->numRefSeqs; ++i) {
      config->refSeq[i] = (hashTableSeq_t*)(config->readBuf + bufOffset);
      bufOffset += sizeof(hashTableSeq_t);
    }
  } else {
    char dir[PATH_MAX];
    strcpy(dir, filename);
    sprintf(
        errorMsg,
        "Hash table %s is version %d - version %d is required\n",
        dirname(dir),
        hashTableVersion,
        HT_VERSION);
    if (hashTableVersion > 2 && hashTableVersion < HT_VERSION) {
      // Note: hash table versions 1 and 2 do not have a command line
      regenHashTableError(filename, NULL /* no liftover */);
    } else if ((hashTableVersion > HT_VERSION + 1) || (hashTableVersion < 1)) {
      // If hash table version is more than 1 version ahead, consider the file invalid
      sprintf(errorMsg, "%s is invalid - hash table version is larger than %d", filename, HT_VERSION + 1);
    }
    return 0;
  }

  config->seqName = (char**)malloc(config->hdr->numRefSeqs * sizeof(char*));
  for (i = 0; i < config->hdr->numRefSeqs; ++i) {
    config->seqName[i] = (char*)(config->readBuf + bufOffset);
    bufOffset += strlen(config->seqName[i]) + 1;
  }

  config->hostVersion = (char*)(config->readBuf + bufOffset);
  bufOffset += strlen(config->hostVersion) + 1;
  config->cmdLine = (char*)(config->readBuf + bufOffset);

  bufOffset += strlen(config->cmdLine) + 1;
  config->refInput = (char*)(config->readBuf + bufOffset);
  bufOffset += strlen(config->refInput) + 1;
  config->refOutput = (char*)(config->readBuf + bufOffset);
  bufOffset += strlen(config->refOutput) + 1;
  config->refIdxFname = (char*)(config->readBuf + bufOffset);
  bufOffset += strlen(config->refIdxFname) + 1;
  config->hashFname = (char*)(config->readBuf + bufOffset);
  bufOffset += strlen(config->hashFname) + 1;
  config->extTabFname = (char*)(config->readBuf + bufOffset);
  bufOffset += strlen(config->extTabFname) + 1;
  config->altLiftover = (char*)(config->readBuf + bufOffset);
  bufOffset += strlen(config->altLiftover) + 1;

  // Optional values
  config->popAltContigsFname  = NULL;
  config->popAltLiftoverFname = NULL;
  config->popSnpsInput        = NULL;
  config->popSnpsOutput       = NULL;
  config->maskBed             = NULL;

  // If there is still file content, read all the optional values
  // If the value is empty, leave it as NULL
  if (bufOffset < fileSize) {
    char* p;

    p = (char*)(config->readBuf + bufOffset);
    if (*p) config->maskBed = p;
  }

  config->usedReadBuf = 1;

  fclose(fp);
  return 1;
}

//-------------------------------------------------------------------------------swhitmore
// writeHashCfgTxt - write hash table configuration file as text.
//
void writeHashCfgTxt(const char* filename, const hashTableConfig_t* config)
{
  FILE*              fp    = fopen(filename, "w");
  hashTableHeader_t* hdr   = config->hdr;
  char*              cmdcp = (char*)malloc(strlen(config->cmdLine) + 1);
  strcpy(cmdcp, config->cmdLine);

#ifdef LOCAL_BUILD
  fprintf(fp, "# Automatically generated by build_hash_table (local version).\n");

#else
  char* cmdline = strtok(cmdcp, " ");
  fprintf(fp, "# Automatically generated by %s version %s.\n", basename(cmdline), config->hostVersion);
#endif
  fprintf(fp, "#    Command line: %s", config->cmdLine);
  fprintf(fp, "\n");
  fprintf(fp, "#    Hash table version %d\n", config->hdr->hashTableVersion);
  fprintf(fp, "#\n");
  fprintf(fp, "# Do not modify.\n");
  fprintf(fp, "\n");

  fprintf(fp, "reference_source     = '%s'\n", config->refInput);
  if (config->altLiftover) {
    fprintf(fp, "alt_liftover         = '%s'\n", config->altLiftover);
  } else {
    fprintf(fp, "alt_liftover         = ''\n");
  }
  if (config->maskBed) {
    fprintf(fp, "mask_bed             = '%s'\n", config->maskBed);
  }
  if (config->popAltContigsFname) {
    fprintf(fp, "pop_alt_contigs      = '%s'\n", config->popAltContigsFname);
  }
  if (config->popAltLiftoverFname) {
    fprintf(fp, "pop_alt_liftover     = '%s'\n", config->popAltLiftoverFname);
  }
  if (config->popSnpsInput) {
    fprintf(fp, "pop_snps_source      = '%s'\n", config->popSnpsInput);
  }
  fprintf(fp, "reference_name       = '%s'\n", config->refOutput);
  fprintf(fp, "reference_index      = '%s'\n", config->refIdxFname);
  fprintf(fp, "reference_sequences  = %lld\n", config->hdr->numRefSeqs);
  fprintf(fp, "reference_len        = %lld\n", (long long int)hdr->refSeqLen);
  fprintf(fp, "reference_len_raw    = %lld\n", (long long int)hdr->refLenRaw);
  fprintf(fp, "reference_len_not_n  = %lld\n", (long long int)hdr->refLenNotN);
  fprintf(fp, "reference_alt_seed   = %lu\n", (long unsigned)hdr->refAltSeed);
  fprintf(fp, "reference_alt_start  = %lld\n", (long long int)hdr->refAltStart);
  fprintf(fp, "hash_table           = '%s'\n", config->hashFname);
  fprintf(fp, "extend_table         = '%s'\n", config->extTabFname);
  fprintf(fp, "hash_table_bytes     = %lld\n", (long long int)hdr->hashTableBytes);
  fprintf(fp, "extend_table_records = %lu\n", (long unsigned)hdr->extTabRecs);
  if (config->popSnpsInput && config->popSnpsOutput) {
    fprintf(fp, "pop_snps             = '%s'\n", config->popSnpsOutput);
  }

  fprintf(fp, "digest_type          = %lu\n", (long unsigned)hdr->digestType);
  fprintf(fp, "digest               = 0x%08lX\n", (long unsigned)hdr->digest);
  fprintf(fp, "ref_digest           = 0x%08lX\n", (long unsigned)hdr->refDigest);
  fprintf(fp, "ref_index_digest     = 0x%08lX\n", (long unsigned)hdr->refIndexDigest);
  fprintf(fp, "hash_digest          = 0x%08lX\n", (long unsigned)hdr->hashDigest);
  fprintf(fp, "liftover_digest      = 0x%08lX\n", (long unsigned)hdr->liftoverDigest);
  fprintf(fp, "extend_table_digest  = 0x%08lX\n", (long unsigned)hdr->extTabDigest);
  if (hdr->popSnpsDigest != 0) {
    fprintf(fp, "pop_snps_digest      = 0x%08lX\n", (long unsigned)hdr->popSnpsDigest);
  }
  if (config->maskBed) {
    fprintf(fp, "mask_bed_digest      = 0x%08lX\n", (long unsigned)config->maskBedDigest);
  }

  fprintf(fp, "pri_seed_bases       = %u\n", hdr->priSeedBases);
  fprintf(fp, "max_seed_bases       = %u\n", hdr->maxSeedBases);
  fprintf(fp, "max_ext_increment    = %u\n", hdr->maxExtIncrement);
  fprintf(fp, "ref_seed_interval    = %g\n", hdr->refSeedInterval);
  fprintf(fp, "table_addr_bits      = %u\n", hdr->tableAddrBits);
  fprintf(fp, "table_size_64ths     = %u\n", hdr->tableSize64ths);
  fprintf(fp, "max_seed_freq        = %u\n", hdr->maxSeedFreq);
  fprintf(fp, "pri_max_seed_freq    = %u\n", hdr->priMaxSeedFreq);
  fprintf(fp, "max_seed_freq_len    = %u\n", hdr->maxSeedFreqLen);
  fprintf(fp, "target_seed_freq     = %g\n", hdr->targetSeedFreq);
  if (hdr->maxMultBaseSeeds > 0) {
    fprintf(fp, "max_mult_base_seeds  = %u\n", hdr->maxMultBaseSeeds);
  }
  if (hdr->liftMatchSeedInt > 0) {
    fprintf(fp, "lift_match_seed_int  = %u\n", hdr->liftMatchSeedInt);
  }
  fprintf(fp, "min_freq_to_extend   = %u\n", hdr->minFreqToExtend);
  fprintf(fp, "thinning_freq_cap    = %g\n", hdr->thinningFreqCap);
  fprintf(fp, "max_thinning_factor  = %u\n", hdr->thinningPeriod);
  fprintf(fp, "pri_crc_bits         = %u\n", hdr->priCrcBits);
  fprintf(fp, "sec_crc_bits         = %u\n", hdr->secCrcBits);
  fprintf(fp, "seed_len_cost        = %g\n", hdr->seedLenCost);
  fprintf(fp, "seed_freq_cost       = %g\n", hdr->seedFreqCost);
  fprintf(fp, "extension_cost       = %g\n", hdr->extensionCost);
  fprintf(fp, "ext_step_cost        = %g\n", hdr->extStepCost);
  fprintf(fp, "ext_rec_cost         = %g\n", hdr->extRecCost);
  fprintf(fp, "repair_strategy      = %u\n", hdr->repairStrategy);
  fprintf(fp, "min_repair_prob      = %g\n", hdr->minRepairProb);
  fprintf(fp, "anchor_bin_bits      = %u\n", hdr->anchorBinBits);
  fprintf(fp, "hi_freq_rand_hit     = %u\n", hdr->hiFreqRandHit);
  fprintf(fp, "ext_rand_hit_freq    = %u\n", hdr->extRandHitFreq);
  fprintf(fp, "pri_crc_poly         = %llX\n", (long long unsigned)((uint64_t*)hdr->priCrcPoly)[0]);
  fprintf(fp, "sec_crc_poly         = %llX\n", (long long unsigned)((uint64_t*)hdr->secCrcPoly)[0]);
  fflush(fp);

  unsigned int i;
  for (i = 0; i < config->hdr->numRefSeqs; ++i) {
    fprintf(fp, "reference_sequence%-2u    = '%s'\n", i, config->seqName[i]);
    fprintf(fp, "reference_start%-2u       = %llu\n", i, config->refSeq[i]->seqStart);
    fprintf(fp, "reference_beg_trim%-2u    = %lu\n", i, config->refSeq[i]->begTrim);
    fprintf(fp, "reference_end_trim%-2u    = %lu\n", i, config->refSeq[i]->endTrim);
    fprintf(fp, "reference_len%-2u         = %u\n", i, config->refSeq[i]->seqLen);
  }

  fflush(fp);

  free(cmdcp);
}
