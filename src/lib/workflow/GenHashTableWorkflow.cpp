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
//
//#include <fstream>
//#include <limits>
//#include <cerrno>
//#include <cstring>
//
//#include "boost/iostreams/filter/gzip.hpp"
//
//#include <boost/iostreams/device/back_inserter.hpp>
//#include <boost/iostreams/filtering_stream.hpp>
//#include <boost/filesystem.hpp>
//
#include "common/Debug.hpp"
#include "options/DragenOsOptions.hpp"
//
#include "common/hash_generation/gen_hash_table.h"
#include "common/hash_generation/hash_table_compress.h"
#include "workflow/GenHashTableWorkflow.hpp"

#include "common/public/linux_utils.hpp"

namespace dragenos {
namespace workflow {
std::string GetFullPath(const std::string& p)
{
  boost::filesystem::path path(p);
  boost::filesystem::path cpath(p);
  return cpath.string();
}

/// Return whether we should be strict about disabling pipelines that are not valid for the detected
/// reference.
bool ReferenceAutoDetectValidate(const options::DragenOsOptions& opts)
{
  return false;
}

/// Get the directory that has the reference autodetection config files.
std::string GetReferenceAutoDetectDir(const options::DragenOsOptions& opts)
{
  // Always get the path relative to the build directory
  std::string       dir;
  const std::string autodetectDir = "reference_autodetect/";
  getBuildPath(dir, autodetectDir);
  return dir;
}

//-------------------------------------------------------------------------------swhitmore
// getRelativeLiftoverPath - Returns the liftover file location relative to the path to
// the dragen binary.
//
void getRelativeLiftoverPath(const std::string& oldPath, std::string& newPath)
{
  char buf[PATH_MAX];

  ssize_t len = ::readlink("/proc/self/exe", buf, sizeof(buf));
  if (len != -1) {
    buf[len]  = '\0';
    char* ptr = std::strrchr(buf, '/');
    if (ptr != NULL) {
      ptr[0] = '\0';
    }
  } else {
    buf[0] = '\0';
  }

  boost::filesystem::path dragen_exe_dir(buf);
  boost::filesystem::path p(oldPath);
  newPath = dragen_exe_dir.parent_path().string() + "/liftover/" + p.filename().string();
}
//-------------------------------------------------------------------------------swhitmore
// SetBuildHashTableOptions - Populate the hashTableConfig_t struct from config
// file options.
// Every allocation here should have a corresponding free in FreeBuildHashTableOptions()
//
void SetBuildHashTableOptions(
    const options::DragenOsOptions& opts, hashTableConfig_t* config, HashTableType hashTableType)
{
  std::string ref  = GetFullPath(opts.htReference_);
  config->refInput = strdup(ref.c_str());

  config->altLiftover = NULL;

  config->decoyFname     = NULL;
  std::string decoysFile = GetFullPath(opts.htDecoys_);
  std::cerr << "Supressing decoys" << std::endl;

  if (!decoysFile.empty()) {
    assert(boost::filesystem::exists(decoysFile));
    std::cerr << "Using decoys file " << decoysFile << std::endl;
    config->decoyFname = strdup(decoysFile.c_str());
  }

  config->popAltContigsFname = NULL;

  config->popAltLiftoverFname = NULL;

  config->popSnpsInput = NULL;

  config->maskBed     = NULL;
  std::string maskBed = GetFullPath(opts.htMaskBed_);
  if (!maskBed.empty()) {
    config->maskBed = strdup(maskBed.c_str());
  }

  config->hdr->refSeedInterval = opts.htRefSeedInterval_;
  config->hdr->priSeedBases    = opts.htSeedLen_;
  if (opts.htMaxExtSeedLen_) {
    config->hdr->maxSeedBases = opts.htMaxExtSeedLen_;
  } else {
    // Value is 0 - set to 999 to trigger default behavior
    config->hdr->maxSeedBases = 999;
  }
  config->hdr->maxSeedFreq     = opts.htMaxSeedFreq_;
  config->hdr->priMaxSeedFreq  = opts.htPriMaxSeedFreq_;
  config->hdr->maxSeedFreqLen  = opts.htMaxSeedFreqLen_;
  config->hdr->targetSeedFreq  = opts.htTargetSeedFreq_;
  config->hdr->thinningFreqCap = opts.htSoftSeedFreqCap_;
  config->hdr->thinningPeriod  = opts.htMaxDecFactor_;
  config->hdr->seedLenCost     = opts.htCostCoeffSeedLen_;
  config->hdr->seedFreqCost    = opts.htCostCoeffSeedFreq_;
  config->hdr->extensionCost   = opts.htCostPenalty_;
  config->hdr->extStepCost     = opts.htCostPenaltyIncr_;
  config->hdr->repairStrategy  = opts.htRepairStrategy_;
  config->hdr->anchorBinBits   = opts.htAnchorBinBits_;
  config->hdr->minRepairProb   = opts.htMinRepairProb_;
  config->hdr->hiFreqRandHit   = opts.htRandHitHifreq_;
  config->hdr->extRandHitFreq  = opts.htRandHitExtend_;
  config->hdr->extRecCost      = opts.htExtRecCost_;
  if (opts.exists("ht-max-multi-base-seeds")) {
    config->hdr->maxMultBaseSeeds = opts.htMaxMultiBaseSeeds_;
  } else {
    if (config->popSnpsInput) {
      config->hdr->maxMultBaseSeeds = 64;
    } else {
      config->hdr->maxMultBaseSeeds = 0;
    }
  }

  config->maxThreads =
      opts.exists("ht-num-threads") ? opts.htNumThreads_ : boost::thread::hardware_concurrency();
  config->sizeStr    = strdup(opts.htSize_.c_str());
  config->sjSizeStr  = opts.exists("ht-sj-size") ? strdup(opts.htSjSize_.c_str()) : NULL;
  config->memSizeStr = strdup(opts.htMemLimit_.c_str());
  if (opts.exists("ht-max-table-chunks")) {
    config->maxGB = opts.htMaxTableChunks_;
  }
  config->hdr->maxExtIncrement = opts.htMaxExtIncr_;
  config->extTableAlloc        = opts.htExtTableAlloc_;
  config->priPolyIndex         = opts.htCrcPrimary_;
  config->secPolyIndex         = opts.htCrcExtended_;
  config->methylatedConv       = 0;
  config->overrideCheck        = opts.htOverrideSizeCheck_;
  config->testOnly             = opts.htTestOnly_;
  config->showIntParams        = opts.htDumpIntParams_;
  config->writeCompFile        = 1;
  config->writeHashFile        = opts.htWriteHashBin_;
  config->altContigValidate    = false;
  config->autoDetectValidate   = ReferenceAutoDetectValidate(opts);
  config->autoDetectDir        = strdup(GetReferenceAutoDetectDir(opts).c_str());

  // Override parameters which are set based on the hash table type
  // Note: a default value should have been specified above
  if (hashTableType == HT_TYPE_METHYL_G_TO_A or hashTableType == HT_TYPE_METHYL_C_TO_T) {
    config->methylatedConv = hashTableType;
  }
  if (hashTableType == HT_TYPE_ANCHORED) {
    config->hdr->priSeedBases  = 11;
    config->hdr->maxSeedBases  = 11;
    config->hdr->maxSeedFreq   = 1;
    config->hdr->anchorBinBits = 16;
  }
}

void FreeBuildHashTableOptions(hashTableConfig_t* config)
{
  freeHashParams(config);

  if (config->refInput) {
    free(config->refInput);
    config->refInput = NULL;
  }

  if (config->altLiftover) {
    free(config->altLiftover);
    config->altLiftover = NULL;
  }

  if (config->decoyFname) {
    free(config->decoyFname);
    config->decoyFname = NULL;
  }

  if (config->popAltContigsFname) {
    free(config->popAltContigsFname);
    config->popAltContigsFname = NULL;
  }

  if (config->popAltLiftoverFname) {
    free(config->popAltLiftoverFname);
    config->popAltLiftoverFname = NULL;
  }

  if (config->popSnpsInput) {
    free(config->popSnpsInput);
    config->popSnpsInput = NULL;
  }

  if (config->maskBed) {
    free(config->maskBed);
    config->maskBed = NULL;
  }

  if (config->sizeStr) {
    free((char*)config->sizeStr);
    config->sizeStr = NULL;
  }

  if (config->sjSizeStr) {
    free((char*)config->sjSizeStr);
    config->sjSizeStr = NULL;
  }

  if (config->memSizeStr) {
    free((char*)config->memSizeStr);
    config->memSizeStr = NULL;
  }

  if (config->autoDetectDir) {
    free((char*)config->autoDetectDir);
    config->autoDetectDir = NULL;
  }
}

//-------------------------------------------------------------------------------swhitmore
// uncompressHashComp - Uncompress hash_table.cmp (also uses reference.bin) to
// hash_table.bin.  For debugging purposes only.
//
void uncompressHashCmp(const options::DragenOsOptions& opts)
{
  std::string refdir = opts.htReference_;
  std::cout << "==================================================================\n";
  std::cout << "Uncompressing " << refdir << " hash_table.cmp to hash_table.bin & extend_table.bin\n";
  std::cout << "==================================================================\n";
  std::string refPath      = refdir + "/reference.bin";
  std::string hashCmpPath  = refdir + "/hash_table.cmp";
  std::string hashBinPath  = refdir + "/hash_table.bin";
  std::string extTablePath = refdir + "/extend_table.bin";
  int         numThreads   = boost::thread::hardware_concurrency();
  if (!decompAndWriteHashTable(
          refPath.c_str(), hashCmpPath.c_str(), hashBinPath.c_str(), extTablePath.c_str(), numThreads)) {
    BOOST_THROW_EXCEPTION(std::logic_error(std::string("Could not decompress ") << hashCmpPath));
  }
}

void buildHashTable(const options::DragenOsOptions& opts)
{
  if (opts.htUncompress_) {
    uncompressHashCmp(opts);
    return;
  }

  if (!opts.buildHashTable_) {
    return;
  }

  DRAGEN_OS_THREAD_CERR << "Version: " << common::Version::string() << std::endl;
  DRAGEN_OS_THREAD_CERR << "argc: " << opts.argc() << " argv: " << opts.getCommandLine() << std::endl;

  std::cout << "==================================================================\n";
  std::cout << "Building hash table from " << opts.htReference_ << "\n";
  std::cout << "==================================================================\n";

  // Determine the type of hash table to generate
  std::vector<HashTableType> hashTableTypes;
  //  if (opts->BuildMethylatedHashtables()) {
  //    hashTableTypes.push_back(HT_TYPE_METHYL_G_TO_A);
  //    hashTableTypes.push_back(HT_TYPE_METHYL_C_TO_T);
  //    hashTableTypes.push_back(HT_TYPE_NORMAL);
  //  } else {
  hashTableTypes.push_back(HT_TYPE_NORMAL);

  // If the user specifies to build RNA tables
  //    if (opts->BuildRNAHashTable()) {
  //      hashTableTypes.push_back(HT_TYPE_ANCHORED);
  //    }
  //  }

  for (auto it = hashTableTypes.begin(); it != hashTableTypes.end(); ++it) {
    hashTableConfig_t bhtConfig;
    memset(&bhtConfig, 0, sizeof(bhtConfig));
    // Set parameters shared with build_hash_table that the user cannot change
    boost::filesystem::path p(opts.outputDirectory_);
    std::string             outdir = boost::filesystem::canonical(p).string();
    setDefaultHashParams(&bhtConfig, outdir.c_str(), (*it));
    SetBuildHashTableOptions(opts, &bhtConfig, *it);

    char* errorMsg = generateHashTable(&bhtConfig, opts.argc(), const_cast<char**>(opts.argv()));
    if (errorMsg) {
      std::cerr << "Hash table generation failed: " << errorMsg;
      BOOST_THROW_EXCEPTION(
          std::logic_error(std::string("ERROR: Hash table generation failed - ") << errorMsg));
    }

    FreeBuildHashTableOptions(&bhtConfig);
  }
}
}  // namespace workflow
}  // namespace dragenos
