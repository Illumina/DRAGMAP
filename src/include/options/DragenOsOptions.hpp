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

#ifndef OPTIONS_DRAGEN_OS_OPTIONS_HPP
#define OPTIONS_DRAGEN_OS_OPTIONS_HPP

#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include <string>
#include <thread>

#include "common/Program.hpp"
#include "common/hash_generation/gen_hash_table.h"

namespace dragenos {
namespace options {

class DragenOsOptions : public common::Options {
public:
  DragenOsOptions();

private:
  std::string usagePrefix() const { return "dragenos -r <reference> -b <base calls> [optional arguments]"; }
  void        postProcess(boost::program_options::variables_map& vm);
  void        SetBuildHashTableOptions(hashTableConfig_t* config, HashTableType hashTableType);

public:
  std::string             description_;
  boost::filesystem::path refDir_;
  bool                    mmapReference_ = false;
  bool                    loadReference_ = false;
  std::string             inputFile1_;
  std::string             inputFile2_;
  std::string             outputDirectory_  = "";
  std::string             outputFilePrefix_ = "";

  std::string rgid_ = "1";
  std::string rgsm_ = "none";

  bool interleaved_ = false;
  //bool mapperCigar_;
  bool mapOnly_;
  int  swAll_ = 0;  // Aligner.sw-all

  std::string methodSmithWatermanDeprecated_;
  std::string methodSmithWaterman_ =
      "mengyao";  // "mengyao" : vectorized SW library /   "dragen" for legacy code

  bool        samplingEnabled_       = true;  // sampling-enabled
  double      alignerPeMeanInsert_   = 0.0;   // Aligner.pe-stat-mean-insert
  double      alignerPeStddevInsert_ = 0.0;   // Aligner.pe-stat-stddev-insert
  std::string alignerPeQuartilesInsert_;
  int         alignerPeq25_ = 0, alignerPeq50_ = 0, alignerPeq75_ = 0;  // Aligner.pe-stat-quartiles-insert
  double      alignerPeMeanReadLen_ = 0;                                // Aligner.pe-stat-mean-read-len

  uint64_t peStatsIntervalSize_     = 25000;   // pe-stats-interval-size
  uint64_t peStatsSampleSize_       = 100000;  // pe-stats-sample-size
  uint8_t  peStatsIntervalMemory_   = 10;      // pe-stats-interval-memory
  uint8_t  peStatsIntervalDelay_    = 5;       // pe-stats-interval-delay
  bool     peStatsContinuousUpdate_ = false;   // pe-stats-continuous-update
  bool     peStatsUpdateLogOnly_    = false;   // pe-stats-update-log-only

  double mapperFilterLenRatio_ = 4.0;  // Mapper.filter-len-ratio

  uint32_t alignerPeOrientation_    = 0;    // Aligner.pe-orientation
  double   alignerResqueSigmas_     = 0;    //2.5;     // Aligner.rescue-sigmas
  double   alignerResqueCeilFactor_ = 3.0;  // Aligner.rescue-ceil-factor
  uint32_t alignerResqueMinIns_     = 0;    // Aligner.rescue-min-ins
  uint32_t alignerResqueMaxIns_     = 0;    // Aligner.rescue-max-ins

  uint32_t alignerMapqMinLen_  = 50;  // Aligner.mapq-min-len
  uint32_t alignerMapqMax_     = 60;  // Aligner.mapq-max
  uint32_t alignerUnpairedPen_ = 80;  // Aligner.unpaired-pen
  int      alignerXsPairPen_   = 25;  // Aligner.xs-pair-penalty
  int      alignerSampleMapq0_ = 1;   // Aligner.sample-map0

  int  alignerSecAligns_     = 0;      // Aligner.sec-aligns
  int  alignerSecScoreDelta_ = 0;      // Aligner.sec-score-delta
  int  alignerSecPhredDelta_ = 0;      // Aligner.sec-phred-delta
  bool alignerSecAlignsHard_ = false;  // Aligner.sec-aligns-hard
  int  mapperNumThreads_ =
      0;  // Maximum worker threads for map /align. If not defined, then use maximum on system.
  const int matchScore_       = 1;
  const int mismatchScore_    = -4;
  const int gapExtendPenalty_ = 1;
  const int gapInitPenalty_   = 7;
  const int unclipScore_      = 5;
  int       alnMinScore_      = 22 * matchScore_;
  int       suppMinScoreAdj_  = 8;
  int       mapq_min_len_     = 50;

  char inputQnameSuffixDelim_ = ' ';
  int  fastqOffset_           = 33;

  bool preserveMapAlignOrder_ = false;

  bool        verbose_        = false;
  bool        buildHashTable_ = false;
  bool        htUncompress_   = false;
  std::string htReference_;

  // Number of positions per reference seed
  double htRefSeedInterval_ = 1.0;
  // Initial seed length to store in hash table
  int htSeedLen_ = 21;
  // Maximum extended seed length.  When set to 0, host software will automatically
  // calculate ht-max-ext-seed-len as ht-seed-len + 128.
  int htMaxExtSeedLen_ = 0;
  // Maximum allowed frequency for a seed match after extension attempts (1-256)
  int htMaxSeedFreq_ = 16;

  // Target seed frequency for seed extension
  double htTargetSeedFreq_ = 4.0;
  // Soft seed frequency cap for thinning
  double htSoftSeedFreqCap_ = 12.0;
  // Maximum decimation factor for seed thinning (1-16)
  int htMaxDecFactor_ = 1;
  // Cost coefficient of extended seed length
  double htCostCoeffSeedLen_ = 1.0;
  // Cost coefficient of extended seed frequency
  double htCostCoeffSeedFreq_ = 0.5;
  // Cost penalty to extend a seed by any number of bases
  double htCostPenalty_ = 0.0;
  // Cost penalty to incrementally extend a seed another step
  double htCostPenaltyIncr_ = 0.7;
  // Include a random hit with each HIFREQ record
  int htRandHitHifreq_ = 1;
  // Include a random hit with each EXTEND record of this frequency
  int htRandHitExtend_ = 8;
  // Whether to automatically generate C->T and G->A converted hashtables
  bool htMethylated_ = false;
  // Maximum worker threads. If not defined, then use maximum on system.
  int htNumThreads_ = 8;
  // Size of hash table, units B|KB|MB|GB.  When set to 0, automatically calculated by host software
  std::string htSize_ = "0GB";
  // Reserved space for RNA annotated SJs, units B|KB|MB|GB.  When undefined, automatically calculated by host
  // software
  std::string htSjSize_;
  // Memory limit (hash table + reference), units B|KB|MB|GB.  When set to 0, automatically
  // calculated by the host software - 32GB for the whole human genome
  std::string htMemLimit_ = "0GB";

  // Path to decoys file (FASTA format) - defaults to FASTA in /opt/edico/liftover
  std::string htDecoys_;

  std::string htMaskBed_;  //--ht-mask-bed_

  // Maximum frequency for a primary seed match (0 => use maxSeedFreq)
  int htPriMaxSeedFreq_ = 0;
  // Ramp from priMaxSeedFreq reaches maxSeedFreq at this seed length
  int htMaxSeedFreqLen_ = 98;
  // Maximum bases to extend a seed by in one step
  int htMaxExtIncr_ = 12;
  // 8-byte records to reserve in extend_table.bin (0=automatic)
  int htExtTableAlloc_ = 0;
  // Index of CRC polynomial for hashing primary seeds
  int htCrcPrimary_ = 0;
  // Index of CRC polynomial for hashing extended seeds
  int htCrcExtended_ = 0;
  // Override hash table size check [0|1]
  int htOverrideSizeCheck_ = 0;
  // Testing - show user parameters, but don't do anything [0|1]
  bool htTestOnly_ = false;
  // Show internal parameters [0|1]
  bool htDumpIntParams_ = false;
  // Maximum ~1GB thread table chunks in memory at once.  If set to 0, software will
  // automatically calculate value.
  int htMaxTableChunks_ = 0;
  // Seed extension repair strategy: 0=none, 1=best, 2=rand
  int htRepairStrategy_ = 0;
  // Minimum probability of success to attempt extended seed repair
  double htMinRepairProb_ = 0.2;
  // Bits defining reference bins for anchored seed search. Set to 0 to disable
  int htAnchorBinBits_ = 0;
  // If uncompressed hash_table.bin should be written during hash table build
  bool htWriteHashBin_ = false;
  // Cost penalty for each EXTEND or INTERVAL record
  double htExtRecCost_ = 4.0;
  // Maximum seeds populated at multi-base codes
  int htMaxMultiBaseSeeds_ = 0;
};

}  // namespace options
}  // namespace dragenos

#endif  // #ifndef OPTIONS_DRAGEN_OS_OPTIONS_HPP
