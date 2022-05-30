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

#include <algorithm>
#include <fstream>
#include <ostream>
#include <string>
#include <typeinfo>
#include <vector>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/regex.hpp>
#include <boost/assign.hpp>
#include <boost/foreach.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/thread.hpp>

#include "common/Exceptions.hpp"
#include "common/Version.hpp"
#include "options/DragenOsOptions.hpp"

namespace dragenos {
namespace options {

namespace bpo = boost::program_options;
namespace bfs = boost::filesystem;
using boost::format;
using common::InvalidOptionException;

DragenOsOptions::DragenOsOptions() : refDir_("./"), inputFile1_(""), inputFile2_(""), mapOnly_(false)
{
  // deprecated command line options. Still valid but will error when conflict with official ones or warning
  // when the corresponding official is not being used instead.
  unnamedOptions_.add_options()(
      "Aligner.sw-method",
      bpo::value<std::string>(&methodSmithWatermanDeprecated_),
      "Smith Waterman implementation (dragen / mengyao  default = dragen)");

  namedOptions_.add_options()(
      "ref-dir,r",
      bpo::value<decltype(refDir_)>(&refDir_),
      "directory with reference and hash tables. Must contain the uncompressed hashtable.")(
      "fastq-file1,1", bpo::value<std::string>(&inputFile1_), "FASTQ file to send to card (may be gzipped)")(
      "fastq-file2,2",
      bpo::value<std::string>(&inputFile2_),
      "Second FASTQ file with paired-end reads (may be gzipped)")(
      "bam-input,b", bpo::value<std::string>(&inputFile1_), "Input BAM file")(
      "interleaved",
      bpo::value<bool>(&interleaved_)->default_value(interleaved_)->implicit_value(true),
      "Interleaved paired-end reads in single bam or FASTQ")
      //("mapper_cigar"   , bpo::value<bool>(&mapperCigar_),
      //        "no real alignment, produces alignment information based on seed chains only -- dragen
      //        legacy")
      ("mmap-reference",
       bpo::value<bool>(&mmapReference_)->default_value(mmapReference_),
       "memory-map reference data instead of pre-loading. This allows for quicker runs when only a "
       "handful of reads need to be aligned")(
          "RGID", bpo::value<std::string>(&rgid_)->default_value(rgid_), "Read Group ID")(
          "RGSM", bpo::value<std::string>(&rgsm_)->default_value(rgsm_), "Read Group Sample")(
          "output-directory",
          bpo::value<std::string>(&outputDirectory_)->default_value(outputDirectory_),
          "Output directory")(
          "output-file-prefix",
          bpo::value<std::string>(&outputFilePrefix_)->default_value(outputFilePrefix_),
          "Output filename prefix")(
          "ref-load-hash-bin",
          bpo::value<bool>(&loadReference_)->default_value(loadReference_),
          "Expect to find uncompressed hash table in the reference directory.")(
          "fastq-offset",
          bpo::value<int>(&fastqOffset_)->default_value(fastqOffset_),
          "FASTQ quality offset value. Set to 33 or 64")(
          "input-qname-suffix-delimiter",
          bpo::value<char>(&inputQnameSuffixDelim_)->default_value(inputQnameSuffixDelim_),
          "Character that delimits input qname suffixes, e.g. / for /1")(
          "Aligner.sw-all",
          bpo::value<int>(&swAll_)->default_value(swAll_),
          "Value of 1 forces smith waterman on all candidate alignments")(
          "Aligner.smith-waterman-method",
          bpo::value<std::string>(&methodSmithWaterman_),
          "Smith Waterman implementation (dragen / mengyao  default = dragen)")(
          "map-only",
          bpo::value<bool>(&mapOnly_)->default_value(mapOnly_),
          "no real alignment, produces alignment information based on seed chains only")(
          "enable-sampling",
          bpo::value<bool>(&samplingEnabled_)->default_value(samplingEnabled_),
          "Automatically detect paired-end parameters by running a sample through the mapper/aligner")(
          "Aligner.pe-stat-mean-insert",
          bpo::value<double>(&alignerPeMeanInsert_)->default_value(alignerPeMeanInsert_),
          "Expected mean of the insert size")(
          "Aligner.pe-stat-stddev-insert",
          bpo::value<double>(&alignerPeStddevInsert_)->default_value(alignerPeStddevInsert_),
          "Expected standard deviation of the insert size")(
          "Aligner.pe-stat-quartiles-insert",
          bpo::value<std::string>(&alignerPeQuartilesInsert_),
          "Q25, Q50, and Q75 quartiles for the insert size")(
          "Aligner.pe-stat-mean-read-len",
          bpo::value<double>(&alignerPeMeanReadLen_)->default_value(alignerPeMeanReadLen_),
          "When setting paired-end insert stats, the expected mean read length")(
          "pe-stats-interval-size",
          bpo::value<uint64_t>(&peStatsIntervalSize_)->default_value(peStatsIntervalSize_),
          "Number of reads to include in each interval of updating paired-end stats")(
          "pe-stats-sample-size",
          bpo::value<uint64_t>(&peStatsSampleSize_)->default_value(peStatsSampleSize_),
          "Number of most recent pairs to include in each update of the paired-end stats")(
          "pe-stats-interval-memory",
          bpo::value<uint8_t>(&peStatsIntervalMemory_)->default_value(peStatsIntervalMemory_),
          "Include reads from up to this many stats intervals in paired-end stats calculations")(
          "pe-stats-interval-delay",
          bpo::value<uint8_t>(&peStatsIntervalDelay_)->default_value(peStatsIntervalDelay_),
          "Number of intervals of lag between sending reads and using resulting stats")(
          "Mapper.filter-len-ratio",
          bpo::value<double>(&mapperFilterLenRatio_)->default_value(mapperFilterLenRatio_),
          "Ratio for controlling seed chain filtering")(
          "Aligner.pe-orientation",
          bpo::value<unsigned>(&alignerPeOrientation_),
          "Expected paired-end orientation: 0=FR, 1=RF, 2=FF")(
          "Aligner.rescue-sigmas",
          bpo::value<double>(&alignerResqueSigmas_)->default_value(alignerResqueSigmas_),
          "For tuning the rescue scan window")(
          "Aligner.rescue-ceil-factor",
          bpo::value<double>(&alignerResqueCeilFactor_)->default_value(alignerResqueCeilFactor_),
          "For tuning the rescue scan window maximum ceiling")(
          "num-threads",
          bpo::value<decltype(mapperNumThreads_)>(&mapperNumThreads_)
              ->default_value(std::thread::hardware_concurrency()),
          "Worker threads for mapper/aligner (default = maximum available on system)")

          ("Aligner.sec-aligns",
           bpo::value<int>(&alignerSecAligns_)->default_value(alignerSecAligns_),
           "Maximum secondary (suboptimal) alignments to report per read")(
              "Aligner.sec-score-delta",
              bpo::value<int>(&alignerSecScoreDelta_)->default_value(alignerSecScoreDelta_),
              "Secondary aligns allowed with pair score no more than this far below primary")(
              "Aligner.sec-phred-delta",
              bpo::value<int>(&alignerSecPhredDelta_)->default_value(alignerSecPhredDelta_),
              "Only secondary alignments with likelihood within this phred of the primary")(
              "Aligner.sec-aligns-hard",
              bpo::value<bool>(&alignerSecAlignsHard_)->default_value(alignerSecAlignsHard_),
              "Set to force unmapped when not all secondary alignments can be output")

              ("preserve-map-align-order",
               bpo::value<bool>(&preserveMapAlignOrder_)->default_value(preserveMapAlignOrder_),
               "Preserve the order of mapper/aligner output to produce deterministic results. Impacts performance")

                  ("build-hash-table",
                   bpo::value<bool>(&buildHashTable_)->default_value(buildHashTable_),
                   "Generate a reference/hash table")(
                      "ht-reference",
                      bpo::value<decltype(htReference_)>(&htReference_),
                      "Reference in FASTA format")

                      ("ht-ref-seed-interval",
                       bpo::value<decltype(htRefSeedInterval_)>(&htRefSeedInterval_),
                       "Number of positions per reference seed")(
                          "ht-seed-len",
                          bpo::value<decltype(htSeedLen_)>(&htSeedLen_)->default_value(htSeedLen_),
                          "Initial seed length to store in hash table")(
                          "ht-max-ext-seed-len",
                          bpo::value<decltype(htMaxExtSeedLen_)>(&htMaxExtSeedLen_)
                              ->default_value(htMaxExtSeedLen_),
                          "Maximum extended seed length")(
                          "ht-max-seed-freq",
                          bpo::value<decltype(htMaxSeedFreq_)>(&htMaxSeedFreq_)
                              ->default_value(htMaxSeedFreq_),
                          "Maximum allowed frequency for a seed match after extension attempts")(
                          "ht-target-seed-freq",
                          bpo::value<decltype(htTargetSeedFreq_)>(&htTargetSeedFreq_)
                              ->default_value(htTargetSeedFreq_),
                          "Target seed frequency for seed extension")(
                          "ht-soft-seed-freq-cap",
                          bpo::value<decltype(htSoftSeedFreqCap_)>(&htSoftSeedFreqCap_)
                              ->default_value(htSoftSeedFreqCap_),
                          "Soft seed frequency cap for thinning")(
                          "ht-max-dec-factor",
                          bpo::value<decltype(htMaxDecFactor_)>(&htMaxDecFactor_)
                              ->default_value(htMaxDecFactor_),
                          "Maximum decimation factor for seed thinning")(
                          "ht-cost-coeff-seed-len",
                          bpo::value<decltype(htCostCoeffSeedLen_)>(&htCostCoeffSeedLen_)
                              ->default_value(htCostCoeffSeedLen_),
                          "Cost coefficient of extended seed length")(
                          "ht-cost-coeff-seed-freq",
                          bpo::value<decltype(htCostCoeffSeedFreq_)>(&htCostCoeffSeedFreq_)
                              ->default_value(htCostCoeffSeedFreq_),
                          "Cost coefficient of extended seed frequency")(
                          "ht-cost-penalty",
                          bpo::value<decltype(htCostPenalty_)>(&htCostPenalty_)
                              ->default_value(htCostPenalty_),
                          "Cost penalty to extend a seed by any number of bases")(
                          "ht-cost-penalty-incr",
                          bpo::value<decltype(htCostPenaltyIncr_)>(&htCostPenaltyIncr_)
                              ->default_value(htCostPenaltyIncr_),
                          "Cost penalty to incrementally extend a seed another step")(
                          "ht-rand-hit-hifreq",
                          bpo::value<decltype(htRandHitHifreq_)>(&htRandHitHifreq_)
                              ->default_value(htRandHitHifreq_),
                          "Include a random hit with each HIFREQ record")(
                          "ht-rand-hit-extend",
                          bpo::value<decltype(htRandHitExtend_)>(&htRandHitExtend_)
                              ->default_value(htRandHitExtend_),
                          "Include a random hit with each EXTEND record of this frequency")(
                          "ht-methylated",
                          bpo::value<decltype(htMethylated_)>(&htMethylated_)->default_value(htMethylated_),
                          "If set to true, generate C->T and G->A converted pair of hashtables")(
                          "ht-num-threads",
                          bpo::value<decltype(htNumThreads_)>(&htNumThreads_),
                          "Worker threads for generating hash table")(
                          "ht-size",
                          bpo::value<decltype(htSize_)>(&htSize_)->default_value(htSize_),
                          "Size of hash table, units B|KB|MB|GB")(
                          "ht-sj-size",
                          bpo::value<decltype(htSjSize_)>(&htSjSize_),
                          "Reserved space for RNA annotated SJs, units B|KB|MB|GB")(
                          "ht-mem-limit",
                          bpo::value<decltype(htMemLimit_)>(&htMemLimit_)->default_value(htMemLimit_),
                          "Memory limit (hash table + reference)")(
                          "ht-decoys",
                          bpo::value<decltype(htDecoys_)>(&htDecoys_)->default_value(htDecoys_),
                          "Path to decoys file (FASTA format)")(
                          "ht-mask-bed",
                          bpo::value<decltype(htMaskBed_)>(&htMaskBed_)->default_value(htMaskBed_),
                          "Bed file for base masking")

      // Hashtable internal parameters
      ("ht-pri-max-seed-freq",
       bpo::value<decltype(htPriMaxSeedFreq_)>(&htPriMaxSeedFreq_)->default_value(htPriMaxSeedFreq_),
       "Maximum frequency for a primary seed match (0 => use maxSeedFreq)")(
          "ht-max-seed-freq-len",
          bpo::value<decltype(htMaxSeedFreqLen_)>(&htMaxSeedFreqLen_)->default_value(htMaxSeedFreqLen_),
          "Ramp from priMaxSeedFreq reaches maxSeedFreq at this seed length")(
          "ht-max-ext-incr",
          bpo::value<decltype(htMaxExtIncr_)>(&htMaxExtIncr_)->default_value(htMaxExtIncr_),
          "Maximum bases to extend a seed by in one step")(
          "ht-ext-table-alloc",
          bpo::value<decltype(htExtTableAlloc_)>(&htExtTableAlloc_)->default_value(htExtTableAlloc_),
          "8-byte records to reserve in extend_table.bin (0=automatic)")(
          "ht-crc-primary",
          bpo::value<decltype(htCrcPrimary_)>(&htCrcPrimary_)->default_value(htCrcPrimary_),
          "Index of CRC polynomial for hashing primary seeds")(
          "ht-crc-extended",
          bpo::value<decltype(htCrcExtended_)>(&htCrcExtended_)->default_value(htCrcExtended_),
          "Index of CRC polynomial for hashing extended seeds")(
          "ht-override-size-check",
          bpo::value<decltype(htOverrideSizeCheck_)>(&htOverrideSizeCheck_)
              ->default_value(htOverrideSizeCheck_),
          "Override hash table size check")(
          "ht-test-only",
          bpo::value<decltype(htTestOnly_)>(&htTestOnly_)->default_value(htTestOnly_),
          "Testing - show user parameters, but don't do anything")(
          "ht-dump-int-params",
          bpo::value<decltype(htDumpIntParams_)>(&htDumpIntParams_)->default_value(htDumpIntParams_),
          "Testing - dump internal parameters")(
          "ht-max-table-chunks",
          bpo::value<decltype(htMaxTableChunks_)>(&htMaxTableChunks_),
          "Maximum ~1GB thread table chunks in memory at once")(
          "ht-repair-strategy",
          bpo::value<decltype(htRepairStrategy_)>(&htRepairStrategy_)->default_value(htRepairStrategy_),
          "Seed extension repair strategy: 0=none, 1=best, 2=rand")(
          "ht-min-repair-prob",
          bpo::value<decltype(htMinRepairProb_)>(&htMinRepairProb_)->default_value(htMinRepairProb_),
          "Minimum probability of success to attempt extended seed repair")(
          "ht-anchor-bin-bits",
          bpo::value<decltype(htAnchorBinBits_)>(&htAnchorBinBits_)->default_value(htAnchorBinBits_),
          "Bits defining reference bins for anchored seed search")(
          "ht-write-hash-bin",
          bpo::value<decltype(htWriteHashBin_)>(&htWriteHashBin_)->default_value(htWriteHashBin_),
          "Write decompressed hash_table.bin and extend_table.bin (0/1)")(
          "ht-uncompress",
          bpo::value<decltype(htUncompress_)>(&htUncompress_)->default_value(htUncompress_),
          "Uncompress hash_table.cmp to hash_table.bin and extend_table.bin (standalone option)")(
          "ht-ext-rec-cost",
          bpo::value<decltype(htExtRecCost_)>(&htExtRecCost_)->default_value(htExtRecCost_),
          "Cost penalty for each EXTEND or INTERVAL record")(
          "ht-max-multi-base-seeds",
          bpo::value<decltype(htMaxMultiBaseSeeds_)>(&htMaxMultiBaseSeeds_),
          "Maximum seeds populated at multi-base codes")

          ("verbose,v", bpo::bool_switch(&verbose_)->default_value(verbose_), "Be talkative");
}

template <typename OptionType>
void checkWarnDeprecatedOption(
    bpo::variables_map& vm,
    const std::string&  officialName,
    OptionType&         official,
    const std::string&  deprecatedName,
    const OptionType&   deprecated)
{
  if (!vm.count(officialName) && vm.count(deprecatedName) && official != deprecated) {
    std::cerr << "WARNING: option '" << deprecatedName
              << "' is deprecated. While still valid, it is recommended to use '" << officialName
              << "' instead." << std::endl;
    official = deprecated;
    return;
  }

  if (vm.count(officialName) && vm.count(deprecatedName) && official != deprecated) {
    BOOST_THROW_EXCEPTION(InvalidOptionException(
        std::string("value specified for '") + officialName + "' conflicts with '" + deprecatedName +
        "'. Please resolve command line option conflict."));
  }
}

void DragenOsOptions::postProcess(bpo::variables_map& vm)
{
  if (vm.count("help") || version_) {
    return;
  }

  if (buildHashTable_ || htUncompress_) {
    return;
  }

  checkWarnDeprecatedOption(
      vm,
      "Aligner.smith-waterman-method",
      methodSmithWaterman_,
      "Aligner.sw-method",
      methodSmithWatermanDeprecated_);

  if (inputFile1_.empty() || boost::filesystem::is_directory(inputFile1_)) {
    BOOST_THROW_EXCEPTION(
        InvalidOptionException("fastq-file1 or bam-input must point to an existing fastq file"));
  }

  if (!inputFile2_.empty() && boost::filesystem::is_directory(inputFile2_)) {
    BOOST_THROW_EXCEPTION(InvalidOptionException("fastq-file2 must point to an existing fastq file"));
  }

  if (!alignerPeQuartilesInsert_.empty()) {
    std::vector<std::string> split;
    boost::split(split, alignerPeQuartilesInsert_, boost::is_space());
    if (3 != split.size()) {
      BOOST_THROW_EXCEPTION(InvalidOptionException(
          "ERROR: 3 values required for configuration parameter pe-stat-quartiles-insert"));
    }
    alignerPeq25_ = std::stoi(split.at(0));
    alignerPeq50_ = std::stoi(split.at(1));
    alignerPeq75_ = std::stoi(split.at(2));
  }

  if (alignerPeMeanInsert_ || alignerPeStddevInsert_ || !alignerPeQuartilesInsert_.empty() ||
      alignerPeMeanReadLen_) {
    if (samplingEnabled_) {
      std::cerr
          << "WARNING: automatic insert-size detection is disabled due to manual override of associated settings.";
    }
    samplingEnabled_ = false;
  }

  if (!outputDirectory_.empty()) {
    if (outputFilePrefix_.empty()) {
      BOOST_THROW_EXCEPTION(InvalidOptionException(
          "ERROR: Output file prefix (--output-file-prefix) is required when --output-directory is given"));
    }
  }

  alnMinScore_ = 22 * matchScore_;
}

}  // namespace options
}  // namespace dragenos
