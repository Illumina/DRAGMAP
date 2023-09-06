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

#ifndef REFERENCE_HASHTABLE_CONFIG_HPP
#define REFERENCE_HASHTABLE_CONFIG_HPP

#include <array>
#include <cinttypes>
#include <string>
#include <utility>
#include <vector>

#include "common/Exceptions.hpp"
#include "sequences/CrcPolynomial.hpp"

namespace dragenos {

namespace reference {

namespace detail {
// IMPORTANT: the version should only be incremented when the hash table is not
// backwords compatible.  This is used to tell user's when they need to regenerate
// hash tables and by QA to load the correct hash table in their automation framework.
//#define HT_VERSION 7
//static const int HASH_HEADER_SIZE = 512;

// Padding for reference.bin (also used by host software to generate virtual "global"
// position)
#define REF_SEQ_ALIGN_BASES 1024
#define REF_SEQ_MIN_PAD_BASES 256
#define REF_SEQ_END_PAD_BASES 163840

// Hash table types for generation.
typedef enum {
  HT_TYPE_NORMAL = 0,
  HT_TYPE_METHYL_G_TO_A,
  HT_TYPE_METHYL_C_TO_T,
  HT_TYPE_ANCHORED,
  HT_TYPE_NUM_MAX,
} HashtableType;

// Algorithm used to generate digest
typedef enum {
  DIGEST_CRC32  = 0,
  DIGEST_CRC32C = 1,
} DigestType;

// Trie for liftover alignments
#define LIFTOVER_BASE_WIDTH_LOG2 12
#define LIFTOVER_BASE_WIDTH (1 << LIFTOVER_BASE_WIDTH_LOG2)
#define LIFTOVER_CHILDREN_LOG2 2
#define LIFTOVER_CHILDREN (1 << LIFTOVER_CHILDREN_LOG2)
#define LIFT_STATUS_UNALIGNED 0
#define LIFT_STATUS_MATCH 1
#define LIFT_STATUS_INSERT 2
#define LIFT_STATUS_MIXED 3

typedef struct {
  uint64_t pos;
  uint16_t width;
  int8_t   dir;
  uint8_t  status;
} liftoverNode_t;

#define DEFAULT_MEM_SIZE_STR "32GB"
#define MAX_MEM_SIZE_STR "64GB"
#define REF_BASES_PER_IDX_RESERVE_BYTE 64
#define DRAM_FILE_ALIGN_BYTES 256

typedef struct {
  uint32_t hashtableVersion;  // Version number of this hash table
  uint64_t hashtableBytes;    // #bytes in references hash table
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
  // V8 additions
  uint32_t extTabRecs;       // Number of 8-byte records in extend_table.bin
  uint32_t extTabDigest;     // Digest of extend_table.bin
  double   extRecCost;       // Cost penalty for each EXTEND or INTERVAL record
  uint32_t minFreqToExtend;  // Minimum seed frequency eligible for seed extension
  uint8_t  padding[276];     // Reserved for future use
} __attribute__((packed)) hashtableHeader_t;

typedef struct {
  uint64_t seqStart;  // Sequence start offset in reference.bin (untrimmed portion)
  uint32_t begTrim;   // Portion of seqLen leading 'N's trimmed from sequence in reference.bin
  uint32_t endTrim;   // Portion of seqLen trailing 'N's trimmed from sequence in reference.bin
  uint32_t seqLen;    // Reference sequence len (original length, including trimmed portions)
} __attribute__((packed)) hashtableSeq_t;

// For reading version 4 hash_table.cfg.bin files and reporting an error
typedef struct {
  uint64_t seqStart;  // Sequence start offset in reference.bin (untrimmed portion)
  uint32_t seqLen;    // Reference sequence len
} __attribute__((packed)) hashtableSeqv4_t;

struct hashtableConfig_t {
  hashtableHeader_t* hdr;
  hashtableSeq_t**   refSeq;
  char**             seqName;

  int         maxThreads;      // Maximum worker threads
  int         maxGB;           // Maxiumum ~1GB thread table chunks in memory at once
  int         writeHashFile;   // Boolean - write uncompressed hash_table.bin
  int         writeCompFile;   // Boolean - write compressed hash_table.cmp
  const char* sizeStr;         // Size of hash table, units B|KB|MB|GB
  const char* memSizeStr;      // Memory limit (hash table + reference) units B|KB|MB|GB
  const char* sjSizeStr;       // Space to reserve for RNA annotated SJs, NULL/empty for automatic
  uint32_t    methylatedConv;  // For methlated DNA processing, convert G to A or C to T
  int         priPolyIndex;    // Index of CRC polynomial for hashing primary seeds
  int         secPolyIndex;    // Index of CRC polynomial for hashing extended seeds

  char*    refInput;           // Name of FASTA input file
  char*    altLiftover;        // Name of SAM format liftover of alternate contigs in refInput
  char*    configFname;        // Name of hash table configuration file (txt)
  char*    configBinFname;     // Name of hash table configuration file (bin)
  char*    hashFname;          // Name of hash table output file
  char*    compFname;          // Name of compressed hash table output file
  char*    refOutput;          // Name of reference output file
  char*    refIdxFname;        // Name of reference index output file
  char*    repMaskFname;       // Name of repeat mask bitmap output file
  char*    strFname;           // Name of short tandem repeats file
  char*    statsFname;         // Name of hash table stats output file
  char*    decoyFname;         // Name of additional decoy FASTA file
  char*    hostVersion;        // Software version string
  char*    cmdLine;            // Command line used to generate hash table
  int      overrideCheck;      // Override hash tables size check
  int      testOnly;           // Testing - show user parameters, but don't do anything
  int      showIntParams;      // Testing - show internal parameters
  uint8_t* readBuf;            // Buffer used when reading binary config file
  int      usedReadBuf;        // 1 if readBuf used, 0 otherwise
  int      altContigValidate;  // If false, skips hg38/hg19 alt-contigs validation (dragen only)
};

void setDefaultHashParams(hashtableConfig_t* defConfig, const char* dir, HashtableType hashtableType);

//char* generateHashTable(hashtableConfig_t* config, int argc, char* argv[]);

void freeHashParams(hashtableConfig_t* config);

// Query trie for liftover alignments
int liftoverPosition(
    hashtableConfig_t* config, liftoverNode_t* liftover, uint64_t altPos, uint64_t* priPos, int* revComp);
int liftoverSeedIndex(
    hashtableConfig_t* config,
    liftoverNode_t*    liftover,
    uint32_t           altSeedIdx,
    uint32_t*          priSeedIdx,
    int*               revComp);

}  // namespace detail

class HashtableConfig {
public:
  typedef detail::hashtableHeader_t Header;
  // the sequence information and sequence name must be reordered independently - the id_ will be used for
  // this
  struct Sequence : public detail::hashtableSeq_t {
    unsigned id_;

    friend std::ostream& operator<<(std::ostream& os, const Sequence& seq)
    {
      return os << "Sequence(" << seq.id_ << "id " << seq.seqStart << "ss " << seq.seqLen << "sl "
                << seq.begTrim << "bt " << seq.endTrim << "et)";
    }
  };

  typedef std::array<uint64_t, 2> Region;
  /**
   ** \brief constructor from the raw data as it would be in the binary file
   **
   **/
  HashtableConfig(const char* const data, const size_t size);
  uint32_t getHashtableVersion() const { return header_.hashtableVersion; }
  uint64_t getHashtableBytes() const { return header_.hashtableBytes; }  // #bytes in references hash table
  uint64_t getHashtableRecordCount() const { return getHashtableBytes() / 8; }
  uint64_t getHashtableBucketCount() const { return getHashtableBytes() / 64; }
  uint64_t getExtendTableRecordCount() const { return header_.extTabRecs; }
  uint64_t getExtendTableBytes() const { return getExtendTableRecordCount() * 8; }
  uint32_t getMinimunFrequencyToExtend() const { return header_.minFreqToExtend; }
  uint32_t getMaxSeedFrequency() const { return header_.maxSeedFreq; }
  /// total number of bases in the reference sequence file - including padding
  uint64_t                     getReferenceSequenceLength() const { return header_.refSeqLen; }
  unsigned long                getReferenceLength() const;
  std::vector<Region>          getTrimmedRegions() const;
  unsigned                     getNumberOfSequences() const { return header_.numRefSeqs; }
  const std::vector<Sequence>& getSequences() const { return sequences_; }
  // sequence name in the original order
  const std::vector<std::string>& getSequenceNames() const { return sequenceNames_; }
  // sequence name for the sequence at the given offset in sequences_
  const std::string& getSequenceName(const int sequenceOffset) const
  {
    return sequenceNames_.at(sequences_.at(sequenceOffset).id_);
  }
  unsigned long getTableSize64Ths() const { return header_.tableSize64ths; }
  unsigned long getPrimaryCrcBits() const { return header_.priCrcBits; }
  unsigned long getSecondaryCrcBits() const { return header_.secCrcBits; }
  unsigned long getCrcBits(const bool primary) const
  {
    return (primary) ? getPrimaryCrcBits() : getSecondaryCrcBits();
  }
  uint8_t                          getPriCrcPoly(size_t i) const { return header_.priCrcPoly[i]; }
  uint8_t                          getSecCrcPoly(size_t i) const { return header_.secCrcPoly[i]; }
  const uint8_t*                   getPriCrcPoly() const { return header_.priCrcPoly; }
  const uint8_t*                   getSecCrcPoly() const { return header_.secCrcPoly; }
  typedef sequences::CrcPolynomial CrcPolynomial;
  CrcPolynomial getPrimaryPolynomial() const { return CrcPolynomial(getPrimaryCrcBits(), getPriCrcPoly()); }
  CrcPolynomial getSecondaryPolynomial() const
  {
    return CrcPolynomial(getSecondaryCrcBits(), getSecCrcPoly());
  }
  unsigned getPrimarySeedBases() const { return header_.priSeedBases; }
  /// convert a position from the hashtable into reference coordinates (sequence id, offset)
  std::pair<size_t, int64_t> convertToReferenceCoordinates(uint64_t position) const;
  /// get the range of positions for a given sequence
  static std::pair<uint64_t, uint64_t> getPositionRange(const Sequence& sequence)
  {
    const size_t trimmedLength = sequence.seqLen - sequence.begTrim - sequence.endTrim;
    return std::pair<uint64_t, uint64_t>(sequence.seqStart, sequence.seqStart + trimmedLength);
  }
  bool beyondLastCfgSequence(uint64_t position) const
  {
    const auto& sequence = sequences_.back();
    return position > (sequence.seqStart - sequence.begTrim + sequence.seqLen);
  }

private:
  const Header header_;
  // sorted by increasing reference position
  const std::vector<Sequence> sequences_;
  // sorted in the original input order
  const std::vector<std::string> sequenceNames_;
  const std::string              hostVersion_;
  const std::string              commandLine_;
  const std::string              refInput_;
  const std::string              refOutput_;
  const std::string              refIdxFname_;
  const std::string              refFname_;
  const std::string              hashFname_;
  const std::string              altLiftover_;
  // extract the header from the raw data
  const Header* header(const char* const data, const size_t size) const
  {
    if (sizeof(Header) > size) {
      BOOST_THROW_EXCEPTION(
          common::InvalidParameterException("The size provided is smaller than the Header size"));
    }
    const Header* const ret = reinterpret_cast<const Header*>(data);
    if (nullptr == ret) {
      throw common::InvalidParameterException("Failed to contruct HashtableHeader from provided data");
    }
    return ret;
  }
  std::vector<Sequence> sequences(const char* const data, const size_t size) const
  {
    std::vector<Sequence> ret;
    if (sizeof(Header) + header_.numRefSeqs * sizeof(Sequence) > size) {
      BOOST_THROW_EXCEPTION(common::InvalidParameterException(
          "The size provided is smaller than the size of the header and expected reference sequences"));
    }
    const char* current = data + sizeof(header_);
    while (ret.size() < header_.numRefSeqs) {
      detail::hashtableSeq_t tmp = *reinterpret_cast<const detail::hashtableSeq_t*>(current);
      const unsigned         id  = ret.size();
      ret.push_back(Sequence());
      ret.back().seqStart = tmp.seqStart;
      ret.back().begTrim  = tmp.begTrim;
      ret.back().endTrim  = tmp.endTrim;
      ret.back().seqLen   = tmp.seqLen;
      ret.back().id_      = id;
      current += sizeof(detail::hashtableSeq_t);
    }
    std::sort(ret.begin(), ret.end(), [](const Sequence& lhs, const Sequence& rhs) {
      return lhs.seqStart < rhs.seqStart;
    });
    return ret;
  }
  std::vector<std::string> sequenceNames(const char* const data, const size_t size) const
  {
    std::vector<std::string> ret;
    const char* current = data + sizeof(header_) + header_.numRefSeqs * sizeof(detail::hashtableSeq_t);
    while (ret.size() < header_.numRefSeqs) {
      ret.emplace_back(current);
      current += ret.back().size() + 1;
      if ((ssize_t)(size + 1) < (current - data))  // account for positionning 'current' past the '\0'
      {
        BOOST_THROW_EXCEPTION(common::InvalidParameterException(
            "The size provided is smaller than the size of the header and expected reference sequences with their sequence names"));
      }
    }
    assert(ret.size() == sequences_.size());
    return ret;
  }
  std::string refIdxFname(const char* const data, const size_t /*size*/) const
  {
    size_t bufOffset = sizeof(header_) + header_.numRefSeqs * sizeof(detail::hashtableSeq_t);
    for (uint32_t i = 0; i < header_.numRefSeqs; ++i) {
      bufOffset += sequenceNames_[i].size() + 1;
    }
    const char* hostVersion = data + bufOffset;
    bufOffset += strlen(hostVersion) + 1;
    const char* cmdLine = data + bufOffset;
    bufOffset += strlen(cmdLine) + 1;
    const char* refInput = data + bufOffset;
    bufOffset += strlen(refInput) + 1;
    const char* refOutput = data + bufOffset;
    bufOffset += strlen(refOutput) + 1;
    const char* refIdxFname = data + bufOffset;
    return std::string(refIdxFname);
  }
  const char* stringAddress(const char*, const std::string*)
  {
    static const char* m = "NOT IMPLEMENTED";
    return m;
  }
};

}  // namespace reference
}  // namespace dragenos

#endif  // #ifndef REFERENCE_HASHTABLE_CONFIG_HPP
