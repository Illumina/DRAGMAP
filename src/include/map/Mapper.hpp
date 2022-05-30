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

#ifndef MAP_MAPPER_HPP
#define MAP_MAPPER_HPP

#include <bitset>
#include <vector>

#include <boost/format.hpp>

#include "BestIntervalTracker.hpp"
#include "common/Exceptions.hpp"
#include "reference/HashRecord.hpp"
#include "reference/Hashtable.hpp"
#include "sequences/Read.hpp"

namespace dragenos {
namespace map {

/**
 ** \brief A class that can find relevant chains of positions for reads
 **
 **/
class Mapper {
public:
  typedef reference::HashRecord          HashRecord;
  typedef reference::Hashtable           Hashtable;
  typedef reference::ExtendTableInterval ExtendTableInterval;
  typedef reference::ExtendTableRecord   ExtendTableRecord;
  typedef sequences::Read                Read;
  typedef sequences::Seed                Seed;
  static constexpr unsigned              MAX_EXTENSION_STEP = 12;
  static constexpr unsigned EXTENSION_ID_BIN_SHIFT = HashRecord::EXTENSION_ID_BITS + 2 * MAX_EXTENSION_STEP;
  static constexpr unsigned MAX_HIFREQ_HITS        = 16;

  explicit Mapper(const Hashtable* hashtable)
    : hashtable_(hashtable),
      extensionIdBinMask_(generateExtensionIdBinMask(hashtable)),
      addressSegmentMask_(generateAddressSegmentMask(hashtable))
  {
  }
  /**
   ** \brief find all the relevant chains of positions for the given read
   **
   ** A position is a seed position and implemented as a pair (Seed, Position).
   ** A chain of positions is a range of positions and is implemented as a pair
   ** providing the begin and end iterator of the range.
   **
   ** The position chains are returned via an output parameter instead of a return
   ** value because this method will be called at a high frequency and it is worth
   ** letting the client code managing the output buffer to avoid spurious memory
   ** allocations.
   **
   **/
  void getPositionChains(const Read& read, ChainBuilder& chainBuilder) const;
  /**
   ** \brief find all the relevant hash records for a given seed and add them to the
   ** position chains
   **
   ** This is the method implementing the management of the seeds and their extensions.
   ** It is responsible for coordinating the generation of the hash for the primary seed
   ** and for all the extensions. This method will also delegate the call to the
   ** associated hashtable to retrieve the relevant hash records. Finally, it will handle
   ** appropriately the different types of records and ultimately dispatch all the hit
   ** records to the chain builder.
   **/
  void addToPositionChains(
      const Seed&                       seed,
      ChainBuilder&                     chainBuilder,
      std::vector<BestIntervalTracker>& globalBestIntvls,
      uint32_t&                         intvl_non_sample_longest,
      uint32_t&                         num_extension_failure) const;
  void addRandomSamplesToPositionChains(
      const Seed&                    seed,
      const bool                     seedIsReverseComplement,
      const unsigned                 halfExtension,
      const std::vector<HashRecord>& hashRecords,
      ChainBuilder&                  chainBuilder) const;
  void getRandomSamplesFromMatchInterval(
      const Seed&                     seed,
      const uint32_t                  sampleSize,
      const uint32_t                  intvl_start,
      const uint32_t                  intvl_len,
      std::vector<ExtendTableRecord>& hashRecords) const;
  void addExtraIntervalSamplesToPositionChains(
      const BestIntervalTracker& globalBestIntvl, ChainBuilder& chainBuilder) const;
  /**
   ** \brief Center the correctly oriented wings into a 24 bits vector (up to 6 bases wings) and concatenate
   *with the extension id
   **
   ** \param seed the seed to extend
   ** \param hash the hash to use to generate the extension ID bin
   ** \param extendRecord the hash record that triggered the extension
   ** \param fromHalfExtension the size of each extension wing on the previous query
   ** \param seedIsReverseComplement flag indicating if the extension wings should be reverse complemented
   **/
  uint64_t getExtendedKey(
      const Seed&       seed,
      const uint64_t    hash,
      const HashRecord& extendRecord,
      const unsigned    fromHalfExtension,
      const bool        seedIsReverseComplement) const;
  /// Select the primary hash MSB that don't overlap bit positions from secondary hashes
  static uint64_t generateAddressSegmentMask(const Hashtable* hashtable)
  {
    const uint64_t primaryMask   = (static_cast<uint64_t>(1) << hashtable->getPrimaryCrcBits()) - 1;
    const uint64_t secondaryMask = (static_cast<uint64_t>(1) << hashtable->getSecondaryCrcBits()) - 1;
    return primaryMask ^ secondaryMask;
  }
  /// mask used to select the correct number of LSB from the previous hash to extract the extension Id bin
  static uint64_t generateExtensionIdBinMask(const Hashtable* hashtable)
  {
    const unsigned secondaryCrcBits   = hashtable->getSecondaryCrcBits();
    const unsigned extensionIdBinBits = secondaryCrcBits > EXTENSION_ID_BIN_SHIFT ? secondaryCrcBits - 42 : 0;
    return (static_cast<uint64_t>(1) << extensionIdBinBits) - 1;
  }
  /// extract the extension id bin from the hash and shift it to the correct bit position
  uint64_t getShiftedExtensionIdBin(const uint64_t hash) const
  {
    return (hash & extensionIdBinMask_) << EXTENSION_ID_BIN_SHIFT;
  }
  const Hashtable* getHashtable() const { return hashtable_; }
  /**
   ** \brief output the seed chains without mapping
   **
   ** The information is stuffed into the CIGAR, according to this format:
   ** Header (64 bytes)
   ** [2]	Exact seed attempts
   ** [2]	Edited seed attempts
   ** [2]	Primary accesses:  first
   ** [2]	Primary accesses:  probe
   ** [2]	Primary accesses:  chain
   ** [2]	Secondary accesses:  first
   ** [2]	Secondary accesses:  probe
   ** [2]	Secondary accesses:  chain
   ** [2]	Primary results:  miss
   ** [2]	Primary results:  hit
   ** [2]	Primary results:  hi-freq
   ** [2]	Primary results:  extend
   ** [2]	Secondary results:  miss
   ** [2]	Secondary results:  hit
   ** [2]	Largest miss interval
   ** [2]	Longest seed extension
   ** [2]	Sum of seed extensions
   ** [2]	Seeds with 1 primary hit
   ** [2]	Seeds with 2 primary hits
   ** [2]	Seeds with 3 primary hits
   ** [2]	Seeds with primary freq. 4+
   ** [2]	Seeds with primary freq. 6+
   ** [2]	Seeds with primary freq. 8+
   ** [2]	Seeds with primary freq. 12+
   ** [2]	Seeds with primary freq. 16+
   ** [2]	Seeds with primary freq. 24+
   ** [2]	Seeds with primary freq. 32+
   ** [2]	Seeds with primary freq. 48+
   ** [2]	Seeds with primary freq. 64+
   ** [2]	Seeds with primary freq. 96+
   ** [2]	Seeds with primary freq. 128+
   ** [2]	Chain count
   ** Chain Record (16 bytes)
   ** [4]	Start position in reference (chromosome lookup not performed)
   ** [4]	End position in reference (End > Start)
   ** [2]	Start position in read
   ** [2]	End position in read (End > Start)
   ** [2]	Seed count
   ** [2]	Flags:
   ** 0x1 = Reverse-complemented
   ** 0x2 = Chain filtered out
   ** 0x4 = Edited seed present
   ** 0x8 = Perfect alignment
   ** 0x80 = (simulation hack) REF_ID in Flags[15:8]; reference positions relative
   **
   **/
  void generateMapperCigar(std::ostream& os, ChainBuilder& chainBuilder) const;

private:
  const Hashtable* hashtable_;
  const uint64_t   extensionIdBinMask_;
  /// mask used to keep all related hashes (probing region, chains, extensions) in the same memory segment
  const uint64_t addressSegmentMask_;

  /// local vectors used in addToPositionChains, member variables for malloc optimization
  mutable std::vector<HashRecord>          hashRecords_;
  mutable std::vector<ExtendTableInterval> extendTableIntervals_;
};

}  // namespace map
}  // namespace dragenos

#endif  // #ifndef MAP_MAPPER_HPP
