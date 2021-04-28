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

#ifndef REFERENCE_HASHTABLE_HPP
#define REFERENCE_HASHTABLE_HPP

#include <boost/assert.hpp>
#include <boost/format.hpp>
#include <memory>
#include <vector>

#include "common/Bits.hpp"
#include "reference/Bucket.hpp"
#include "reference/ExtendTableInterval.hpp"
#include "reference/ExtendTableRecord.hpp"
#include "reference/HashRecord.hpp"
#include "reference/HashtableConfig.hpp"
#include "reference/HashtableTraits.hpp"
#include "sequences/CrcHasher.hpp"

namespace dragenos {
namespace reference {

/**
 ** \brief Implementation of a hashtable
 **
 ** Common implementation for versions V7 and V8. The hashtable binary headers are
 ** fully backward compatible (the full binary size is 512 bytes and the new fields
 ** introduced in V8 are set to 0 in V7). The most fundamental difference is the
 ** introduction of the "Extend Table". For V7 it will simply be not present, it will
 ** be set to a nullptr with a reported count of 0 records and it won't be accessed
 ** anyway because the hashtable will not have any record describing any interval set.
 **/
class Hashtable {
public:
  //typedef sequences::Hash Hash;
  typedef uint64_t                   Hash;
  typedef sequences::CrcHasher       CrcHasher;
  typedef reference::HashtableConfig HashtableConfig;
  typedef HashtableTraits            Traits;
  /**
   ** \brief consttructor
   **
   ** \param pointer to the hashtable  co
   **/
  Hashtable(const reference::HashtableConfig* config, const uint64_t* table, const uint64_t* extendTable);
  constexpr static uint64_t getRecordsPerBucket() { return 1UL << Traits::BUCKET_RECORDS_LOG2; }
  constexpr static uint64_t getBytesPerRecord() { return 1UL << Traits::HASH_RECORD_BYTES_LOG2; }
  constexpr static uint64_t getBytesPerBucket() { return getRecordsPerBucket() * getBytesPerRecord(); }
  unsigned                  getPrimaryCrcBits() const { return config_->getPrimaryCrcBits(); }
  unsigned                  getSecondaryCrcBits() const { return config_->getSecondaryCrcBits(); }
  uint64_t getNumberOfBlocks() const { return 1UL << (getPrimaryCrcBits() - getSecondaryCrcBits()); }
  //  uint64_t         getBytesPerBlock() const { return config_->getHashtableBytes() / getNumberOfBlocks(); }
  uint64_t getBytesPerBlock() const
  {
    return (uint32_t)((1 << Traits::MAX_WRAP_BYTES_LOG2) * squeezeFactor_);
  }
  uint64_t         getBucketsPerBlock() const { return getBytesPerBlock() / getBytesPerBucket(); }
  uint32_t         getMinimunFrequencyToExtend() const { return config_->getMinimunFrequencyToExtend(); }
  uint32_t         getMaxSeedFrequency() const { return config_->getMaxSeedFrequency(); }
  const CrcHasher* getPrimaryHasher() const { return primaryHasher_.get(); }
  const CrcHasher* getSecondaryHasher() const { return secondaryHasher_.get(); }
  /**
   ** \brief  calculate the "virtual byte address" of the record for the given hash.
   **
   ** The address is calculated by taking the ADDRESS_START LSB of the hash and multiplying by the
   ** squeeze factor of the hash table.
   ** This is NOT the actual address of the record but the "address" that will then be used to get
   ** the initial bucket and the thread Id for the hash.
   **/
  uint64_t getVirtualByteAddress(const Hash& hash) const
  {
    constexpr static unsigned ADDRESS_START = 19;
    constexpr static unsigned ADDRESS_BITS =
        35;  // up to 35 address bits out of 64 bits - hasher might use less
    return common::bits::getBits<ADDRESS_START, ADDRESS_BITS>(hash) * squeezeFactor_;
  }
  uint64_t getBucketIndex(const uint64_t& virtualByteAddress) const
  {
    return virtualByteAddress >> Traits::HASH_BUCKET_BYTES_LOG2;
  }
  size_t   getHashRecordCount() const { return config_->getHashtableBytes() / sizeof(HashRecord); }
  unsigned getPrimarySeedBases() const { return config_->getPrimarySeedBases(); }
  bool     followChain(const HashRecord& rec, const Hash hash) const;

  /**
   ** \brief process bucket
   **
   ** This is the processing of the initial bucket, before any probing or chaining.
   ** CHAIN_CON and subsequent records must be ignored
   ** A CHAIN_BEG record indicate chaining and is pushed as the last element of hits.
   **/
  bool processInitialBucket(
      const Bucket&                     bucket,
      const Hash&                       hash,
      const uint64_t                    matchBits,
      const uint8_t                     hashThreadId,
      std::vector<HashRecord>&          hits,
      std::vector<ExtendTableInterval>& extenTableIntervals,
      bool                              trace) const;
  /**
   ** \brief when the initial bucket overflows, the default strategy is to store excess hits in the neighbor
   *buckets
   **
   ** The neighborhood is up to 8 consecutive buckets with the caveat that there is a mudulo operation to
   *force
   ** all the buckets to belong to the same hashtable block. These hashtable blocks are identified as having
   *the
   ** same 5 MSB (5 is inferred from the difference between the number of bits from the primary CRC and the
   *number
   ** of bits from the secondary CRC).
   **
   ** \param bucketIndex the index of the initial bucket
   **/
  void probeNeighborBuckets(
      const uint64_t                    initialBucketIndex,
      const uint64_t                    matchBits,
      const uint8_t                     hashThreadId,
      std::vector<HashRecord>&          hits,
      std::vector<ExtendTableInterval>& extenTableIntervals,
      bool                              trace) const;
  /**
   ** \brief probe a single bucket from the probing neighborhood
   **
   ** return true if a relevant record with LF was found. Otherwise returns false.
   **/
  bool probeBucket(
      const Bucket&                     bucket,
      const uint64_t                    matchBits,
      const uint8_t                     hashThreadId,
      std::vector<HashRecord>&          hits,
      std::vector<ExtendTableInterval>& extenTableIntervals,
      bool                              trace) const;
  /**
   ** \brief look for additional in the specified chained bucket
   **
   ** Only records after a CHAIN_CON_* record are taken into consideration.
   ** Completion of the chaining is detected either when the "Last in thread" flag as usual, or when the
   ** MASK or FILTER of the CHAIN_CON_* record returns 0 for the specified hash key. When the mask or filter
   ** returns true, it means that the chain must continue to another bucket, identified by the chain pointer
   ** in the CHAIN_CON_* record. In that case, the CHAIN_CON_* record is stored as the last element of the
   ** hits vector.
   **
   ** \return true if the chaining is complete (either "Last in thread" or "Last in chain" detected
   **/
  bool chainBucket(
      const Bucket&                     bucket,
      const Hash&                       hash,
      const uint64_t                    matchBits,
      const uint8_t                     hashThreadId,
      std::vector<HashRecord>&          hits,
      std::vector<ExtendTableInterval>& extenTableIntervals,
      bool                              trace) const;
  /**
   ** \brief find all relevant HIT and EXTEND records for the given hash key.
   **
   ** The relevant records are further processed to produce only
   ** - regular HIT records that can readily be used by the mapper
   ** - a list of extend table intervals appropriate for various sampling strategies.
   **   in that list, interval at position i includes all intervals a positions j>i.
   **
   ** Prior to that reduction, the relevant records can be:
   ** - Regular HIT records
   ** - Random  Samples (V7 only): HIT records after a HISEQ or EXTEND record
   ** - INTERVAL sets (V8 onwards): a set f 1, 2 or 3 records of one of the INTERVAL_
   **   types describing an interval in the associated extend table
   ** - EXTEND records: used to extend the seed and access a smaller set of relevant  positions
   ** - HIFREQ records (V7 only):
   **
   ** The HIT records related to a specific hash key can be found in several locations:
   ** - the bucket uniquely identified by the hash key
   ** - the buckets within the probing region of the initial bucket
   ** - the the buckets chained from original bucket
   ** - the extended table for extended and high frequency seeds
   ** Chaining and probing are mutually exclusive. Chaining is triggered by the presence of a
   ** CHAIN_BEG_MASK or CHAIN_BEG_LIST record in the initial bucket that matches the query
   ** hash key. Probing is triggered when reaching the end of a bucket without finding any
   ** relevant record with the "Last in thread" flag set and without chaining.
   **/
  void getHits(
      const Hash&                       hash,
      bool                              isExtended,
      std::vector<HashRecord>&          hits,
      std::vector<ExtendTableInterval>& extenTableIntervals,
      bool                              trace = false) const;
  /**
   ** \brief calculate the thread Id of a hash value, for correct matching of the hash records from the
   *hashtable.
   **
   ** The thrad is made of bits [8:3] of the bucket address. The lower 3 bits identify the hash record in the
   *bucket
   ** (the 3 LSB of the address are discarded because a hash record is 8=2^3 bytes). The higher 3 bits enable
   *unique
   ** identification of a hash when probing within a neighborhood of 8 buckets.
   **/
  uint8_t getThreadIdFromVirtualByteAddress(const uint64_t virtualByteAddress) const
  {
    return common::bits::getBits<3, Traits::HASH_RECORD_THREAD_BITS>(virtualByteAddress);
  }

  static uint32_t getListValue(const Hash& hash)
  {
    return common::bits::getBits<0, HashRecord::FILTER_LIST_VALUE_BITS>(hash);
  }
  static uint8_t getMaskValue(const Hash& hash)
  {
    return common::bits::getBits<0, HashRecord::FILTER_MASK_VALUE_BITS>(hash);
  }
  uint64_t getHashBits(const Hash& hash) const
  {
    return common::bits::getBits<0, Traits::HASH_RECORD_HASH_BITS>(hash);
  }
  /// The 30 MSB of the hash record, including thread Id (6 bits), hash bits (23 bits) and EX flag
  uint64_t getMatchBits(const Hash& hash, const bool isExtended) const
  {
    return (((getThreadIdFromVirtualByteAddress(getVirtualByteAddress(hash))
              << Traits::HASH_RECORD_HASH_BITS) |
             getHashBits(hash))
            << 1) |
           static_cast<uint64_t>(isExtended);
  }
  const ExtendTableRecord& getExtendTableRecord(size_t i) const
  {
    return *reinterpret_cast<const ExtendTableRecord*>(extendTable_ + i);
  }

private:
  const HashtableConfig* const config_;
  const uint64_t* const        table_;
  const uint64_t* const        extendTable_;
  const Bucket*                buckets_;
  const size_t                 bucketsCount_;
  std::unique_ptr<CrcHasher>   primaryHasher_;
  std::unique_ptr<CrcHasher>   secondaryHasher_;
  double                       squeezeFactor_;
};

typedef Hashtable PrimaryHashtable;

}  // namespace reference
}  // namespace dragenos

#endif  // ifndef REFERENCE_HASHTABLE_HPP
