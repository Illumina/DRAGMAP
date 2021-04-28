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

#include "reference/Hashtable.hpp"

namespace dragenos {
namespace reference {

Hashtable::Hashtable(const HashtableConfig* config, const uint64_t* table, const uint64_t* extendTable)
  : config_(config),
    table_(table),
    extendTable_(extendTable),
    buckets_(reinterpret_cast<const Bucket*>(table_)),
    bucketsCount_(config_->getHashtableBytes() / sizeof(Bucket)),
    primaryHasher_(new CrcHasher(config_->getPrimaryPolynomial())),
    secondaryHasher_(new CrcHasher(config_->getSecondaryPolynomial())),
    squeezeFactor_(static_cast<double>(config_->getTableSize64Ths()) / 64)
{
}

bool Hashtable::processInitialBucket(
    const Bucket&                     bucket,
    const Hash&                       hash,
    const uint64_t                    matchBits,
    const uint8_t                     hashThreadId,
    std::vector<HashRecord>&          hits,
    std::vector<ExtendTableInterval>& extendTableIntervals,
    const bool                        trace) const
{
  // TODO: check if it is required to have at least one hash record in the thread Id to start probing
  bool       lastInThread = true;
  bool       fullBucket   = true;  // probably not useful, other than sanity check
  bool       chaining     = false;
  HashRecord chainBeginRecord;
  for (auto& hashRecord : bucket) {
    //////////
    if (trace) {
      std::cerr << "    Hashtable::processInitialBucket: record: " << std::hex << std::setw(8)
                << std::setfill('0') << (hashRecord.getValue() >> 32) << " " << std::setw(8)
                << (hashRecord.getValue() & 0xFFFFFFFF) << " : " << std::setw(2)
                << ((hashRecord.getValue() >> 24) & 0xFF) << " hashThreadId: " << (unsigned)hashThreadId
                << " hashRecord thread Id: " << std::setw(2) << (unsigned)hashRecord.getThreadId()
                << std::setfill(' ') << " LF: " << hashRecord.isLastInThread() << std::dec
                << " recordType: " << (int)hashRecord.getType() << std::endl;
    }
    //////////

    const auto type = hashRecord.getType();
    if ((HashRecord::HIT == type) || (HashRecord::HIFREQ == type) || (HashRecord::EXTEND == type) ||
        (HashRecord::INTERVAL_SL == type) || (HashRecord::INTERVAL_SLE == type) ||
        (HashRecord::INTERVAL_S == type) || (HashRecord::INTERVAL_L == type)) {
      if (hashThreadId == hashRecord.getThreadId()) {
        if (matchBits == hashRecord.getMatchBits()) {
          hits.push_back(hashRecord);
        }
        lastInThread = hashRecord.isLastInThread();
        if (lastInThread) {
          break;
        }
      }
      continue;
    } else if ((HashRecord::CHAIN_BEG_MASK == type) || (HashRecord::CHAIN_BEG_LIST == type)) {
      if (followChain(hashRecord, hash)) {
        chaining         = true;
        lastInThread     = false;
        chainBeginRecord = hashRecord;
      }
      continue;
    } else if (HashRecord::EMPTY == type) {
      fullBucket = false;  // probably not useful because there should be LF for all threads in the bucket
      // in principle this would imply no probing and no chaining;
      // NOTE: there could be more records after the EMPTY record therefore we can't break
      continue;
    } else if ((HashRecord::CHAIN_CON_MASK == type) || (HashRecord::CHAIN_CON_LIST == type)) {
      // all subsequent records are only relevant to chaining into this bucket
      break;
    } else if (HashRecord::REPAIR == type) {
      BOOST_THROW_EXCEPTION(std::invalid_argument("REPAIR type obsolete for Hash Records"));
    }
    boost::format message =
        boost::format("Unknown Hash Record type: %x: record: %x") % type % hashRecord.getValue();
    BOOST_THROW_EXCEPTION(std::invalid_argument(message.str()));
  }
  // if there are EMPTY records then there should be no chaining and no probing
  BOOST_ASSERT(fullBucket || lastInThread);
  BOOST_ASSERT(fullBucket || !chaining);
  // if lastInChain then chaining should not be possible
  BOOST_ASSERT(!(lastInThread && chaining));
  if (chaining) {
    BOOST_ASSERT(!lastInThread);
    hits.push_back(chainBeginRecord);
  }
  return lastInThread;
}

bool Hashtable::probeBucket(
    const Bucket&                     bucket,
    const uint64_t                    matchBits,
    const uint8_t                     hashThreadId,
    std::vector<HashRecord>&          hits,
    std::vector<ExtendTableInterval>& extendTableIntervals,
    const bool                        trace) const
{
  auto hashRecord = bucket.cbegin();
  // alternate exit on relevant hashRecord->isLastInThread()
  while (hashRecord != bucket.cend() && (HashRecord::CHAIN_CON_MASK != hashRecord->getType()) &&
         (HashRecord::CHAIN_CON_LIST != hashRecord->getType())) {
    //////////
    if (trace) {
      std::cerr << "    Hashtable::processInitialBucket(probeBucket): record: " << std::hex << std::setw(8)
                << std::setfill('0') << (hashRecord->getValue() >> 32) << " " << std::setw(8)
                << (hashRecord->getValue() & 0xFFFFFFFF) << " : " << std::setw(2)
                << ((hashRecord->getValue() >> 24) & 0xFF) << " hashThreadId: " << (unsigned)hashThreadId
                << " hashRecord thread Id: " << std::setw(2) << (unsigned)hashRecord->getThreadId()
                << std::setfill(' ') << " LF: " << hashRecord->isLastInThread() << std::dec
                << " recordType: " << (int)hashRecord->getType() << std::endl;
    }
    //////////

    const auto type = hashRecord->getType();
    if ((HashRecord::HIT == type) || (HashRecord::HIFREQ == type) || (HashRecord::EXTEND == type) ||
        (HashRecord::INTERVAL_SL == type) || (HashRecord::INTERVAL_SLE == type) ||
        (HashRecord::INTERVAL_S == type) || (HashRecord::INTERVAL_L == type)) {
      if (hashThreadId == hashRecord->getThreadId()) {
        if (matchBits == hashRecord->getMatchBits()) {
          hits.push_back(*hashRecord);
        }
        if (hashRecord->isLastInThread()) {
          return true;
        }
      }
    }
    ++hashRecord;
  }
  return false;
}

void Hashtable::probeNeighborBuckets(
    const uint64_t                    initialBucketIndex,
    const uint64_t                    matchBits,
    const uint8_t                     hashThreadId,
    std::vector<HashRecord>&          hits,
    std::vector<ExtendTableInterval>& extendTableIntervals,
    const bool                        trace) const
{
  // probing is forced to be constrained to a single block with a modulo operation
  const uint64_t blockStartBucketIndex = initialBucketIndex - (initialBucketIndex % getBucketsPerBlock());

  if (trace)
    std::cerr << std::hex << " bucket address: " << &buckets_[blockStartBucketIndex]
              << "\tbucket index: " << blockStartBucketIndex << std::dec
              << "\tgetBucketsPerBlock() :" << getBucketsPerBlock() << std::endl;
  // In practice, should exit on an LF=true record
  for (unsigned i = 1; Traits::MAX_PROBES > i; ++i) {
    if (trace) std::cerr << " probe cnt:" << i << std::endl;
    const uint64_t currentbucketIndex =
        blockStartBucketIndex + ((initialBucketIndex + i) % getBucketsPerBlock());
    if (probeBucket(
            buckets_[currentbucketIndex], matchBits, hashThreadId, hits, extendTableIntervals, trace)) {
      return;
    }
  }
  // TODO: check if the LF is expected to be always set when probing
  // BOOST_THROW_EXCEPTION(std::invalid_argument("Probing completed through all neighbor buckets without finding any hash record with LF=true"));
}

bool Hashtable::followChain(const HashRecord& record, const Hash hash) const
{
  const auto recordType = record.getType();
  if ((recordType == HashRecord::CHAIN_BEG_MASK) || (recordType == HashRecord::CHAIN_CON_MASK)) {
    // Check if the correct bit is turned on in the filter for the given hash
    const uint32_t filterMask = record.getFilterMask();
    return ((filterMask >> getMaskValue(hash)) & 1);
  }
  if ((recordType == HashRecord::CHAIN_BEG_LIST) || (recordType == HashRecord::CHAIN_CON_LIST)) {
    // Check if the correct target-hash bits are present in at least one of the CHAIN_LIST sub-records
    const uint32_t listValue  = getListValue(hash);
    uint32_t       filterList = record.getFilterMask();
    for (unsigned i = 0; i < HashRecord::FILTER_LIST_COUNT; i++) {
      if (listValue == common::bits::getBits<0, HashRecord::FILTER_LIST_VALUE_BITS>(filterList)) {
        return true;
      }
      filterList = (filterList >> HashRecord::FILTER_LIST_VALUE_BITS);
    }
    return false;
  }
  BOOST_ASSERT_MSG(false, "followChain invoked on a record that is not a chain record");
  return false;
}

bool Hashtable::chainBucket(
    const Bucket&                     bucket,
    const Hash&                       hash,
    const uint64_t                    matchBits,
    const uint8_t                     hashThreadId,
    std::vector<HashRecord>&          hits,
    std::vector<ExtendTableInterval>& extendTableIntervals,
    const bool                        trace) const
{
  //bool fullBucket = true;
  auto hashRecord = bucket.cbegin();
  // skip all hash records until the first CHAIN_CON_{MASK,LIST}
  while ((hashRecord != bucket.cend()) && !hashRecord->isChainCon()) {
    //////////
    if (trace) {
      std::cerr << "    Hashtable::chainBucket: discarding: record: " << std::hex << std::setw(8)
                << std::setfill('0') << (hashRecord->getValue() >> 32) << " " << std::setw(8)
                << (hashRecord->getValue() & 0xFFFFFFFF) << " : " << std::setw(2)
                << ((hashRecord->getValue() >> 24) & 0xFF) << " hashThreadId: " << (unsigned)hashThreadId
                << " hashRecord thread Id: " << std::setw(2) << (unsigned)hashRecord->getThreadId()
                << std::setfill(' ') << " LF: " << hashRecord->isLastInThread() << std::dec
                << " recordType: " << (int)hashRecord->getType() << std::endl;
    }
    //////////

    // if (HashRecord::EMPTY == hashRecord->getType())
    //{
    //  fullBucket = false; // probably not useful because there should be LF for all threads in the bucket
    //}
    ++hashRecord;
  }
  if (hashRecord == bucket.cend()) {
    BOOST_THROW_EXCEPTION(std::invalid_argument("Failed to fing CHAIN_CON record in bucket"));
  }

  //////////
  if (trace) {
    std::cerr << "    Hashtable::chainBucket: CHAIN_CON record: " << std::hex << std::setw(8)
              << std::setfill('0') << (hashRecord->getValue() >> 32) << " " << std::setw(8)
              << (hashRecord->getValue() & 0xFFFFFFFF) << " : " << std::setw(2)
              << ((hashRecord->getValue() >> 24) & 0xFF) << " hashThreadId: " << (unsigned)hashThreadId
              << " hashRecord thread Id: " << std::setw(2) << (unsigned)hashRecord->getThreadId()
              << std::setfill(' ') << " LF: " << hashRecord->isLastInThread() << std::dec
              << " recordType: " << (int)hashRecord->getType() << std::endl;
  }
  //////////

  const auto chainConRecord = *hashRecord;
  // if this is the last bucket in the chain for the hash, it will also be the Last in Thread
  bool lastInThread = false;
  ++hashRecord;
  while (hashRecord != bucket.cend()) {
    //////////
    if (trace) {
      std::cerr << "    Hashtable::chainBucket: considering record: " << std::hex << std::setw(8)
                << std::setfill('0') << (hashRecord->getValue() >> 32) << " " << std::setw(8)
                << (hashRecord->getValue() & 0xFFFFFFFF) << " : " << std::setw(2)
                << ((hashRecord->getValue() >> 24) & 0xFF) << " hashThreadId: " << (unsigned)hashThreadId
                << " hashRecord thread Id: " << std::setw(2) << (unsigned)hashRecord->getThreadId()
                << std::setfill(' ') << " LF: " << hashRecord->isLastInThread() << std::dec
                << " recordType: " << (int)hashRecord->getType() << std::endl;
    }
    //////////

    const auto type = hashRecord->getType();
    if ((HashRecord::HIT == type) || (HashRecord::HIFREQ == type) || (HashRecord::EXTEND == type) ||
        (HashRecord::INTERVAL_SL == type) || (HashRecord::INTERVAL_SLE == type) ||
        (HashRecord::INTERVAL_S == type) || (HashRecord::INTERVAL_L == type)) {
      if (hashThreadId == hashRecord->getThreadId()) {
        if (matchBits == hashRecord->getMatchBits()) {
          hits.push_back(*hashRecord);
        }
        if (hashRecord->isLastInThread()) {
          lastInThread = true;
          break;
        }
      }
    }
    // else if (HashRecord::EMPTY == type)
    //{
    //  fullBucket = false; // probably not useful because there should be LF for all threads in the bucket
    //}
    else if (HashRecord::EMPTY != type)
    // else
    {
      boost::format message = boost::format("Unexpected record type while chaining: %1x: %8x %8x") % type %
                              (hashRecord->getValue() >> 32) % (hashRecord->getValue() & 0xFFFFFFFF);
      BOOST_THROW_EXCEPTION(std::invalid_argument(message.str()));
    }
    ++hashRecord;
  }
  // TODO: verify consistency between lastInThread, fullBucket and followChain
  //BOOST_ASSERT(fullBucket || lastInThread);
  // if (lastInThread && followChain(chainConRecord, hash))
  //{
  //  BOOST_THROW_EXCEPTION(std::invalid_argument("Both 'Last in thread' and 'Follow chain' are true"));
  //}
  // if (followChain(chainConRecord, hash))
  if ((!lastInThread) && followChain(chainConRecord, hash)) {
    hits.push_back(chainConRecord);
    return false;
  } else {
    return true;
  }
}

void Hashtable::getHits(
    const Hash&                       hash,
    const bool                        isExtended,
    std::vector<HashRecord>&          hits,
    std::vector<ExtendTableInterval>& extendTableIntervals,
    const bool                        trace) const
{
  hits.clear();
  const auto matchBits          = getMatchBits(hash, isExtended);
  const auto virtualByteAddress = getVirtualByteAddress(hash);
  const auto bucketIndex        = getBucketIndex(virtualByteAddress);
  const auto hashThreadId       = getThreadIdFromVirtualByteAddress(virtualByteAddress);
  hits.clear();

  if (trace)
    std::cerr << std::hex << "hash: " << hash << " bucket index: " << bucketIndex << std::dec
              << " initial interval count: " << extendTableIntervals.size() << std::endl;

  bool lastInThread = processInitialBucket(
      buckets_[bucketIndex], hash, matchBits, hashThreadId, hits, extendTableIntervals, trace);

  if (trace)
    std::cerr << "found hits:" << hits.size() << " interval count: " << extendTableIntervals.size()
              << " lastInThread: " << lastInThread << std::endl;

  if (!lastInThread) {
    const bool probing = hits.empty() || (!hits.back().isChainBegin());
    if (probing) {
      probeNeighborBuckets(bucketIndex, matchBits, hashThreadId, hits, extendTableIntervals, trace);
    } else /* chaining */
    {
      const auto baseBucketIndex = (bucketIndex >> HashRecord::CHAIN_POINTER_BITS)
                                   << HashRecord::CHAIN_POINTER_BITS;
      while (!lastInThread) {
        BOOST_ASSERT(!hits.empty());
        const HashRecord chainingRecord = hits.back();
        hits.pop_back();
        BOOST_ASSERT(chainingRecord.isChainRecord());
        const uint64_t chainPointer = chainingRecord.getChainPointer();
        lastInThread                = chainBucket(
            buckets_[baseBucketIndex + chainPointer],
            hash,
            matchBits,
            hashThreadId,
            hits,
            extendTableIntervals,
            trace);
      }
    }
  }
  // TODO: convert interval sets into extend table intervals, if any
  // ASSUMPTION: the interval records, if any are at the back

  ////////////////////
  if (trace) std::cerr << " final hit count: " << hits.size() << std::endl;
  ////////////////////

  auto begin = hits.end();
  while ((hits.begin() != begin) && ((HashRecord::INTERVAL_SL == (begin - 1)->getType()) ||
                                     (HashRecord::INTERVAL_SLE == (begin - 1)->getType()) ||
                                     (HashRecord::INTERVAL_S == (begin - 1)->getType()) ||
                                     (HashRecord::INTERVAL_L == (begin - 1)->getType()))) {
    --begin;
  }
  if (hits.end() != begin) {
    extendTableIntervals.push_back(ExtendTableInterval(begin, hits.end()));
    hits.erase(begin, hits.end());
  }
}

}  // namespace reference
}  // namespace dragenos
