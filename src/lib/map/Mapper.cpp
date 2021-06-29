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

#include <boost/format.hpp>
#include <boost/range/adaptor/reversed.hpp>
#include <numeric>
#include <unordered_set>

#include "common/Crc32Hw.hpp"
#include "common/DragenLogger.hpp"
#include "map/Mapper.hpp"

namespace dragenos {
namespace map {

// debug variables
const char* typeString[17] = {"EMPTY",
                              "HIFREQ",
                              "EXTEND",
                              "REPAIR",  // obsolete
                              "CHAIN_BEG_MASK",
                              "CHAIN_BEG_LIST",
                              "CHAIN_CON_MASK",
                              "CHAIN_CON_LIST",
                              "INTERVAL_SL",
                              "INTERVAL_SLE",
                              "INTERVAL_S",
                              "INTERVAL_L",
                              "",
                              "",
                              "",
                              "",
                              "HIT"};
// end of debug variables
#ifdef TRACE_SEED_CHAINS
static uint32_t seedOffset = 0;
#endif
void Mapper::getPositionChains(const Read& read, ChainBuilder& chainBuilder) const
{
  chainBuilder.clear();
  const unsigned seedLength = hashtable_->getPrimarySeedBases();
  chainBuilder.setFilterConstant(seedLength);
  const unsigned                   readLength         = read.getLength();
  constexpr int32_t                SEED_PERIOD        = 2;
  constexpr uint32_t               SEED_PATTERN       = 0x01;
  constexpr uint8_t                FORCE_LAST_N_SEEDS = 0;  //1;
  std::vector<BestIntervalTracker> globalBestIntvls;
  globalBestIntvls.reserve(readLength);
  uint32_t num_extension_failure      = 0;
  uint32_t longest_nonsample_seed_len = 0;
  uint32_t num_non_sample_seed_chains = 0;
  // TODO: check the cost of the underlying memory allocations and cace the seed positions buffer if needed
  const auto seedOffsets =
      Seed::getSeedOffsets(readLength, seedLength, SEED_PERIOD, SEED_PATTERN, FORCE_LAST_N_SEEDS);
#ifdef TRACE_SEED_CHAINS
  for (size_t i = 0; i != seedOffsets.size(); ++i) {
    seedOffset = seedOffsets[i];
#else
  for (const auto& seedOffset : seedOffsets) {
#endif

#ifdef TRACE_SEED_CHAINS
    std::cerr << "\n------------------------\nMapper::getPositionChains: seed offset: " << seedOffset
              << std::endl;
    std::cerr << "--------------------------longest_nonsample_seed_len:" << longest_nonsample_seed_len
              << std::endl;
#endif

    if (Seed::isValid(read, seedOffset, seedLength)) {
      const sequences::Seed seed(&read, seedOffset, seedLength);
      assert(
          seed.isValid(0));  // getSeedOffset is supposed to produce offsets only for valid non-extended seeds
      addToPositionChains(
          seed, chainBuilder, globalBestIntvls, longest_nonsample_seed_len, num_extension_failure);
    } else {
#ifdef TRACE_SEED_CHAINS
      std::cerr << "Seed validation failed(either contains N or longer than read length)." << std::endl;
#endif
    }
  }
  // random sampling from extra interval
  if (!globalBestIntvls.empty()) {
    BestIntervalTracker globalBestIntvl = globalBestIntvls.front();

    for (const auto& item : globalBestIntvls) {
      if (globalBestIntvl.isWorseThan(item)) {
#ifdef TRACE_SEED_CHAINS
        std::cerr << "Replace globalBestIntvl:\t" << globalBestIntvl.getSeed()
                  << "\tstart:" << globalBestIntvl.getStart() << "\tlength:" << globalBestIntvl.getLength()
                  << "\tseedLenght:" << globalBestIntvl.getSeedLength() << "\twith\t" << item.getSeed()
                  << "\tstart:" << item.getStart() << "\tlength:" << item.getLength()
                  << "\tseedLength:" << item.getSeedLength() << std::endl;
#endif

        globalBestIntvl = item;
      }
    }

    num_non_sample_seed_chains = std::count_if(chainBuilder.begin(), chainBuilder.end(), [](SeedChain item) {
      return not item.hasOnlyRandomSamples();
    });
    if (globalBestIntvl.isValidExtra(num_non_sample_seed_chains, longest_nonsample_seed_len)) {
#ifdef TRACE_SEED_CHAINS
      std::cerr << "Sampling from global best interval:\t" << globalBestIntvl.getStart() << ":"
                << globalBestIntvl.getLength()
                << "\tlongest non-sample seed length:" << longest_nonsample_seed_len << std::endl;
#endif
      addExtraIntervalSamplesToPositionChains(globalBestIntvl, chainBuilder);
    } else {
#ifdef TRACE_SEED_CHAINS
      std::cerr << "Had global best intervals but not valid extra:\t" << globalBestIntvl.getStart() << ":"
                << globalBestIntvl.getLength()
                << "\tlongest non-sample seed length:" << longest_nonsample_seed_len << std::endl;
#endif
    }
  }
  chainBuilder.filterChains();
}

void Mapper::addRandomSamplesToPositionChains(
    const Seed&                    seed,
    const bool                     seedIsReverseComplement,
    const unsigned                 halfExtension,
    const std::vector<HashRecord>& hashRecords,
    ChainBuilder&                  chainBuilder) const
{
  const bool isRandomSample = true;
  // skip the first record as it is either HIFREQ or EXTEND
  for (auto record = hashRecords.begin() + 1; hashRecords.end() != record; ++record) {
    if (record->isHit()) {
      const bool orientation = (seedIsReverseComplement ^ record->isReverseComplement());
      chainBuilder.addSeedPosition(
          SeedPosition(seed, record->getPosition(), halfExtension), orientation, isRandomSample);
    } else {
      boost::format message =
          boost::format("Expected HIT hash record but got record of type: %i") % record->getType();
      BOOST_THROW_EXCEPTION(common::PostConditionException(message.str()));
    }
  }
}

uint64_t Mapper::getExtendedKey(
    const Seed&       seed,
    const uint64_t    hash,
    const HashRecord& extendRecord,
    const unsigned    fromHalfExtension,
    const bool        seedIsReverseComplement) const
{
  const auto     shiftedExtensionIdBin = getShiftedExtensionIdBin(hash);
  const unsigned halfExtension         = extendRecord.getExtensionLength() / 2;
  const unsigned toHalfExtension       = fromHalfExtension + halfExtension;
  const unsigned halfPaddingBits       = (6 - halfExtension) * 2;
  const auto extendBases = seed.getExtendedData(fromHalfExtension, toHalfExtension, seedIsReverseComplement)
                           << halfPaddingBits;
  const auto shiftedExtensionId = extendRecord.getExtensionId() << (2 * MAX_EXTENSION_STEP);
  return shiftedExtensionIdBin | shiftedExtensionId | extendBases;
}

void Mapper::addToPositionChains(
    const Seed&                       seed,
    ChainBuilder&                     chainBuilder,
    std::vector<BestIntervalTracker>& globalBestIntvls,
    uint32_t&                         longest_nonsample_seed_len,
    uint32_t&                         num_extension_failure) const
{
  //////////
  //std::cerr << "Mapper::addToPositionChains" << std::endl;
  //////////

  // TODO: check the cost of the underlying memory allocations and cache the hash records buffer if needed
  std::vector<HashRecord>          hashRecords;
  std::vector<ExtendTableInterval> extendTableIntervals;
  unsigned                         fromHalfExtension       = 0;  // all seeds start as primary seeds
  const auto                       forwardData             = seed.getPrimaryData(false);
  const auto                       reverseData             = seed.getPrimaryData(true);
  const bool                       seedIsReverseComplement = (reverseData < forwardData);
  const auto                       primaryData = seedIsReverseComplement ? reverseData : forwardData;
  const auto                       hash        = getHashtable()->getPrimaryHasher()->getHash64(primaryData);
  getHashtable()->getHits(hash, false, hashRecords, extendTableIntervals);

  ////////////////
  // std::cerr << "Mapper::addToPositionChains: found " << hashRecords.size() << " hash records:";
  // for (const auto &record: hashRecords ) std::cerr << "  " << std::hex << std::setfill('0') << std::setw(8) << (record.getValue() >> 32) << ":" << std::setw(8) << (record.getValue() & 0xFFFFFFFF) << std::setfill(' ')<< std::dec<<"\t"<<typeString[record.getType()];
  // std::cerr << std::endl;
  ////////////////

  if (hashRecords.empty() and extendTableIntervals.empty()) {
    return;
  }
  // DEPRECATED - V7 only
  if (HashRecord::HIFREQ == hashRecords.front().getType()) {
    addRandomSamplesToPositionChains(
        seed, seedIsReverseComplement, fromHalfExtension, hashRecords, chainBuilder);
    return;
  }
  // at this point, it's either HIT or EXTEND records. First EXTEND as needed
  uint64_t extensionHash   = hash;
  bool     extensionFailed = false;

  // initialize variables for best interval tracking during extension
  uint32_t            nextStart           = 0;
  uint32_t            lastExtendTableSize = 0;  // defend against orphan EXTEND record
  BestIntervalTracker localBestIntvl(seed, 0, 0, 0);

  if (extendTableIntervals.size() != lastExtendTableSize) {  // new intervals seen
    if (extendTableIntervals.size() != 1) {
      boost::format message =
          boost::format("Expected size:1 extendTableIntervals, but got %i") % extendTableIntervals.size();
      BOOST_THROW_EXCEPTION(common::PostConditionException(message.str()));
    }
    nextStart = extendTableIntervals.front().getStart();
    localBestIntvl =
        BestIntervalTracker(seed, nextStart, extendTableIntervals.front().getLength(), fromHalfExtension);
    lastExtendTableSize = extendTableIntervals.size();
#ifdef TRACE_SEED_CHAINS
    std::cerr << "Init localBestIntvlPtr with\t" << localBestIntvl.getSeed()
              << "\tseedLength:" << localBestIntvl.getSeedLength()
              << "\thalfExtension:" << localBestIntvl.getHalfExtension()
              << "\tstart: " << localBestIntvl.getStart() << "\tlength: " << localBestIntvl.getLength()
              << std::endl;
#endif
  }

  while (!hashRecords.empty() && (HashRecord::EXTEND == hashRecords.front().getType())) {
    const uint64_t addressSegment = hash & addressSegmentMask_;
    const auto     extendRecord   = hashRecords.front();
    // DEPRECATED - V7 only
    addRandomSamplesToPositionChains(
        seed, seedIsReverseComplement, fromHalfExtension, hashRecords, chainBuilder);

    ////////////////
    // std::cerr << "Mapper::addToPositionChains: extented seed from " << fromHalfExtension << " to " << (fromHalfExtension + extendRecord.getExtensionLength() / 2) << ":";
    // for (const auto &record: hashRecords ) std::cerr << "  " << std::hex << std::setfill('0') <<
    // std::setw(8)
    // << (record.getValue() >> 32) << ":" << std::setw(8) << (record.getValue() & 0xFFFFFFFF) <<
    // std::setfill(' ')
    // << std::dec<<"\t"<<typeString[record.getType()];
    // std::cerr << std::endl;
    ////////////////

    hashRecords.clear();
    if (seed.isValid(fromHalfExtension + extendRecord.getExtensionLength() / 2)) {
      const auto extendedKey =
          getExtendedKey(seed, extensionHash, extendRecord, fromHalfExtension, seedIsReverseComplement);
      const auto extendedHash = addressSegment | getHashtable()->getSecondaryHasher()->getHash64(extendedKey);
      getHashtable()->getHits(extendedHash, true, hashRecords, extendTableIntervals);
      fromHalfExtension += extendRecord.getExtensionLength() / 2;
      extensionHash = extendedHash;

      // if extension failed, i.e. neither HIT nor INTERVAL
      if (hashRecords.empty() and lastExtendTableSize == extendTableIntervals.size()) extensionFailed = true;
      // local best interval tracking, process if extendTableIntervals is updated
      else if (lastExtendTableSize != extendTableIntervals.size()) {
        nextStart += extendTableIntervals.back().getStart();
        BestIntervalTracker currentIntvl(
            seed, nextStart, extendTableIntervals.back().getLength(), fromHalfExtension);
        if (!localBestIntvl.isValidInterval()) {
          if (extendTableIntervals.size() != 1) {
            boost::format message = boost::format("Expected size:1 extendTableIntervals, but got %i") %
                                    extendTableIntervals.size();
            BOOST_THROW_EXCEPTION(common::PostConditionException(message.str()));
          }
          localBestIntvl = currentIntvl;
#ifdef TRACE_SEED_CHAINS
          std::cerr << "Init(midway) localBestIntvlPtr with\t" << currentIntvl.getSeed()
                    << "\tseedLength:" << currentIntvl.getSeedLength()
                    << "\thalfExtension:" << currentIntvl.getHalfExtension()
                    << "\tstart: " << currentIntvl.getStart() << "\tlength: " << currentIntvl.getLength()
                    << std::endl;
#endif
        } else if (localBestIntvl.isWorseThan(currentIntvl)) {
#ifdef TRACE_SEED_CHAINS
          std::cerr << "Replace localBestIntvlPtr: Seed:" << localBestIntvl.getSeed()
                    << "\tseedLength:" << localBestIntvl.getSeedLength()
                    << "\thalfExtension:" << localBestIntvl.getHalfExtension()
                    << "\tstart: " << localBestIntvl.getStart() << "\tlength: " << localBestIntvl.getLength()
                    << "\twith\t" << currentIntvl.getSeed() << "\tseedLength:" << currentIntvl.getSeedLength()
                    << "\thalfExtension:" << currentIntvl.getHalfExtension()
                    << "\tstart: " << currentIntvl.getStart() << "\tlength: " << currentIntvl.getLength()
                    << std::endl;
#endif
          localBestIntvl = currentIntvl;
        } else {
#ifdef TRACE_SEED_CHAINS
          std::cerr << "Retain localBestIntvlPtr: Seed:" << localBestIntvl.getSeed()
                    << "\tseedLength:" << localBestIntvl.getSeedLength()
                    << "\thalfExtension:" << localBestIntvl.getHalfExtension()
                    << "\tstart: " << localBestIntvl.getStart() << "\tlength: " << localBestIntvl.getLength()
                    << "\tand skip\t" << currentIntvl.getSeed()
                    << "\tseedLength:" << currentIntvl.getSeedLength()
                    << "\thalfExtension:" << currentIntvl.getHalfExtension()
                    << "\tstart: " << currentIntvl.getStart() << "\tlength: " << currentIntvl.getLength()
                    << std::endl;
#endif
        }
        lastExtendTableSize = extendTableIntervals.size();
      }
      // else we see orphan EXTEND record or only HIT
    } else
      extensionFailed = true;
  }
  if (!extendTableIntervals.empty() &&
      extendTableIntervals.back().getLength() > getHashtable()->getMaxSeedFrequency())
    extensionFailed = true;
  // global best interval tracking
  if (localBestIntvl.isValidInterval()) globalBestIntvls.push_back(localBestIntvl);

#ifdef TRACE_SEED_CHAINS
  std::cerr << "Mapper::addToPositionChains: found " << extendTableIntervals.size()
            << " extendTableIntervals records:";
  for (const auto& record : extendTableIntervals)
    std::cerr << "  " << std::setfill('0') << std::setw(8) << (record.getStart()) << ":" << std::setw(8)
              << (record.getLength()) << std::setfill(' ') << "\thex:" << std::hex << std::setfill('0')
              << std::setw(8) << (record.getStart()) << ":" << std::setw(8) << (record.getLength())
              << std::dec;
  std::cerr << std::endl;
#endif

  // at this point there should be only HIT records (possibly 0) - add them to the seed chains
  if (!hashRecords.empty()) {
    bool isHitSeen = false;
    for (const auto& record : boost::adaptors::reverse(hashRecords)) {
      // special case encountered in alt-aware hashtables
      if (record.isDummyHit()) {
        continue;
      }
      if (!record.isHit()) {
        boost::format message =
            boost::format("Expected HIT hash record but got record of type: %i") % record.getType();
        BOOST_THROW_EXCEPTION(common::PostConditionException(message.str()));
      }
      isHitSeen                 = true;
      const bool isRandomSample = false;
      const bool orientation    = (seedIsReverseComplement ^ record.isReverseComplement());
      chainBuilder.addSeedPosition(
          SeedPosition(seed, record.getPosition(), fromHalfExtension), orientation, isRandomSample);
      // global best interval tracking
      if (longest_nonsample_seed_len < seed.getPrimaryLength() + 2 * fromHalfExtension) {
        longest_nonsample_seed_len = seed.getPrimaryLength() + 2 * fromHalfExtension;
      }
    }
    if (isHitSeen) return;
  } else if (!extendTableIntervals.empty()) {
    const uint32_t start = std::accumulate(
        extendTableIntervals.begin(),
        extendTableIntervals.end(),
        0UL,
        [](const uint32_t current, const ExtendTableInterval& interval) {
          return current + interval.getStart();
        });
    const uint32_t length = extendTableIntervals.back().getLength();

    if (extensionFailed) {
      num_extension_failure++;
      if (num_extension_failure > MAX_HIFREQ_HITS) return;

      std::vector<ExtendTableRecord> sampledExtendHashRecords;
      getRandomSamplesFromMatchInterval(seed, 1, start, length, sampledExtendHashRecords);
      if (!sampledExtendHashRecords.empty()) {
        const auto& record         = sampledExtendHashRecords.front();
        const bool  isRandomSample = true;
        const bool  orientation    = (seedIsReverseComplement ^ record.isReverseComplement());
        chainBuilder.addSeedPosition(
            SeedPosition(seed, record.getPosition(), fromHalfExtension), orientation, isRandomSample);
      }
    } else {
      assert(extendTableIntervals.back().getLength() <= getHashtable()->getMaxSeedFrequency());
      for (uint32_t i = start + length - 1; i != start - 1; --i) {
        const auto& record         = getHashtable()->getExtendTableRecord(i);
        const bool  isRandomSample = false;
        const bool  orientation    = (seedIsReverseComplement ^ record.isReverseComplement());
        chainBuilder.addSeedPosition(
            SeedPosition(seed, record.getPosition(), fromHalfExtension), orientation, isRandomSample);
      }
      // global best interval tracking
      if (longest_nonsample_seed_len < seed.getPrimaryLength() + 2 * fromHalfExtension) {
        longest_nonsample_seed_len = seed.getPrimaryLength() + 2 * fromHalfExtension;
      }
      // TODO: process extliftover hits if any
    }
  } else {  // could be orphan EXTEND
    //    boost::format message =
    //        boost::format("Expected non-empty extendTableIntervals, but got %i") % extendTableIntervals.size();
    //    BOOST_THROW_EXCEPTION(common::PostConditionException(message.str()));
  }
}

void Mapper::addExtraIntervalSamplesToPositionChains(
    const BestIntervalTracker& bestIntvl, ChainBuilder& chainBuilder) const
{
  Seed       seed                    = bestIntvl.getSeed();
  unsigned   fromHalfExtension       = bestIntvl.getHalfExtension();
  const auto forwardData             = seed.getPrimaryData(false);
  const auto reverseData             = seed.getPrimaryData(true);
  const bool seedIsReverseComplement = (reverseData < forwardData);

#ifdef TRACE_SEED_CHAINS
  std::cerr << "Global Best Interval sampling:"
            << "\tPrimarySeed:" << bestIntvl.getSeed() << "\tintervalStart:" << bestIntvl.getStart()
            << "\tintervalLength:" << bestIntvl.getLength()
            << "\thalfExtension:" << bestIntvl.getHalfExtension() << "\n";
#endif
  int nonExtraChainWaterMark = chainBuilder.size();
  if (bestIntvl.getLength() <= BestIntervalTracker::intvl_max_hits)  // fetch all
  {
    // the absolute positions must be calculated from the first interval
    const uint32_t start  = bestIntvl.getStart();
    const uint32_t length = bestIntvl.getLength();
    for (uint32_t i = start + length - 1; i != start - 1; --i) {
      const auto& record         = getHashtable()->getExtendTableRecord(i);
      const bool  isRandomSample = false;
      const bool  orientation    = (seedIsReverseComplement ^ record.isReverseComplement());
      chainBuilder.addSeedPosition(
          SeedPosition(seed, record.getPosition(), fromHalfExtension), orientation, isRandomSample);
    }
  } else  // random sample
  {
    std::vector<ExtendTableRecord> sampledExtendHashRecords;
    getRandomSamplesFromMatchInterval(
        seed,
        BestIntervalTracker::intvl_sample_hits,
        bestIntvl.getStart(),
        bestIntvl.getLength(),
        sampledExtendHashRecords);

    for (const auto& record : sampledExtendHashRecords) {
      const bool isRandomSample = true;
      const bool orientation    = (seedIsReverseComplement ^ record.isReverseComplement());
      chainBuilder.addSeedPosition(
          SeedPosition(seed, record.getPosition(), fromHalfExtension), orientation, isRandomSample);
    }
  }
  // set newly added seedChain as extra
  for (auto iter = chainBuilder.begin() + nonExtraChainWaterMark; iter != chainBuilder.end(); ++iter) {
    iter->setExtra(true);
  }
}

void Mapper::getRandomSamplesFromMatchInterval(
    const Seed&                     seed,
    const uint32_t                  sampleSize,
    const uint32_t                  intvl_start,
    const uint32_t                  intvl_len,
    std::vector<ExtendTableRecord>& hashRecords) const
{
  auto read     = seed.getRead();
  auto readName = read->getName();

  uint32_t sampledIndex = 0;
  uint32_t SEED         = 0;  // 32-bit SEED for random sampling
  uint32_t maxRounds    = 0;  // failsafe, maximum sampling rounds

  std::bitset<0x4000>          hitVector;
  std::unordered_set<uint32_t> alreadyFetchedIntvlPos;

  // Calculate SEED for random sampling
  // For 1 random sample after failed seed extension:
  if (1 == sampleSize) {
    maxRounds  = 1;
    uint32_t A = common::crc32c_hw(0, reinterpret_cast<const uint8_t*>(readName.data()), readName.size());
    uint32_t B = read->getPosition() << 31;
    B |= 0x40000000;
    B |= seed.getReadPosition() & 0x3FFFFFFF;
    SEED = A + B;
  }
  // For K random samples from a (single) extra interval:
  else if (1 < sampleSize) {
    maxRounds = 0x4000;
    SEED      = common::crc32c_hw(0, reinterpret_cast<const uint8_t*>(readName.data()), readName.size());
    SEED      = read->getPosition() ? (SEED ^(1<<31)) : SEED;
  } else {
    boost::format message = boost::format("Expected positive target sampleSize: %d") % sampleSize;
    BOOST_THROW_EXCEPTION(common::PostConditionException(message.str()));
  }

#ifdef TRACE_SEED_CHAINS
  std::cerr << "SEED:" << std::hex << SEED << "\t";
#endif

  uint32_t C = 0;
  uint32_t Z = 0;
  uint32_t K = 0;

  for (uint32_t X = SEED; K < sampleSize and X < SEED + maxRounds; X++) {
    C                  = common::crc32c_hw(C, reinterpret_cast<const uint8_t*>(&X), 4);
    sampledIndex       = std::floor(intvl_len * double(C) / std::pow(2, 32));
    const auto& record = getHashtable()->getExtendTableRecord(intvl_start + sampledIndex);

#ifdef TRACE_SEED_CHAINS
    std::cerr << "CRC(C):" << std::hex << C << std::dec << "\tX:" << X << "\tsampledIndex:" << sampledIndex
              << "\treadPosition:" << std::hex << record.getPosition() << std::dec << "\t";
#endif

    // filter record
    if (alreadyFetchedIntvlPos.find(record.getPosition()) != alreadyFetchedIntvlPos.end() ||
        ExtendTableRecord::LiftCode::ALT == record.getLiftCode() ||
        ExtendTableRecord::LiftCode::DIF_PRI == record.getLiftCode()) {
#ifdef TRACE_SEED_CHAINS
      std::cerr << "Failed" << std::endl;
#endif
      continue;
    }
    if (0x4000 < intvl_len) {
      Z            = SEED + sampledIndex;
      sampledIndex = common::crc32c_hw(0, reinterpret_cast<const uint8_t*>(&Z), 4) & 0x3FFF;
    }
    if (hitVector.test(sampledIndex)) {
#ifdef TRACE_SEED_CHAINS
      std::cerr << "Failed" << std::endl;
#endif

      continue;
    }

#ifdef TRACE_SEED_CHAINS
    std::cerr << "Survived" << std::endl;
#endif

    hitVector.set(sampledIndex);
    alreadyFetchedIntvlPos.insert(record.getPosition());

    hashRecords.push_back(record);
    K++;
  }
}

}  // namespace map
}  // namespace dragenos
