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

#ifndef MAP_SEED_CHAIN_HPP
#define MAP_SEED_CHAIN_HPP

#include <boost/io/ios_state.hpp>
#include <iomanip>
#include <limits>
#include <map>

#include "map/SeedPosition.hpp"

namespace dragenos {
namespace map {

/**
 ** \brief A class that aggregates seed mappings into coherent scaffolds for read alignment
 **
 ** Seeds are stored in the chain by strictly increasing read position of the last base of the seed
 ** (i.e. seed read position plus seed length plus half-extension).
 ** Each seed has an "age" and a "diagonal". The age is the distance along the read between the last
 ** base of the seed and the last base of the current (or planned) head of the chain. The diagonal is
 ** the projected reference position of the first base of the read assuming a gapless alignment from
 ** the seed reference position.
 ** The ages are classified in three categories: active, old and ancient. The diagonal is used to
 ** calculate the diameter and radius of the chain:
 ** - diameter: distance between the leftmost and rightmost diagonals across relevant seeds
 ** - radius: increase in diameter between consecutive seeds
 ** The diameter and radius are used to control which seeds can be added to the chain. The diameter
 ** test discards the ancient seeds and the radius test only uses the active seeds. This enables
 ** some drift in the wavefront used later in the smith waterman alignment.
 **
 **/
class SeedChain {
public:
  typedef sequences::Seed                      Seed;
  typedef map::SeedPosition::ReferencePosition ReferencePosition;
  static constexpr unsigned                    LARGE_QUANTIZER = 16;
  static constexpr unsigned                    SMALL_QUANTIZER = 4;
  static constexpr unsigned                    MAX_DIAMETER    = 8;  // 4 base, SMALL_QUANTIZED
  static constexpr unsigned                    MAX_RADIUS      = 5;
  static constexpr unsigned                    OLD             = 9;  // 16 base, LARGE_QUANTIZED
  static constexpr unsigned                    ANCIENT         = 31;

  static constexpr unsigned DIAG_HALF_BIN_SIZE = 256;
  /// default constructor: not reverse complement, assume only random samples
  SeedChain() : reverseComplement_(false), randomSamplesOnly_(true) {}

  void clear()
  {
    seedPositions_.clear();
    diagonalTable_.clear();
    reverseComplement_ = false;
    randomSamplesOnly_ = true;
    initialDiagonal_   = 0;
    perfectAlignment_  = true;
    filtered_          = false;
    needRescue_        = true;
    extra_             = false;
    coverage_          = 0;
    firstReadBase_     = std::numeric_limits<uint32_t>::max();
    lastReadBase_      = std::numeric_limits<uint32_t>::min();
    firstChainRefBase_ = std::numeric_limits<uint32_t>::max();
    lastChainRefBase_  = std::numeric_limits<uint32_t>::min();

    firstSeedReferencePosition_ = std::numeric_limits<uint32_t>::max();
    lastSeedReferencePosition_  = std::numeric_limits<uint32_t>::min();
  }
  void setReverseComplement(bool reverseComplement) { reverseComplement_ = reverseComplement; }
  bool isReverseComplement() const { return reverseComplement_; }
  bool hasOnlyRandomSamples() const { return randomSamplesOnly_; }
  bool setRandomSamplesOnly(bool randomSamplesOnly) { return randomSamplesOnly_ = randomSamplesOnly; }
  /**
   ** \brief check if the chain would accept a seed mapping to the reference with that position and
   *orientation
   **
   **/
  bool accepts(const SeedPosition& seedPosition, const bool rcFlag) const;
  /**
   ** \brief  add a new element to the chain. It is the responsibility of the client code to check acceptance
   **
   ** \param rsFlag flag indicating if the mapping is a random sample (from high frequency seed or seed that
   ** needs an extension
   **/
  //void add(const Seed &seed, const uint32_t referencePosition, bool rsFlag);
  void addSeedPosition(const SeedPosition& seedPosition, bool rsFlag)
  {
    if (empty()) {
      initialDiagonal_ = getDiagonal(seedPosition);
    } else {
      if (getDiagonal(seedPosition) != initialDiagonal_) {
        perfectAlignment_ = false;
      }
    }
    updateSeedChainInfo(seedPosition);
    seedPositions_.emplace_back(seedPosition);
    randomSamplesOnly_ &= rsFlag;
  }
  /** \brief calculate diagonal from seed read position, reference position, seed length and orientation
   **
   ** The diagonal is defined as the reference position of the leftmost base in the read (the base with the
   ** smallest reference position). This would be:
   ** - Forward chain: reference position of the seed minus the seed read position
   ** - Reverse chain: reference position of the seed minus (read length - seed length - seed read position -
   *1)
   **/
  uint32_t getDiagonal(const Seed& seed, const uint32_t referencePosition) const;
  uint32_t getDiagonal(const SeedPosition& seedPosition) const
  {
    return getDiagonal(seedPosition.getSeed(), seedPosition.getReferencePosition());
  }
  /// checks that the new diagonal does not conflict with rc flag
  bool passesInversionTest(const SeedPosition& seedPosition) const;
  /// checks that the new diagonal would fit within the diameter of all non-ancient diagonals
  bool passesDiameterTest(const unsigned newStartOffset, const uint32_t diagonal) const;
  /// checks that the new diagonal would fit within the radius od all active diagonals
  bool passesRadiusTest(const uint32_t diagonal) const;
  /// checks if the new seed would make ancient all existing elements of the chain
  bool terminates(const Seed& seed) const;

  ReferencePosition firstSeedReferencePosition() const
  {
    return reverseComplement_ ? lastSeedReferencePosition_ : firstSeedReferencePosition_;
  }
  ReferencePosition lastSeedReferencePosition() const
  {
    return reverseComplement_ ? firstSeedReferencePosition_ : lastSeedReferencePosition_;
  }

  /// smallest diagonal across all seeds in the chain
  ReferencePosition firstReferencePosition() const { return firstChainRefBase_; }
  /// rightmost base position projected across all the seeds in the chain
  ReferencePosition lastReferencePosition() const { return lastChainRefBase_; }
  /// offset of the closest base to the beginning of the read that is covered by a seed (primary or extension
  /// wings)
  unsigned getLength() const
  {
    return lastReferencePosition() > firstReferencePosition()
               ? lastReferencePosition() - firstReferencePosition()
               : 0;
  }

  /// offset of the closest base to the beginning of the read that is covered by a seed (primary or extension
  /// wings)
  unsigned firstReadBase() const { return firstReadBase_; }
  /// offset of the closest base to the end of the read that is covered by a seed (primary or extension wings)
  unsigned lastReadBase() const { return lastReadBase_; }
  /// length of the seed chain spanning over the read
  unsigned getReadSpanLength() const
  {
    if (lastReadBase() < firstReadBase()) {
      std::cerr << "WARNING: lastReadBase()=" << lastReadBase() << ", firstReadBase()=" << firstReadBase()
                << ", SeedChain.size()=" << size() << std::endl;
      return 0;
    }
    return lastReadBase() - firstReadBase();
  }
  /// coverage of the seed chain spanning over the read
  unsigned getReadCovLength() const { return coverage_; }

  /// update first and last read base
  void updateRefBase(const SeedPosition& a);
  /// update first and last read base
  void updateReadBase(const SeedPosition& a);

  void updateDiagonalTable(const SeedPosition& seedPosition);
  /// call all updates functions
  void updateSeedChainInfo(const SeedPosition& seedPosition)
  {
    updateReadBase(seedPosition);
    updateRefBase(seedPosition);
    updateDiagonalTable(seedPosition);
  }

  friend std::ostream& operator<<(std::ostream& os, const SeedChain& chain)
  {
    boost::io::ios_flags_saver ifs(os);
    return os << "CHAIN:"
              << " RC=" << chain.isReverseComplement() << ", qryPos1=" << chain.firstReadBase()
              << ", qryPos2=" << chain.lastReadBase() << ", refPos1=0x" << std::hex << std::uppercase
              << std::setw(9) << std::setfill('0')
              << chain.firstSeedReferencePosition()  //<< ", " << chain.firstChainRefBase_
              << ", refPos2=0x" << std::setw(9)
              << chain.lastSeedReferencePosition()  //<< ", " << chain.lastChainRefBase_
              << ", perf=" << chain.isPerfect() << ", diag=" << chain.initialDiagonal_
              << ", filt=" << chain.isFiltered() << ", resc=" << chain.needRescue() << ", cov=" << std::dec
              << chain.getReadCovLength();
  }

  bool isPerfect() const { return perfectAlignment_; }
  void setPerfect(bool perfect) { perfectAlignment_ = perfect; }

  bool isFiltered() const { return filtered_; }
  void setFiltered(bool status) { filtered_ = status; }

  bool needRescue() const { return needRescue_; }
  void setNeedRescue(bool status) { needRescue_ = status; }

  bool isExtra() const { return extra_; }
  void setExtra(bool status) { extra_ = status; }

private:
  /// true if the seed chain is Reverse-Complement
  bool reverseComplement_;
  /// false if there is at least one non-random-sample in the chain
  bool                         randomSamplesOnly_;
  std::vector<SeedPosition>    seedPositions_;
  std::map<uint32_t, uint32_t> diagonalTable_;  // k:diagonal v:lastSeedOffset

  uint32_t initialDiagonal_  = 0;
  bool     perfectAlignment_ = true;
  bool     filtered_         = false;
  bool     needRescue_       = true;
  bool     extra_            = false;

  uint32_t coverage_ = 0;

  /// tmp variables
  uint32_t firstReadBase_     = std::numeric_limits<uint32_t>::max();
  uint32_t lastReadBase_      = std::numeric_limits<uint32_t>::min();
  uint32_t firstChainRefBase_ = std::numeric_limits<uint32_t>::max();
  uint32_t lastChainRefBase_  = std::numeric_limits<uint32_t>::min();

  uint32_t firstSeedReferencePosition_ = std::numeric_limits<uint32_t>::max();
  uint32_t lastSeedReferencePosition_  = std::numeric_limits<uint32_t>::min();

public:
  auto size() const -> decltype(seedPositions_.size()) { return seedPositions_.size(); }
  auto empty() const -> decltype(seedPositions_.empty()) { return seedPositions_.empty(); }
  auto front() const -> decltype(seedPositions_.front()) { return seedPositions_.front(); }
  auto back() const -> decltype(seedPositions_.back()) { return seedPositions_.back(); }
  auto begin() const -> decltype(seedPositions_.begin()) { return seedPositions_.begin(); }
  auto end() const -> decltype(seedPositions_.end()) { return seedPositions_.end(); }
  auto rbegin() const -> decltype(seedPositions_.rbegin()) { return seedPositions_.rbegin(); }
  auto rend() const -> decltype(seedPositions_.rend()) { return seedPositions_.rend(); }
};

}  // namespace map
}  // namespace dragenos

#endif  // #ifndef MAP_SEED_CHAIN_HPP
