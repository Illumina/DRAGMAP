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

#include "map/SeedChain.hpp"

#include <boost/assert.hpp>
#include <cstdlib>
#include <limits>

namespace dragenos {
namespace map {

bool SeedChain::accepts(const SeedPosition& seedPosition, const bool reverseComplement) const
{
  if (empty()) {
    return true;
  }

  if ((reverseComplement != reverseComplement_) || terminates(seedPosition.getSeed())) {
    return false;
  }
  const auto diagonal = getDiagonal(seedPosition.getSeed(), seedPosition.getReferencePosition());

  // 512+ away tests
  auto diag_ph = diagonal + DIAG_HALF_BIN_SIZE;
  auto diag_mh = diagonal - DIAG_HALF_BIN_SIZE;
  if (seedPosition.getLastBaseReadPosition() < lastReadBase_) {
    diag_ph = diagonal + MAX_RADIUS * SMALL_QUANTIZER;
    diag_mh = diagonal - MAX_RADIUS * SMALL_QUANTIZER;
  }
  if (initialDiagonal_ <= diag_mh or initialDiagonal_ >= diag_ph) return false;

  // Reject if the new seed is further than the "max age" before the original first seed
  if (seedPosition.getFirstBaseReadPosition() / LARGE_QUANTIZER + ANCIENT < firstReadBase_ / LARGE_QUANTIZER)
    return false;

  if (not passesInversionTest(seedPosition)) return false;
#ifdef TRACE_SEED_CHAINS
  std::cerr << "  passesInversionTest " << std::endl;
#endif
  // we have a hit in the right orientation at a suitable position along the read
  return passesDiameterTest(seedPosition.getSeed().getReadPosition(), diagonal) && passesRadiusTest(diagonal);
}

uint32_t SeedChain::getDiagonal(const Seed& seed, const uint32_t referencePosition) const
{
  return referencePosition + (reverseComplement_ ? seed.getReadPosition() : -seed.getReadPosition());
}

bool SeedChain::passesInversionTest(const SeedPosition& seedPosition) const
{
  // original values
  uint64_t firstBaseRefPos = firstSeedReferencePosition_;
  uint64_t lastBaseRefPos  = lastSeedReferencePosition_;

  // update with new if the seed extends the chain
  if (firstReadBase() > seedPosition.getFirstBaseReadPosition()) {
    firstBaseRefPos = reverseComplement_ ? seedPosition.getLastBaseReferencePosition()
                                         : seedPosition.getFirstBaseReferencePosition();
  }

  if (lastReadBase() < seedPosition.getLastBaseReadPosition()) {
    lastBaseRefPos = reverseComplement_ ? seedPosition.getFirstBaseReferencePosition()
                                        : seedPosition.getLastBaseReferencePosition();
  }

  // check that we have not inverted the seed chain
  if (reverseComplement_) {
    if (firstBaseRefPos < lastBaseRefPos) {
      //      std::cerr << "INVERSION: " << seedPosition <<
      //        "first_ref_pos " << std::hex << firstBaseRefPos <<
      //        " < " << lastBaseRefPos << " last_ref_pos" << *this << std::endl;
      return false;
    }
  } else {
    if (firstBaseRefPos > lastBaseRefPos) {
      //      std::cerr << "INVERSION: " << seedPosition <<
      //        "first_ref_pos " << std::hex << firstBaseRefPos <<
      //        " > " << lastBaseRefPos << " last_ref_pos" << *this << std::endl;
      return false;
    }
  }
  return true;
}

bool SeedChain::passesDiameterTest(const unsigned newStartOffset, const uint32_t diagonal) const
{
  BOOST_ASSERT(!diagonalTable_.empty());
  for (const auto& kv : diagonalTable_) {
    if (kv.second / LARGE_QUANTIZER + OLD <= newStartOffset / LARGE_QUANTIZER) {
      continue;
    }
    if (std::abs(static_cast<long>(diagonal) / SMALL_QUANTIZER - kv.first / SMALL_QUANTIZER) > MAX_DIAMETER)
      return false;
  }
#ifdef TRACE_SEED_CHAINS
  std::cerr << "  passesDiameterTest " << std::endl;
#endif
  return true;
}

bool SeedChain::passesRadiusTest(const uint32_t diagonal) const
{
  BOOST_ASSERT(!diagonalTable_.empty());
  uint32_t leftmost  = diagonalTable_.begin()->first;
  uint32_t rightmost = diagonalTable_.rbegin()->first;

  return (diagonal / SMALL_QUANTIZER + MAX_RADIUS >= leftmost / SMALL_QUANTIZER) &&
         (rightmost / SMALL_QUANTIZER + MAX_RADIUS >= diagonal / SMALL_QUANTIZER);
}

bool SeedChain::terminates(const Seed& seed) const
{
  return (!seedPositions_.empty() and diagonalTable_.empty());
}

void SeedChain::updateReadBase(const SeedPosition& a)
{
  // update coverage
  int seedLen = a.getSeed().getPrimaryLength() + a.getHalfExtension() * 2;
  if (firstReadBase_ != std::numeric_limits<uint32_t>::max() and
      lastReadBase_ != std::numeric_limits<uint32_t>::min()) {
    if (a.getLastBaseReadPosition() < firstReadBase_)
      coverage_ += seedLen;
    else if (a.getFirstBaseReadPosition() < firstReadBase_ and firstReadBase_ <= a.getLastBaseReadPosition())
      coverage_ += firstReadBase_ - a.getFirstBaseReadPosition();

    if (a.getFirstBaseReadPosition() <= lastReadBase_ and lastReadBase_ < a.getLastBaseReadPosition())
      coverage_ += a.getLastBaseReadPosition() - lastReadBase_;
    else if (lastReadBase_ < a.getFirstBaseReadPosition())
      coverage_ += seedLen;
  } else
    coverage_ += seedLen;
  // update endpoints
  if (a.getFirstBaseReadPosition() < firstReadBase_) {
    firstReadBase_ = a.getFirstBaseReadPosition();
    firstSeedReferencePosition_ =
        reverseComplement_ ? a.getLastBaseReferencePosition() : a.getFirstBaseReferencePosition();
  }
  if (a.getLastBaseReadPosition() > lastReadBase_) {
    lastReadBase_ = a.getLastBaseReadPosition();
    lastSeedReferencePosition_ =
        reverseComplement_ ? a.getFirstBaseReferencePosition() : a.getLastBaseReferencePosition();
  }
}

void SeedChain::updateDiagonalTable(const SeedPosition& seedPosition)
{
  auto lastSeedOffset = seedPosition.getSeed().getReadPosition();
  for (auto it = diagonalTable_.cbegin(); it != diagonalTable_.cend();) {
    if (it->second / LARGE_QUANTIZER + ANCIENT < lastSeedOffset / LARGE_QUANTIZER) {
      diagonalTable_.erase(it++);
    } else {
      ++it;
    }
  }
  diagonalTable_[getDiagonal(seedPosition)] = lastSeedOffset;
}

void SeedChain::updateRefBase(const SeedPosition& a)
{
  if (empty() or a.getSeed().getPrimaryData(false) != back().getSeed().getPrimaryData(false))
    firstChainRefBase_ = std::min(
        firstChainRefBase_,
        isReverseComplement() ? a.getFirstProjection(true) + a.getSeed().getPrimaryLength() - 1
                              : a.getFirstProjection(false));
  lastChainRefBase_ = std::max(
      lastChainRefBase_,
      isReverseComplement() ? a.getLastProjection(true) + a.getSeed().getPrimaryLength() - 1
                            : a.getLastProjection(false));
}

}  // namespace map
}  // namespace dragenos
