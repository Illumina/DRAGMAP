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

#ifndef MAP_SEED_POSITION_HPP
#define MAP_SEED_POSITION_HPP

#include <ostream>

#include "sequences/Seed.hpp"

namespace dragenos {
namespace map {

class SeedPosition {
public:
  typedef uint32_t ReferencePosition;
  SeedPosition(const sequences::Seed& seed, ReferencePosition position, const unsigned halfExtension)
    : seed_(seed), position_(position), halfExtension_(halfExtension)
  {
  }
  ReferencePosition      getReferencePosition() const { return position_; }
  const sequences::Seed& getSeed() const { return seed_; }
  unsigned         getHalfExtension() const { return halfExtension_; }
  /// Projected reference position of the leftmost base of the read
  ReferencePosition getFirstProjection(bool reverseComplement) const;
  /// Projected reference position of the rightmost base of the read
  ReferencePosition getLastProjection(bool reverseComplement) const;

  unsigned getFirstBaseReferencePosition() const { return position_ - halfExtension_; }

  unsigned getLastBaseReferencePosition() const
  {
    return position_ + seed_.getPrimaryLength() + halfExtension_ - 1;
  }

  /// Position of the last base of the extended seed on the read
  unsigned getFirstBaseReadPosition() const { return seed_.getFirstBaseReadPosition(halfExtension_); }

  /// Position of the last base of the extended seed on the read
  unsigned getLastBaseReadPosition() const { return seed_.getLastBaseReadPosition(halfExtension_); }

  friend std::ostream& operator<<(std::ostream& os, const SeedPosition& sp)
  {
    return os << "SeedPosition(" << sp.seed_ << "," << std::hex << sp.position_ << std::dec << ","
              << sp.halfExtension_ << ")";
  }

  bool operator==(const SeedPosition& that) const
  {
    return seed_ == that.seed_ && position_ == that.position_ && halfExtension_ == that.halfExtension_;
  }

private:
  sequences::Seed   seed_;
  ReferencePosition position_;
  unsigned          halfExtension_;
};

}  // namespace map
}  // namespace dragenos

#endif  // #ifndef MAP_SEED_POSITION_HPP
