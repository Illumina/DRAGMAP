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

#ifndef MAP_BEST_INTERVAL_TRACKER_HPP
#define MAP_BEST_INTERVAL_TRACKER_HPP

#include <ostream>

#include "map/ChainBuilder.hpp"
#include "map/SeedChain.hpp"
#include "sequences/Seed.hpp"

namespace dragenos {
namespace map {

/// tracking best interval
class BestIntervalTracker {
public:
  /// Seed Match Interval Processing parameters and defaults
  static const uint32_t intvl_target_hits        = 32;
  static const uint32_t intvl_min_chains         = 8;
  static const uint32_t intvl_seed_longer_length = 8;
  static const uint32_t intvl_seed_length        = 60;
  static const uint32_t intvl_max_hits           = 16;
  static const uint32_t intvl_sample_hits        = 16;

  BestIntervalTracker(
      const sequences::Seed& primary_seed,
      uint32_t               intvl_start,
      uint32_t               intvl_len,
      uint32_t               fromHalfExtension);
  uint32_t        getStart() const;
  uint32_t        getLength() const;
  uint32_t        getSeedLength() const;
  uint32_t        getHalfExtension() const;
  sequences::Seed getSeed() const;
  bool            isWorseThan(const BestIntervalTracker& a) const;
  bool            isValidExtra(uint32_t num_non_sample_seed_chains, uint32_t longest_nonsample_seed_len);
  bool            isValidInterval() const { return intvl_len_ != 0; }

private:
  sequences::Seed primary_seed_;
  uint32_t        intvl_start_       = 0;  // keep the absolute pos of the extendTableRecord
  uint32_t        intvl_len_         = 0;
  uint32_t        fromHalfExtension_ = 0;
  uint32_t        capped_intvl_len_  = 0;
};

}  // namespace map
}  // namespace dragenos

#endif  // #ifndef MAP_BEST_INTERVAL_TRACKER_HPP
