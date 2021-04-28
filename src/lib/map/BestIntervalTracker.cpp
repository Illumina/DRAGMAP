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

#include "map/BestIntervalTracker.hpp"

namespace dragenos {
namespace map {

BestIntervalTracker::BestIntervalTracker(
    const sequences::Seed& primary_seed, uint32_t intvl_start, uint32_t intvl_len, uint32_t fromHalfExtension)
  : primary_seed_(primary_seed),
    intvl_start_(intvl_start),
    intvl_len_(intvl_len),
    fromHalfExtension_(fromHalfExtension)
{
  capped_intvl_len_ = intvl_len_ < intvl_target_hits ? intvl_len_ : intvl_target_hits;
}

uint32_t BestIntervalTracker::getStart() const
{
  return intvl_start_;
}

uint32_t BestIntervalTracker::getLength() const
{
  return intvl_len_;
}

uint32_t BestIntervalTracker::getSeedLength() const
{
  return primary_seed_.getPrimaryLength() + 2 * fromHalfExtension_;
}

uint32_t BestIntervalTracker::getHalfExtension() const
{
  return fromHalfExtension_;
}

sequences::Seed BestIntervalTracker::getSeed() const
{
  return primary_seed_;
}

bool BestIntervalTracker::isWorseThan(const BestIntervalTracker& a) const
{
  if (fromHalfExtension_ == 0 and intvl_len_ < intvl_target_hits and
      getSeed().getReadPosition() == a.getSeed().getReadPosition())
    return false;
  else
    return (capped_intvl_len_ < a.capped_intvl_len_) or
           (capped_intvl_len_ == a.capped_intvl_len_ and getSeedLength() < a.getSeedLength());
}

bool BestIntervalTracker::isValidExtra(
    uint32_t num_non_sample_seed_chains, uint32_t longest_nonsample_seed_len)
{
  return (
      num_non_sample_seed_chains < intvl_min_chains or getSeedLength() >= intvl_seed_length or
      getSeedLength() >= longest_nonsample_seed_len + intvl_seed_longer_length);
}

}  // namespace map
}  // namespace dragenos
