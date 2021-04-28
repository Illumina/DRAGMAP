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

#include "map/SeedPosition.hpp"

namespace dragenos {
namespace map {

SeedPosition::ReferencePosition SeedPosition::getFirstProjection(bool reverseComplement) const
{
  if (reverseComplement) {
    const auto reverseReadPosition = (seed_.getRead()->getLength() > seed_.getReadPosition())
                                         ? seed_.getRead()->getLength() - seed_.getReadPosition() - 1
                                         : 0;
    return position_ > reverseReadPosition ? position_ - reverseReadPosition : 0;
  } else {
    return (position_ > seed_.getReadPosition()) ? position_ - seed_.getReadPosition() : 0;
  }
}

SeedPosition::ReferencePosition SeedPosition::getLastProjection(bool reverseComplement) const
{
  return SeedPosition::getFirstProjection(reverseComplement) + seed_.getRead()->getLength() - 1;
}

}  // namespace map
}  // namespace dragenos
