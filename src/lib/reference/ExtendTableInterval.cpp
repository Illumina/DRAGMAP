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

#include "reference/ExtendTableInterval.hpp"

#include <sstream>

#include "common/Exceptions.hpp"

namespace dragenos {
namespace reference {

template <typename I>
ExtendTableInterval::ExtendTableInterval(I begin, I end)
  : start_(extractStart(begin, end)),
    length_(extractLength(begin, end)),
    extraLiftoverMatchCount_(extractExtraLiftoverMatchCount(begin, end))
{
}

// explicit instantiation for vectors and pointers
template ExtendTableInterval::ExtendTableInterval(
    const HashRecord* begin, const HashRecord* end);
template ExtendTableInterval::ExtendTableInterval(
    std::vector<HashRecord>::const_iterator begin, std::vector<HashRecord>::const_iterator end);
template ExtendTableInterval::ExtendTableInterval(
    std::vector<HashRecord>::iterator begin, std::vector<HashRecord>::iterator end);

template <typename I>
bool ExtendTableInterval::isValid(I begin, I end)
{
  const auto size = end - begin;
  switch (size) {
  case 1: {
    const auto type = begin->getType();
    return ((HashRecord::INTERVAL_SL == type) && !begin->isMsb()) ||
           ((HashRecord::INTERVAL_SLE == type) && (begin->getExlift() > 0) && !begin->isMsb());
    break;
  }
  case 2: {
    const auto& r0    = *begin;
    const auto& r1    = *(begin + 1);
    const auto  type0 = r0.getType();
    const auto  type1 = r1.getType();
    return ((HashRecord::INTERVAL_SL == type0) && r0.isMsb() && (HashRecord::INTERVAL_S == type1)) ||
           ((HashRecord::INTERVAL_SLE == type0) && (r0.getExlift() > 0) && r0.isMsb() &&
            (HashRecord::INTERVAL_S == type1)) ||
           ((HashRecord::INTERVAL_SLE == type0) && (r0.getExlift() == 0) && (!r0.isMsb()) &&
            (HashRecord::INTERVAL_L == type1)) ||
           ((HashRecord::INTERVAL_S == type0) && (HashRecord::INTERVAL_L == type1));
    break;
  }
  case 3: {
    const auto& r0    = *begin;
    const auto& r1    = *(begin + 1);
    const auto& r2    = *(begin + 2);
    const auto  type0 = r0.getType();
    const auto  type1 = r1.getType();
    const auto  type2 = r2.getType();
    return (HashRecord::INTERVAL_SLE == type0) && (0 == r0.getExlift()) && r0.isMsb() &&
           (HashRecord::INTERVAL_S == type1) && (HashRecord::INTERVAL_L == type2);
    break;
  }
  default:
    return false;
    break;
  }
}

void ExtendTableInterval::unexpectedHashRecord(const HashRecord& hashRecord)
{
  std::ostringstream os;
  os << "unexpected hash records for an interval set:" << std::hex << std::setfill('0');
  const auto value = hashRecord.getValue();
  os << " " << std::setw(8) << (value >> 32) << ":" << (value & 0xffffffff);
  BOOST_THROW_EXCEPTION(common::InvalidParameterException(os.str()));
}

template <typename I>
uint32_t ExtendTableInterval::extractStart(const I begin, const I end)
{
  if (!isValid(begin, end)) {
    std::ostringstream os;
    os << "invalid combination of hash records for an interval set:" << std::hex << std::setfill('0');
    auto current = begin;
    while (current != end) {
      const auto value = current->getValue();
      os << " " << std::setw(8) << (value >> 32) << ":" << (value & 0xffffffff);
      ++current;
    }
    // TODO: verify that this test can safely be dropped
    // was thrown for: invalid combination of hash records for an interval set: f3faf5e5:f8101800
    //BOOST_THROW_EXCEPTION(common::InvalidParameterException(os.str()));
  }
  uint32_t start   = 0;
  I        current = begin;
  while (end != current) {
    // force the type resolution for convenience
    const HashRecord& hashRecord = *current;
    const auto        type       = hashRecord.getType();
    switch (type) {
    case HashRecord::INTERVAL_SL: {
      if (!hashRecord.isMsb())  // SL0
      {
        start += hashRecord.getBits<0, 15>();
      } else  // SL1
      {
        start += (hashRecord.getBits<0, 8>() << 24);
      }
      break;
    }
    case HashRecord::INTERVAL_SLE: {
      const unsigned shift = hashRecord.isMsb() ? 24 : 0;
      start += (hashRecord.getBits<0, 8>() << shift);
      break;
    }
    case HashRecord::INTERVAL_S: {
      start += hashRecord.getBits<0, 24>();
      if (hashRecord.hasCarry()) {
        start += (1UL << 24);
      }
      break;
    }
    case HashRecord::INTERVAL_L: {
      // no start, length only
      break;
    }
    default:
      unexpectedHashRecord(*begin);
      break;
    }
    ++current;
  }
  return start;
}

template <typename I>
uint32_t ExtendTableInterval::extractLength(const I begin, const I end)
{
  uint32_t length  = 0;
  auto     current = begin;
  while (end != current) {
    // force the type resolution for convenience
    const HashRecord& hashRecord = *current;
    const auto        type       = hashRecord.getType();
    switch (type) {
    case HashRecord::INTERVAL_SL: {
      if (!hashRecord.isMsb())  // SL0
      {
        length += hashRecord.getBits<15, 9>();
      } else  // SL1
      {
        length += hashRecord.getBits<8, 16>();
      }
      break;
    }
    case HashRecord::INTERVAL_SLE: {
      if (0 < hashRecord.getExlift()) {
        length += hashRecord.getBits<8, 8>();
      } else  // Exlift == 0
      {
        length += (hashRecord.getBits<8, 8>() << 24);
      }
      break;
    }
    case HashRecord::INTERVAL_S: {
      // no length, length only
      break;
    }
    case HashRecord::INTERVAL_L: {
      length += hashRecord.getBits<0, 24>();
      break;
    }
    default:
      unexpectedHashRecord(*begin);
      break;
    }
    ++current;
  }
  return length;
}

template <typename I>
uint8_t ExtendTableInterval::extractExtraLiftoverMatchCount(const I begin, I /* end */)
{
  if (HashRecord::INTERVAL_SLE == begin->getType()) {
    const HashRecord& hashRecord = *begin;
    return hashRecord.getBits<16, 8>();
  }
  return 0;
}

}  // namespace reference
}  // namespace dragenos
