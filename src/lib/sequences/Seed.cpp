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

#include "sequences/Seed.hpp"

#include <cassert>
#include <sstream>

#include "common/Exceptions.hpp"

namespace dragenos {
namespace sequences {

bool Seed::isValid(const Read& read, unsigned readPosition, unsigned primaryLength)
{
  static constexpr Read::Base N = 0xf;
  if (readPosition + primaryLength > read.getLength()) {
    return false;
  }
  for (unsigned i = 0; primaryLength > i; ++i) {
    const Read::Base base = read.getBase4bpb(readPosition + i);
    if ((0 == base) || (N == base)) {
      return false;
    }
  }
  return true;
}

Seed::Seed(const Read* read, unsigned readPosition, unsigned primaryLength)
  : read_(read), readPosition_(readPosition), primaryLength_(primaryLength)
{
  // 2 bits per base - 4 bases per byte
  if (primaryLength_ > sizeof(Data) * 4) {
    BOOST_THROW_EXCEPTION(common::InvalidParameterException("Seed primary length is mimited to 32 bases"));
  }
}

Seed::Data Seed::generateReverseComplement(Seed::Data data, const unsigned baseCount)
{
  Data complement = ~data;
  Data reversed   = 0;
  for (unsigned i = 0; baseCount > i; ++i) {
    reversed <<= 2;
    reversed |= (complement & 3);
    complement >>= 2;
  }
  return reversed;
}

Seed::Data Seed::getPrimaryData(const bool reverseComplement) const
{
  if (read_->getLength() < readPosition_ + primaryLength_) {
    BOOST_THROW_EXCEPTION(common::PreConditionException("Requesting primary data for an invalid seed"));
  }
  Data data = 0;
  auto dest = reinterpret_cast<uint8_t*>(&data);
  for (unsigned i = 0; i < primaryLength_; ++i) {
    const uint8_t base = read_->getBase2bpb(i + readPosition_) & 0x03;
    dest[i / 4] |= base << (2 * (i % 4));
  }
  return reverseComplement ? generateReverseComplement(data, getPrimaryLength()) : data;
}

Seed::Data Seed::getExtendedData(
    const unsigned fromHalfExtension, const unsigned toHalfExtension, const bool reverseComplement) const
{
  if (!isValid(toHalfExtension)) {
    BOOST_THROW_EXCEPTION(
        common::PreConditionException("Requesting extended data for an invalid seed extension"));
  }
  if ((fromHalfExtension > toHalfExtension) || (toHalfExtension - fromHalfExtension > 16)) {
    BOOST_THROW_EXCEPTION(common::InvalidParameterException("Requesting extended data with invalid range"));
  }
  const unsigned wingLength = toHalfExtension - fromHalfExtension;
  Data           data       = 0;
  auto           dest       = reinterpret_cast<uint8_t*>(&data);
  for (unsigned i = 0; wingLength > i; ++i) {
    const uint8_t base = read_->getBase2bpb(readPosition_ - toHalfExtension + i) & 0x03;
    dest[i / 4] |= base << (2 * (i % 4));
  }
  for (unsigned i = 0; wingLength > i; ++i) {
    const uint8_t base = read_->getBase2bpb(readPosition_ + primaryLength_ + fromHalfExtension + i) & 0x03;
    dest[(i + wingLength) / 4] |= base << (2 * ((i + wingLength) % 4));
  }
  return reverseComplement ? generateReverseComplement(data, 2 * wingLength) : data;
}

std::vector<size_t> Seed::getSeedOffsets(
    const size_t   readLength,
    const unsigned seedLength,
    const uint32_t period,
    const uint32_t pattern,
    const uint8_t  force)
{
  std::vector<size_t> seedOffsets;
  size_t              offset = 0;
  while (offset + seedLength <= readLength) {
    const bool forced         = (offset + seedLength + force > readLength);
    const bool matchesPattern = ((pattern >> (offset % period)) & 1);
    if (matchesPattern || forced) {
      seedOffsets.push_back(offset);
    }
    ++offset;
  }
  return seedOffsets;
}

}  // namespace sequences
}  // namespace dragenos
