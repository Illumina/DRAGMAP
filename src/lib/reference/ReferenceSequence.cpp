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

#include "reference/ReferenceSequence.hpp"

#include <boost/format.hpp>

namespace dragenos {
namespace reference {

unsigned char ReferenceSequence::getBase(size_t position) const
{
#if 0
  size_t address = position;
  auto extractBase = [this](size_t address) -> unsigned char
  {
    const unsigned char result = data_[address / 2];
    return ((address % 2) ? (result >> 4) : result) & 0xF;
  };
  for (const auto &trimmedRegion: trimmedRegions_)
  {
    if (position >= trimmedRegion[1])
    {
      // adjust address to account for trimmed bases
      address -= (trimmedRegion[1] - trimmedRegion[0]);
    }
    else if (position >= trimmedRegion[0])
    {
      return 15; // 15 == encoded N
    }
    else
    {
      // return base in sequence before trimmed region
      return extractBase(address);
    }
  }
  // if we get here we are past the last trim region - check if still be in the sequence
  if (address / 2 >= size_)
  {
    BOOST_THROW_EXCEPTION(common::InvalidParameterException("position greater than reference size"));
  }
  return extractBase(address);
#endif
  if (position / 2 >= size_) {
    boost::format message =
        boost::format("position greater than reference size: %i > 2 * %i") % position % size_;
    BOOST_THROW_EXCEPTION(common::InvalidParameterException(message.str()));
  }
  const unsigned char twoBases = data_[position / 2];
  const bool          msb      = (position % 2);  // use the 4 MSB for odd positions
  return msb ? (twoBases >> 4) : (twoBases & 0xF);
}

unsigned char ReferenceSequence::getRcBase(size_t position) const
{
  unsigned char                              b = getBase(position);
  const static std::array<unsigned char, 16> translate{0b0000,   // 0b0000
                                                       0b1000,   // 0b0001
                                                       0b0100,   // 0b0010
                                                       0b1100,   // 0b0011
                                                       0b0010,   // 0b0100
                                                       0b1010,   // 0b0101
                                                       0b0110,   // 0b0110
                                                       0b1110,   // 0b0111
                                                       0b0001,   // 0b1000
                                                       0b1001,   // 0b1001
                                                       0b0101,   // 0b1010
                                                       0b1101,   // 0b1011
                                                       0b0011,   // 0b1100
                                                       0b1011,   // 0b1101
                                                       0b0111,   // 0b1110
                                                       0b1111};  // 0b1111

  return translate[b];
}

char ReferenceSequence::decodeBase(unsigned char base)
{
  const static std::string translate{"PACMGRSVTWYHKDBN"};
  return translate[base & 0xF];
}

unsigned char ReferenceSequence::translateTo2bpb(unsigned char base4bpb)
{
  const static std::array<unsigned char, 16> translate{0, 0, 1, 0, 2, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0};
  return translate[base4bpb & 0xF];
}

unsigned char ReferenceSequence::translateToR2bpb(unsigned char base4bpb)
{
  const static std::array<unsigned char, 16> translate{0, 3, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  return translate[base4bpb & 0xF];
}

}  // namespace reference
}  // namespace dragenos
