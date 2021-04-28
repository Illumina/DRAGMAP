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

#ifndef COMMON_BITS_HPP
#define COMMON_BITS_HPP

namespace dragenos {
namespace common {
namespace bits {

template <unsigned POSITION, typename T>
bool getFlag(const T value)
{
  static_assert(8 * sizeof(T) > POSITION, "POSITION greater than the number of bits in the type T");
  return (value >> POSITION) & 1;
}

template <unsigned BITS, typename T>
constexpr T getMask()
{
  static_assert(
      8 * sizeof(T) >= BITS, "number of BITS in the mask can't be more than number of bits in the type T");
  return (8 * sizeof(T) > BITS) ? (static_cast<T>(1) << BITS) - 1 : ~static_cast<T>(0);
}

template <unsigned START, unsigned BITS, typename T>
T getBits(const T value)
{
  static_assert(8 * sizeof(T) >= START, "START position can't be more than number of bits in the type T");
  static_assert(
      8 * sizeof(T) >= START + BITS,
      "START position plus number of selected BITS can't be more than number of bits in the type T");
  return (value >> START) & getMask<BITS, T>();
}

}  // namespace bits
}  // namespace common
}  // namespace dragenos

#endif  // #ifndef COMMON_BITS_HPP
