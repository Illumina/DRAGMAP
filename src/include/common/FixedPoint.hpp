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

#ifndef COMMON_FIXED_POINT_HPP
#define COMMON_FIXED_POINT_HPP

#include <inttypes.h>
#include <boost/assert.hpp>
#include <cmath>

//#include "dragen_exception.hpp"

namespace dragenos {
namespace common {

class FixedPointNumber {
private:
  typedef int32_t value_type;
  typedef int64_t product_type;

public:
  FixedPointNumber(const uint8_t int_bits, const uint8_t frac_bits)
    : m_IntegerBits(int_bits), m_FractionBits(frac_bits), m_Value(0)
  {
  }

  FixedPointNumber(const uint8_t int_bits, const uint8_t frac_bits, value_type val)
    : m_IntegerBits(int_bits), m_FractionBits(frac_bits), m_Value(val)
  {
  }

  void AssignInt(int32_t x) { m_Value = x << m_FractionBits; }

  void AssignFloat(float f, bool useRound = true)
  {
    const int32_t one = static_cast<int32_t>((1) << m_FractionBits);
    if (useRound) {
      m_Value = static_cast<value_type>(round(f * one));
    } else {
      m_Value = static_cast<value_type>(ceil(f * one));
    }
  }

  void SetValue(value_type x) { m_Value = x; }

  value_type GetValue() const { return m_Value; }

  double AsDouble() const
  {
    const uint64_t one = 1 << m_FractionBits;
    return static_cast<double>(m_Value) / one;
  }

  uint8_t GetFractionBits() const { return m_FractionBits; }

  uint16_t GetNumBits() const { return m_IntegerBits + m_FractionBits; }

  //--------------------------------------------------------------------------------adam
  // AssignLog2 - run Mike's fast approximation of a log2, as described in his
  // architecture notes.
  //
  void AssignLog2(uint32_t n)
  {
    BOOST_ASSERT(n > 0);
    // The integer part of the log is simply the most significant bit.
    // Builtin_clz returns number of leading (from msb-side) zero-bits.
    const uint32_t msb = 8 * sizeof(n) - __builtin_clz(n) - 1;

    // The fraction part is the next 7 bits below the msb
    uint32_t frac         = n & ((static_cast<uint32_t>(1) << msb) - 1);
    uint32_t bits_in_mask = std::min(static_cast<uint32_t>(7), msb);
    frac >>= msb - bits_in_mask;
    m_Value = (msb << m_FractionBits) | (frac << (m_FractionBits - bits_in_mask));
  }

  //--------------------------------------------------------------------------------adam
  // Multiply: fixed-point multiplication, using an intermediate result twice as
  // large as the operands.
  //
  void Multiply(const FixedPointNumber& x, const FixedPointNumber& y)
  {
    BOOST_ASSERT(8 * sizeof(product_type) >= x.GetNumBits() + y.GetNumBits());

    const uint8_t shift   = (x.GetFractionBits() + y.GetFractionBits() - m_FractionBits);
    product_type  product = (x.GetValue() * y.GetValue()) >> shift;

    // TODO: check for overflow
    m_Value = product;
  }

private:
  uint8_t    m_IntegerBits;
  uint8_t    m_FractionBits;
  value_type m_Value;
};

}  // namespace common
}  // namespace dragenos
#endif  // #ifndef COMMON_FIXED_POINT_HPP
