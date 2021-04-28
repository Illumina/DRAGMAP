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

#ifndef SEQUENCES_CRC_POLYNOMIAL_HPP
#define SEQUENCES_CRC_POLYNOMIAL_HPP

#include <array>
#include <cassert>
#include <cstdint>
#include <string>

namespace dragenos {
namespace sequences {

/**
 ** \brief Encoding of a CRC polynomial
 **
 ** Conceptually, a CRC polynmial is just a bit vector. From a practical point
 ** of view, it is valuable to enable byte operations, which in turns raises the
 ** question of byte and bit order.
 ** This implementation is simply an array of 16 unsigned bytes - enabling up to
 ** 128 bits polynomials - with a little endian bit and byte order.
 ** Unset bits are guaranteed to be 0.
 **/
class CrcPolynomial {
private:
  // This is a table of primitive CRC polynomials up to 128 bits long, with 16 polynomials
  // per length, stored in byte arrays as:  CRC_POLYS[length][number][byte].
  // As an argument to crcHash, CRC_POLYS[L][N] can be used, or even CRC_POLYS[L]
  // to use the first polynomial of that length.  Or favorites can be copied elsewhere.
  static const unsigned char CRC_POLYS[129][16][16];

public:
  /// bit count as specified upon creation
  unsigned       getBitCount() const { return bitCount_; }
  unsigned       getByteCount() const { return (bitCount_ + 7) / 8; }
  const uint8_t* getData() const { return data_.data(); }
  const uint8_t* begin() const { return getData(); }
  const uint8_t* end() const { return begin() + getByteCount(); }

  /// Create a CrcPolynomial rom a byte array
  CrcPolynomial(unsigned bitCount, const uint8_t* data)
    : bitCount_(bitCount), data_(generateData(bitCount, data))
  {
  }

  /// Create a CrcPolynomial from the string representing its hexadecimal value
  CrcPolynomial(unsigned bitCount, const std::string data)
    : bitCount_(bitCount), data_(generateData(bitCount, data))
  {
  }

  CrcPolynomial(unsigned bitCount, unsigned polyIndex)
    : bitCount_(bitCount), data_(generateData(bitCount, CRC_POLYS[bitCount][polyIndex]))
  {
  }

  bool operator==(const CrcPolynomial& other) const { return other.data_ == data_; }

  bool operator==(const std::string& s) const
  {
    CrcPolynomial other(bitCount_, s);
    return (*this == other);
  }

private:
  unsigned                       bitCount_;
  std::array<uint8_t, 16>        data_;
  static std::array<uint8_t, 16> generateData(unsigned bitCount, const uint8_t* data);
  static std::array<uint8_t, 16> generateData(unsigned bitCount, const std::string s);
};

}  // namespace sequences
}  // namespace dragenos

#endif
