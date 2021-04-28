#define SEQUENCES_CRC_POLYNOMIAL_HPP

#include <array>
#include <cassert>
#include <cstdlib>
#include <string>

namespace dragenos {
namespace sequences {

class CrcPolynomial {
public:
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

private:
  unsigned                bitCount_;
  std::array<uint8_t, 16> data_;

  std::array<uint8_t, 16> generateData(unsigned bitCount, const uint8_t* data)
  {
    std::array<uint8_t, 16> result;
    *(reinterpret_cast<__uint128_t*>(result.data())) = *(reinterpret_cast<const __uint128_t*>(data));
    return result;
  }

  std::array<uint8_t, 16> generateData(unsigned bitCount, std::string s)
  {
    std::array<uint8_t, 16> result;
    result.fill(0);
    while (!s.empty() && ('0' == s.front())) s.erase(s.begin());
    if (s.empty()) return result;
    assert(s.size() < 16);  // support only up to 64 bits
    __uint128_t data = strtoull(s.c_str(), NULL, 16);
    assert(0 != data);
    return generateData(bitCount, reinterpret_cast<uint8_t*>(&data));
  }
};

}  // namespace sequences
}  // namespace dragenos
