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

#ifndef SEQUENCES_CRC_HASHER_HPP
#define SEQUENCES_CRC_HASHER_HPP

#include <fstream>
#include <iostream>
#include <memory>

#include <inttypes.h>

#include "CrcPolynomial.hpp"

namespace dragenos {
namespace sequences {

class CrcHasher {
public:
  CrcHasher(CrcPolynomial poly);
  unsigned         getBitCount() const { return bitCount_; }
  unsigned         getByteCount() const { return (bitCount_ + 7) / 8; }
  uint64_t         getHash64(uint64_t value) const;
  static void      crcHashSlow(int bitCount, const uint8_t* poly, const uint8_t* data, uint8_t* hash);
  static uint64_t* crcHash64Init(int bitCount, const uint8_t* poly);

private:
  unsigned                    bitCount_;
  std::unique_ptr<uint64_t[]> init64_;
};

}  // namespace sequences
}  // namespace dragenos

#endif  // #ifndef SEQUENCES_CRC_HASHER_HPP
