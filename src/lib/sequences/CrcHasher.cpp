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

#include <memory.h>
#include <stdlib.h>
#include <iomanip>
#include <iostream>

#include "sequences/CrcHasher.hpp"

namespace dragenos {
namespace sequences {

CrcHasher::CrcHasher(CrcPolynomial poly)
  : bitCount_(poly.getBitCount()), init64_(crcHash64Init(poly.getBitCount(), poly.getData()))
{
}

uint64_t* CrcHasher::crcHash64Init(int bits, const uint8_t* poly)
{
  int       bytes = (bits + 7) >> 3, i, j;
  int       bufQw = (1 + 256 * bytes);
  uint64_t* init  = new uint64_t[bufQw];
  std::fill(init, init + bufQw, 0);
  uint64_t* p = init;
  uint64_t  data;

  if (!init) return NULL;
  // Store the byte count in the init buffer
  // TODO: remove the bytes count from here
  //*p++ = bytes;
  // Store the slow-hash of each byte value 0-255 in each byte position in a data buffer
  for (i = 0; i < bytes; i++) {
    for (j = 0; j < 256; j++) {
      data = (uint64_t)j << (i << 3);
      crcHashSlow(bits, poly, reinterpret_cast<uint8_t*>(&data), reinterpret_cast<uint8_t*>(p++));
    }
  }
  return init;
}

uint64_t CrcHasher::getHash64(const uint64_t value) const
{
  const uint64_t* init = init64_.get();
  const uint8_t*  data = reinterpret_cast<const uint8_t*>(&value);
  uint64_t        hash = 0;
  // TODO: get the number of bytes from a more reasonable location
  //uint64_t bytes = *init++;
  unsigned bytes = getByteCount();
  while (bytes--) {
    hash ^= init[*data++];
    init += 256;
    // if(std::distance((const char*)init64_.get(), (const char*)init) > getByteCount())
    //{
    //  std::cout << "std::distance((const char*)init64_.get(), (const char*)init) " << std::distance((const char*)init64_.get(), (const char*)init) << std::endl;
    //  std::cout << "getByteCount() " << getByteCount() << std::endl;
    //}
  }
  return hash;
}

void CrcHasher::crcHashSlow(int bits, const uint8_t* poly, const uint8_t* data, uint8_t* hash)
{
  int bytes = (bits + 7) >> 3, topByte = bytes - 1;
  int topBitMask = (1 << ((bits + 7) % 8)), topByteMask = ((topBitMask << 1) - 1);
  int i, j, subtract;

  // Since data and polynomial are the same length, copy in all the data bytes immediately.
  // This data order doesn't match normal CRC computation, which by processing byte zero first,
  // effectively treats it as the most significant byte of the dividend.  But with odd bit
  // lengths, this works better.
  memcpy(hash, data, bytes);
  // Loop through the bits
  for (i = 0; i < bits; i++) {
    // Plan to subtract the polynomial if the MSB is 1
    subtract = (hash[topByte] & topBitMask);

    // Left-shift the remainder (corresponds to right-shifting the polynomial position)
    for (j = topByte; j > 0; j--) hash[j] = (hash[j] << 1) | (hash[j - 1] >> 7);
    hash[0] <<= 1;

    // Subtract the polynomial if required to cancel the MSB shifted out
    if (subtract)
      for (j = 0; j < bytes; j++) hash[j] ^= poly[j];
  }
  // Mask off unused positions in the top byte
  hash[topByte] &= topByteMask;
}

}  // namespace sequences
}  // namespace dragenos
