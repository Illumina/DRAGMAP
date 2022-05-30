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

#ifndef REFERENCE_REFERENCE_SEQUENCE_HPP
#define REFERENCE_REFERENCE_SEQUENCE_HPP

#include <array>
#include <boost/format.hpp>
#include <memory>
#include <string>
#include <vector>
#ifdef __AVX2__
#include <immintrin.h>
#endif

#include "common/Exceptions.hpp"

namespace dragenos {
namespace reference {

/**
 ** \brief Wrapper for the "reference.bin" produced with DRAGEN hashtables
 **
 ** Values in "reference.bin" can be:
 **   * padding: bytes with 0x00 value
 **   * two bases IUPAC encoded in little endian order
 **
 ** IUPAC encoding is a 4 bits encoding where the bit position for A, C, G and
 ** T are 0, 1, 2 and 3 respectively. This means that we have A:1, C:2, G:4 and
 ** T:8. Other IUPAC values are obtained with a bitwise OR. For instance, N is
 ** A Or C OR G Or T which has the value 0xF.
 **
 ** Litle endian means that:
 **   * the byte at offset j will contain the two bases preceding the bases at
 **     offset j+1
 **   * the 4 least significant bits are for the base preceding the base encoded
 **     in the 4 most significant bits (e.g. "AG" is encoded as 0x41)
 **
 ** For instance, the sequence NNACGTATAGAC will be encoded as the sequence of
 ** bytes {0xFF, 0x21, 0x84, 0x81, 0x41, 0x21}.
 **
 ** The position encoded in the hashtable is simply twice the byte offset in the
 ** file "reference.bin". The parity of the position indicates the position in
 ** the byte (even position for LSB, odd for MSB).
 **
 ** Retrieving a base from the coordinates in the original FASTA reference (i.e.
 ** the 0-based index of the contig and the 0-based position in that contig)
 ** requires additional information regarding the starting position of each
 ** contig in "reference.bin", together with their length and the size of the
 ** trimed regions at the beginning and end of the contigs (sequences with Ns
 ** only).
 **
 ** Converting a hashtable position (or offset in the "reference.bin" file) back
 ** to a position in the original FASTA reference require to take into account
 ** the padding (regions of 0x00 added into "reference.bin") and trimming.
 **/
class ReferenceSequence {
public:
  typedef std::array<uint64_t, 2> Region;
  /**
   ** \brief Constructor
   **
   ** Simply store locally the pointer to the data. Requires the pointer to
   ** stay valid for the lifetime of the ReferenceSequence instance and all its
   ** copies, or until a reset is calld with a different pointer.
   ** It is the responsibility of the calling code to manage the actual data.
   **
   ** It also requires the list of sequences from the Hashtable Config in order
   ** to convert the absolute reference positions into the corresponding address
   ** in the data. The position to address maping needs to know how many N bases
   ** were removed at the beginning and end of each contig.
   **
   ** \param sequences length and trim information for each of the sequences
   ** \param data pointer to the reference sequence encoded as 2 bases per byte
   ** \param size bumber of bytes of data available (half the total number of bases - including padding -
   *stored in the file reference.bin)
   **/
  ReferenceSequence(
      std::vector<Region>  trimmedRegions = std::vector<Region>(),
      const unsigned char* data           = nullptr,
      const size_t         size           = 0)
    : trimmedRegions_(std::move(trimmedRegions)), data_(data), size_(size)
  {
  }

  void reset(std::vector<Region> trimmedRegions, const unsigned char* data, const size_t size)
  {
    trimmedRegions_ = std::move(trimmedRegions);
    data_           = data;
    size_           = size;
  }

  void getBases(size_t beginPosition, size_t endPosition, std::vector<unsigned char>& out) const
  {
    checkPosition(endPosition);

    size_t len = endPosition - beginPosition;
    size_t pos = 0;
    out.clear();
    out.resize(len);

    // handle even position index
    if (beginPosition % 2 == 1) {
      out[pos] = getBaseNoCheck(beginPosition + pos);
      pos++;
    }

#ifdef __AVX2__
    constexpr int  ELEMS_AVX2 = 32;
    unsigned char* dst        = out.data();
    __m128i        mask       = _mm_set1_epi8(0x0F);

    for (; pos + ELEMS_AVX2 <= len; pos += ELEMS_AVX2) {
      __m128i data = _mm_loadu_si128((__m128i*)(&data_[(beginPosition + pos) / 2]));
      __m128i low  = _mm_and_si128(data, mask);
      __m128i high = _mm_and_si128(_mm_srli_epi16(data, 4), mask);

      __m256i lowExt      = _mm256_cvtepu8_epi16(low);
      __m256i highExt     = _mm256_cvtepu8_epi16(high);
      __m256i highShifted = _mm256_slli_si256(highExt, 1);  // shifting out the byte 15 but we don't care

      __m256i resBases = _mm256_or_si256(lowExt, highShifted);
      _mm256_storeu_si256((__m256i*)&dst[pos], resBases);
    }
#endif

    // process remaining bases
    for (; pos != len; pos++) {
      const unsigned char base4bpb = getBaseNoCheck(beginPosition + pos);
      out[pos]                     = base4bpb;
    }
  }

  void getRcBases(size_t beginPosition, size_t endPosition, std::vector<unsigned char>& out) const
  {
    checkPosition(endPosition);

    size_t len = endPosition - beginPosition;
    out.clear();
    out.resize(len);
    for (size_t pos = 0; pos < len; pos++) {
      const unsigned char base4bpb = getRcBaseNoCheck(pos + beginPosition);
      out[pos]                     = base4bpb;
    }
  }

  inline unsigned char getBase(size_t position) const
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
    // has high impact on performance
    checkPosition(position);
    return getBaseNoCheck(position);
  }

  inline unsigned char getRcBase(size_t position) const
  {
    checkPosition(position);
    return getRcBaseNoCheck(position);
  }

  const unsigned char* getData() const { return data_; }
  size_t               getSize() const { return size_; }
  /// decode 4 bits into AIUPAC character using only 4 LSB
  static char decodeBase(unsigned char base);
  /// translate into 2 bits encoding using only 4 LSB
  static unsigned char translateTo2bpb(unsigned char base4bpb);
  /// translate into 2 bits reverse complement encoding using only 4 LSB
  static unsigned char translateToR2bpb(unsigned char base4bpb);

private:
  void checkPosition(size_t position) const
  {
    if (position / 2 >= size_) {
      boost::format message =
          boost::format("position greater than reference size: %i > 2 * %i") % position % size_;
      BOOST_THROW_EXCEPTION(common::InvalidParameterException(message.str()));
    }
  }

  inline unsigned char getBaseNoCheck(size_t position) const
  {
    const unsigned char twoBases = data_[position / 2];
    const bool          msb      = (position % 2);  // use the 4 MSB for odd positions
    return msb ? (twoBases >> 4) : (twoBases & 0xF);
  }

  inline unsigned char getRcBaseNoCheck(size_t position) const
  {
    unsigned char                              b = getBaseNoCheck(position);
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

  std::vector<Region> trimmedRegions_;
  /// raw data from reference.bin, encoded as 2 bases per byte
  const unsigned char* data_;
  /// number of bytes available in data_
  size_t size_;
};

}  // namespace reference
}  // namespace dragenos

#endif  // #ifndef REFERENCE_REFERENCE_SEQUENCE_HPP
