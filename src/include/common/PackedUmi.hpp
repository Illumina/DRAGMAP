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

#ifndef COMMON_PACKED_UMI_HPP
#define COMMON_PACKED_UMI_HPP

#include <algorithm>
#include <array>
#include <boost/assert.hpp>
#include <cinttypes>
#include <cstddef>
#include <cstring>
#include <string>

#pragma pack(1)

namespace dragenos {
namespace common {
//--------------------------------------------------------------------------------adamb
// A class representing a bit-packed UMI, attached to DBAM records in a blob
// for quick access and use in duplicate-marking and read-collapsing
//
class PackedUmi {
public:
  enum { MAX_UMI_LEN = 8 };

public:
  typedef uint32_t PackedSeq_t;

private:
  //---------------------------------------------------------------------------
  // A structure packing the lengths of the pair of UMI sequences into one byte
  struct SeqLen_t {
  public:
    SeqLen_t& operator=(const SeqLen_t& other)
    {
      m_lens = other.m_lens;
      return *this;
    }

  public:
    SeqLen_t(const uint8_t l1, const uint8_t l2) : m_lens((l1 & 0xFF) | ((l2 & 0xFF) << 4)) {}

    uint8_t getLen(const uint16_t whichSeq) const { return (m_lens >> (4 * whichSeq)) & 0x0F; }

  private:
    uint8_t m_lens;  // lengths, packed 4-bytes per length.
  };

public:
  PackedUmi() : m_lens(0, 0), m_seqs{0, 0}, m_isN{0xFF, 0xFF} {};

  PackedUmi(
      const char*  seq1,  // UMI sequence for first end
      const size_t len1,  // length of first-end UMI sequence
      const char*  seq2,  // UMI sequence for second end
      const size_t len2)  // length of second-end UMI sequence
      ;

  PackedUmi(const PackedSeq_t ps1, const PackedSeq_t ps2);

  PackedUmi& operator=(const PackedUmi& other)
  {
    m_lens = other.m_lens;
    m_seqs = other.m_seqs;
    m_isN  = other.m_isN;
    return *this;
  }

  uint8_t getLen(const uint32_t whichSeq) const
  {
    return std::min(static_cast<const uint8_t>(MAX_UMI_LEN), m_lens.getLen(whichSeq));
  }

  void getAsciiSeq(const uint32_t whichSeq, char* ascii_out, const uint8_t maxLen) const
  {
    getAsciiSeq(ascii_out, std::min(maxLen, getLen(whichSeq)), m_seqs[whichSeq], m_isN[whichSeq]);
  }

  PackedSeq_t getPackedSeq(const uint32_t whichSeq) const
  {
    return packSeq(getLen(whichSeq), m_isN[whichSeq], m_seqs[whichSeq]);
  }

  uint64_t asUint64() const
  {
    return (static_cast<uint64_t>(getPackedSeq(0)) << 32) | (static_cast<uint64_t>(getPackedSeq(1)));
  }

  template <class OutputStream>
  void getTag(OutputStream& os) const
  {
    os.write("RXZ", 3);
    print(os);
  }

  template <class OutputStream>
  void print(OutputStream& os) const
  {
    char ascii[MAX_UMI_LEN];
    for (auto i = 0; i < 2; ++i) {
      getAsciiSeq(i, ascii, sizeof(ascii));
      os.write(ascii, getLen(i));
      if (i == 0) os.put('-');
    }
  }

  bool operator==(const PackedUmi& other) const
  {
    return (
        (getLen(0) == other.getLen(0)) && (getLen(1) == other.getLen(1)) &&
        (getPackedSeq(0) == other.getPackedSeq(0)) && (getPackedSeq(1) == other.getPackedSeq(1)));
  }

  // Compare packed UMI's lexicographically.  We need this to mimic ReCo's ordering
  // of families alphabetically by sequence.
  bool operator<(const PackedUmi& other) const
  {
    char mine[MAX_UMI_LEN];
    char theirs[MAX_UMI_LEN];
    for (uint32_t i = 0; i < 2; ++i) {
      memset(mine, 0, MAX_UMI_LEN);
      memset(theirs, 0, MAX_UMI_LEN);
      getAsciiSeq(i, mine, MAX_UMI_LEN);
      other.getAsciiSeq(i, theirs, MAX_UMI_LEN);

      const int comparison = strncmp(mine, theirs, std::min(getLen(i), other.getLen(i)));
      if (comparison < 0)
        return true;
      else if (comparison > 0)
        return false;
    }

    return false;
  }

  bool isReverseOf(const PackedUmi& other) const
  {
    return (
        (getLen(0) == other.getLen(1)) && (getLen(1) == other.getLen(0)) &&
        (getPackedSeq(0) == other.getPackedSeq(1)) && (getPackedSeq(1) == other.getPackedSeq(0)));
  }

  static void packedSeqToAscii(char* ascii_out, const uint8_t maxLen, const PackedSeq_t packedSeq)
  {
    const uint32_t len = (packedSeq >> 24) & 0xFF;
    BOOST_ASSERT_MSG(
        maxLen > len, "ERROR: must reserve enough space to hold the UMI sequence plus a null-terminator");
    const uint16_t seq = packedSeq & 0xFFFF;
    const uint8_t  isN = (packedSeq >> 16) & 0xFF;
    return getAsciiSeq(ascii_out, len, seq, isN);
  }

  static PackedSeq_t asciiToPackedSeq(const std::string& seq)
  {
    uint16_t     packedSeq = 0;
    uint8_t      isN       = 0;
    const size_t len       = seq.size();
    transformSequence(seq.c_str(), len, &packedSeq, &isN);
    return packSeq(len, isN, packedSeq);
  }

  bool isValid() const { return getLen(1) && getLen(2); }

  void reverseComplement();

  bool arrangeEndsLexicographically();

  bool lexicographicallyPrecedes(const PackedUmi& other) const;

private:
  static PackedSeq_t packSeq(const uint16_t len, const uint8_t isN, const uint16_t seq)
  {
    return (static_cast<const PackedSeq_t>(len) << 24) | (static_cast<const PackedSeq_t>(isN) << 16) |
           (static_cast<const PackedSeq_t>(seq));
  }

  void transformSequence(       // transform one of the two member sequences
      const uint16_t whichSeq,  // 0 or 1
      const char*    seq,       // input ASCII FASTQ-format sequence
      const size_t   len)         // how many bytes in #seq#
      ;

  static void transformSequence(  // transform ASCII FASTQ-format sequence to DBAM
      const char*  seq,           // input ASCII FASTQ-format sequence
      const size_t len,           // how many bytes in #seq#
      uint16_t*    dest,          // where to put the packed sequence bits
      uint8_t*     isN)               // where to put the 'isN' bits
      ;

private:
  static void getAsciiSeq(char* ascii_out, const uint8_t maxLen, uint16_t seq, const uint8_t isN)
  {
    static char DBAM2FASTQ_BASES[] = "ACGT";
    uint8_t     i                  = 0;
    for (i = 0; i < maxLen; ++i, seq >>= 2) {
      if (isN & (1 << i))
        ascii_out[i] = 'N';
      else
        ascii_out[i] = DBAM2FASTQ_BASES[seq & 0x3];
    }
    ascii_out[i] = 0;
  }

private:
  SeqLen_t                m_lens;  // how long are the two UMI sequences
  std::array<uint16_t, 2> m_seqs;  // the UMI sequences, in DBAM format
  std::array<uint8_t, 2>  m_isN;   // bit per base, whether it's an N
};
#pragma pack()

struct PackedUmiHashkeyCalculator {
  size_t operator()(const PackedUmi& pu) const { return std::hash<unsigned long long>()(pu.asUint64()); }
};

}  // namespace common
}  // namespace dragenos

#endif  // #ifndef COMMON_PACKED_UMI_HPP
