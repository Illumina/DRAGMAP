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
#include <memory>
#include <string>
#include <vector>

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

  template <typename IT>
  void getBases(size_t beginPosition, size_t endPosition, IT out) const
  {
    for (auto position = beginPosition; position != endPosition; ++position) {
      const unsigned char base4bpb = getBase(position);
      *out++                       = base4bpb;
      //const unsigned char base2bpb = reference::ReferenceSequence::translateTo2bpb(base4bpb);
      //*out++ = base2bpb;
    }
  }

  template <typename IT>
  void getRcBases(size_t beginPosition, size_t endPosition, IT out) const
  {
    for (auto position = beginPosition; position != endPosition; ++position) {
      const unsigned char base4bpb = getRcBase(position);
      *out++                       = base4bpb;
      //const unsigned char base2bpb = reference::ReferenceSequence::translateToR2bpb(base4bpb);
      //*out++ = base2bpb;
    }
  }

  unsigned char        getBase(size_t position) const;
  unsigned char        getRcBase(size_t position) const;
  const unsigned char* getData() const { return data_; }
  size_t               getSize() const { return size_; }
  /// decode 4 bits into AIUPAC character using only 4 LSB
  static char decodeBase(unsigned char base);
  /// translate into 2 bits encoding using only 4 LSB
  static unsigned char translateTo2bpb(unsigned char base4bpb);
  /// translate into 2 bits reverse complement encoding using only 4 LSB
  static unsigned char translateToR2bpb(unsigned char base4bpb);

private:
  std::vector<Region> trimmedRegions_;
  /// raw data from reference.bin, encoded as 2 bases per byte
  const unsigned char* data_;
  /// number of bytes available in data_
  size_t size_;
};

}  // namespace reference
}  // namespace dragenos

#endif  // #ifndef REFERENCE_REFERENCE_SEQUENCE_HPP
