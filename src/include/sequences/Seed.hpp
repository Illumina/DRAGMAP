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

#ifndef SEQUENCES_SEED_HPP
#define SEQUENCES_SEED_HPP

#include <ostream>
#include <vector>

#include "sequences/Read.hpp"

namespace dragenos {
namespace sequences {

/**
 ** \brief A region along a read used to generate hash keys to query the reference
 ** hashtable
 **
 ** Note that the Seed only encapsulates the concept of a region along the read
 ** in a minimalistic way. Particularly, the concepts of orientation and of seed
 ** extension are mostly left outside of this class. The orientation is supported
 ** in the sense that the Seed methods enable generating the sequence of bases both
 ** in the forward and in the reverse direction. Similarly, seed extensions are
 ** supported in the sense that the Seed methods enable generating the sequences
 ** of bases for any extension wing that the user might need. However, it is the
 ** responsibility of the client code to manage the direction and extension sizes.
 **/
class Seed {
public:
  static const uint32_t DEFAULT_PERIOD       = 2;
  static const uint32_t DEFAULT_PATTERN      = 0x01;
  static const uint8_t  DEFAULT_FORCE_LAST_N = 1;
  typedef uint64_t      Data;
  /// check for Ns in the primary seed
  static bool isValid(const Read& read, unsigned readPosition, unsigned primaryLength);
  /// Constructor required 0 < primaryLength <= 32 (for encoding on 64 bits)
  Seed(const Read* read, unsigned readPosition, unsigned primaryLength);
  /**
   ** \brief returns the list of seed offsets for a read of the given read length
   **
   ** TODO: verify that the vector allocation is acceptable, otherwise, pass it as a reference parameter
   */
  static std::vector<size_t> getSeedOffsets(
      const size_t   readLength,
      const unsigned length,
      const uint32_t period     = DEFAULT_PERIOD,
      const uint32_t pattern    = DEFAULT_PATTERN,
      const uint8_t  forceLastN = DEFAULT_FORCE_LAST_N);
  /**
   ** \brief the encoded data for the primary seed in the format required by the CrCHasher
   **
   ** Each base is encoded on 2 bits with C=1, G=2, T=3 and all other bases encoded to 0.
   ** The encoding is little endian which means that bits [2i+1:2i] will hold the encoding
   ** for the base at position readPosition + i.
   ** If the reverseComplement flag is set, the resulting sequence is then reverse complemented
   ** and stays in the 2*primaryLength least significant bits of the resulting value.
   **
   ** \pre isValid() is true - so that all required bases are within the read
   **
   ** \param reverseComplement the orientation of the required data
   **/
  Data getPrimaryData(bool reverseComplement) const;
  /**
   ** \brief encode data for the extention wings in the format required by the CrCHasher
   **
   ** The wings are the regions defined as:
   ** * before: [readPosition - toHalfExtension, readPosition - fromHalfExtension)
   ** * after: [readPosition + primaryLength + fromHalfExtension, readPosition + primaryLength +
   *toHalfExtension)
   ** The bases are encoded as with the primary data (ACGT -> 0123 respectively,
   ** otherwise 0). They are stored in a little endian sequence (wing before then after).
   ** If the reverseComplement flag is set, the resulting sequence is reverse complemented.
   **
   ** \pre isValid(toHalfExtension) is true - so that all required bases are within the read
   **
   ** \param fromHalfExtension the reach of the last extension wings (0 if none)
   ** \param toHalfExtension the reach of the new extension wings
   ** \param reverseComplement the orientation of the required data
   **/
  Data getExtendedData(unsigned fromHalfExtension, unsigned toHalfExtension, bool reverseComplement) const;
  const Read* getRead() const { return read_; }
  unsigned    getReadPosition() const { return readPosition_; }
  unsigned    getPrimaryLength() const { return primaryLength_; }

  /**
   ** \brief get the read position of the last base of the seed when applying the specified half extension
   *length
   **/
  unsigned getFirstBaseReadPosition(const unsigned halfExtensionLength) const
  {
    return getReadPosition() - halfExtensionLength;
  }

  /**
   ** \brief get the read position of the last base of the seed when applying the specified half extension
   *length
   **/
  unsigned getLastBaseReadPosition(const unsigned halfExtensionLength) const
  {
    return getReadPosition() + getPrimaryLength() + halfExtensionLength - 1;
  }
  /**
   ** \brief check that the validity of seed with the specified half extensing length
   **
   ** The seed is valid if all its bases are within the underlying read:
   ** halfExtensionLength <= readPosition
   ** readPosition + primaryLength + halfExtensionLength <= readLength
   **/
  bool isValid(const unsigned halfExtensionLength) const
  {
    return halfExtensionLength <= getReadPosition() &&
           (getReadPosition() + getPrimaryLength() + halfExtensionLength <= getRead()->getLength());
  }

  static Data generateReverseComplement(Data data, const unsigned baseCount);

  bool operator==(const Seed& that) const
  {
    return read_ == that.read_ && readPosition_ == that.readPosition_ &&
           primaryLength_ == that.primaryLength_;
  }

private:
  const Read* read_;
  unsigned    readPosition_;
  unsigned    primaryLength_;

  friend std::ostream& operator<<(std::ostream& os, const Seed& s)
  {
    return os << "Seed(" << s.readPosition_ << "," << s.primaryLength_ << ")";
  }
};

}  // namespace sequences
}  // namespace dragenos

#endif  // #ifndef SEQUENCES_SEED_HPP
