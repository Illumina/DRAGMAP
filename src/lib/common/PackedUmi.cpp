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

#include "common/PackedUmi.hpp"

#include <string.h>
#include <boost/assert.hpp>

//#include "dragen_exception.hpp"

namespace dragenos {
namespace common {

//--------------------------------------------------------------------------------adamb
// constructor -- extract the sequences & lengths into our internal, binary, packed
// format
//
PackedUmi::PackedUmi(const char* seq1, const size_t len1, const char* seq2, const size_t len2)
  : m_lens(len1, len2), m_seqs{0, 0}, m_isN{0, 0}
{
  BOOST_ASSERT(m_lens.getLen(0) == len1);
  BOOST_ASSERT(m_lens.getLen(1) == len2);
  BOOST_ASSERT_MSG((len1 <= 8) && (len2 <= 8), "ERROR: max support UMI length is 8bp");
  transformSequence(0, seq1, len1);
  transformSequence(1, seq2, len2);
}

//--------------------------------------------------------------------------------adamb
// iniatialize from two pre-packed binary UMIs
//
PackedUmi::PackedUmi(const PackedSeq_t ps1, const PackedSeq_t ps2)
  : m_lens(ps1 >> 24, ps2 >> 24),
    m_seqs{static_cast<uint16_t>(ps1 & 0xFFFF), static_cast<uint16_t>(ps2 & 0xFFFF)},
    m_isN{static_cast<uint8_t>((ps1 >> 16) & 0xFF), static_cast<uint8_t>((ps2 >> 16) & 0xFF)}
{
}

//--------------------------------------------------------------------------------adamb
// Take an input ASCII sequence, transform it into our packed representation,
// and store it in this object in one of the two sequence slots.
//
void PackedUmi::transformSequence(
    const uint16_t whichSeq,  // 0 or 1
    const char*    seq,       // input ASCII FASTQ-format sequence
    const size_t   len)         // how many bytes in #seq#
{
  transformSequence(seq, len, &m_seqs[whichSeq], &m_isN[whichSeq]);
}

//--------------------------------------------------------------------------------adamb
// Transform #seq# (which is #len# bases long) into (two bits per base)
// DBAM format in #dest#.
//
void PackedUmi::transformSequence(
    const char*  seq,   // input ASCII FASTQ-format sequence
    const size_t len,   // how many bytes in #seq#
    uint16_t*    dest,  // where to put the packed sequence bits
    uint8_t*     isN)       // where to put the 'isN' bits
{
  // Make sure the destinations are cleared
  *dest = 0;
  *isN  = 0;

  // Assume that the sequence was initialized to zero before this call
  for (size_t i = 0; i < len; ++i) {
    // Here is the translation table we are shooting for:
    //
    //    A (0x41)  =>  0
    //    C (0x43)  =>  1
    //    G (0x47)  =>  2
    //    T (0x54)  =>  3
    //
    // Turns out we can get this like this:
    //
    //      ((x>>1) ^ (x>>2)) & 3
    //
    //   A: ((0x20) ^ (0x10)) & 3 = 0
    //   C: ((0x21) ^ (0x10)) & 3 = 1
    //   G: ((0x23) ^ (0x11)) & 3 = 2
    //   T: ((0x2A) ^ (0x15)) & 3 = 3
    //
    const uint8_t& ascii_base = seq[i];
    const uint16_t dbam_base  = ((ascii_base >> 1) ^ (ascii_base >> 2)) & 3;
    const uint8_t  shift      = i * 2;
    *dest |= (dbam_base << shift);

    if (ascii_base == 'N') (*isN) |= (1 << i);
  }
}

//--------------------------------------------------------------------------------adamb
// Replace the sequences with their reverse-complements.
//
void PackedUmi::reverseComplement()
{
  char s[MAX_UMI_LEN];

  for (size_t i = 0; i < 2; ++i) {
    // Get the original sequence in ASCII format:
    const uint16_t len = getLen(i);
    getAsciiSeq(i, s, MAX_UMI_LEN);

    // Reverse it:
    std::reverse(s, s + len);

    // Pack it back into m_seq and m_isN:
    transformSequence(s, len, &m_seqs[i], &m_isN[i]);

    // And complement it, taking advantage of DBAM format's nice property
    // that a twos-complement of the sequence is the sequence's biological
    // complement.
    m_seqs[i] = (~m_seqs[i]);
  }
}

//--------------------------------------------------------------------------------adamb
// Arrange the two ends of the umi so that the one with the "first" lexicographic
// ASCII representation is sequence 0, and the other one is sequence 1
// Returns true if the ends were swapped, false otherwise.
//
bool PackedUmi::arrangeEndsLexicographically()
{
  char           s0[MAX_UMI_LEN];
  const uint16_t len0 = getLen(0);
  getAsciiSeq(0, s0, MAX_UMI_LEN);

  char           s1[MAX_UMI_LEN];
  const uint16_t len1 = getLen(1);
  getAsciiSeq(1, s1, MAX_UMI_LEN);

  // If the first sequence comes first alphabetically, then do nothing.
  if (strncmp(s0, s1, std::min(len0, len1)) <= 0) return false;

  // The second sequence comes first alphabetically.  Swap the two ends.
  m_seqs[0] = m_seqs[1] = 0;
  transformSequence(s0, len0, &m_seqs[1], &m_isN[1]);
  transformSequence(s1, len1, &m_seqs[0], &m_isN[0]);
  m_lens = SeqLen_t(len1, len0);
  return true;
}

//--------------------------------------------------------------------------------adamb
// Return true if this packed UMI's ASCII representation precedes that of the
// #other# UMI
//
bool PackedUmi::lexicographicallyPrecedes(const PackedUmi& other) const
{
  for (auto i = 0; i < 2; ++i) {
    char           mine[MAX_UMI_LEN];
    const uint16_t mylen = getLen(i);
    getAsciiSeq(i, mine, MAX_UMI_LEN);

    char           theirs[MAX_UMI_LEN];
    const uint16_t theirlen = other.getLen(i);
    other.getAsciiSeq(i, theirs, MAX_UMI_LEN);

    if (strncmp(mine, theirs, std::min(mylen, theirlen)) < 0) {
      return true;
    }
  }

  return false;
}

}  // namespace common
}  // namespace dragenos
