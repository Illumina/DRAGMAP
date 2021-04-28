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

#ifndef IO_CIGAR_BUILDER_HPP
#define IO_CIGAR_BUILDER_HPP

#include <algorithm>  // for std::reverse
#include <iostream>
#include <vector>

#include <boost/assert.hpp>

//#include "dragen_exception.hpp"

namespace dragenos {
namespace io {

//--------------------------------------------------------------------------------adam
// CigarBuilder - a class to build up a cigar string one operation at a time
//
class CigarBuilder {
private:
  static const uint8_t MATCH;
  static const uint8_t INS;
  static const uint8_t DEL;
  static const uint8_t REFSKIP;
  static const uint8_t SOFTCLIP;
  static const uint8_t HARDCLIP;
  static const uint8_t PAD;
  static const uint8_t SEQMATCH;
  static const uint8_t MISMATCH;

public:
  CigarBuilder();

  void AddMatch() { IncrementCurrentOp(MATCH); }

  void AddInsertion() { IncrementCurrentOp(INS); }

  void AddDeletion() { IncrementCurrentOp(DEL); }

  void AddSoftClip() { IncrementCurrentOp(SOFTCLIP); }

  void AddHardClip() { IncrementCurrentOp(HARDCLIP); }

  friend std::ostream& operator<<(std::ostream& os, const CigarBuilder& cigar);

  void Reverse()
  {  // reverse the order of the CIGAR records
    std::reverse(m_Records.begin(), m_Records.end());
  }

  void ConsolidateRecords();

  const std::vector<uint16_t>* GetRecords() { return &m_Records; }

  void SetReadStart(const int32_t read_start) { m_ReadStart = read_start; }

  int32_t GetReadStart() const { return m_ReadStart; }

  void SetReadEnd(const int32_t read_end) { m_ReadEnd = read_end; }

  int32_t GetReadEnd() const { return m_ReadEnd; }

  void SetRefStart(const int64_t ref_start) { m_RefStart = ref_start; }

  int64_t GetRefStart() const { return m_RefStart; }

  void SetRefEnd(const int64_t ref_end) { m_RefEnd = ref_end; }

  int64_t GetRefEnd() const { return m_RefEnd; }

  size_t GetLen() const { return m_Records.size() * sizeof(uint16_t); }

  //  void Trim(const uint32_t read_len)
  //    ;

private:
  void IncrementCurrentOp(const uint8_t op);

  uint8_t GetCurrentOp() const;

  uint16_t GetCurrentCount() const;

private:
  std::vector<uint16_t> m_Records;
  int32_t               m_ReadStart;
  int32_t               m_ReadEnd;
  int64_t               m_RefStart;
  int64_t               m_RefEnd;
};

inline std::ostream& operator<<(std::ostream& os, const CigarBuilder& cigar)
{
  for (std::vector<uint16_t>::const_iterator it = cigar.m_Records.begin(); it != cigar.m_Records.end();
       ++it) {
    uint16_t num = (*it) >> 4;
    uint8_t  op  = (*it) & 0xF;
    os << num;
    if (op == CigarBuilder::MATCH)
      os << "M";
    else if (op == CigarBuilder::INS)
      os << "I";
    else if (op == CigarBuilder::DEL)
      os << "D";
    else if (op == CigarBuilder::REFSKIP)
      os << "N";
    else if (op == CigarBuilder::SOFTCLIP)
      os << "S";
    else if (op == CigarBuilder::HARDCLIP)
      os << "H";
    else if (op == CigarBuilder::PAD)
      os << "P";
    else if (op == CigarBuilder::SEQMATCH)
      os << "=";
    else if (op == CigarBuilder::MISMATCH)
      os << "X";
    else
      BOOST_ASSERT(false);
  }
  //  os << ", TOTAL BASES: " << total;
  return os;
}

}  // namespace io
}  // namespace dragenos

#endif  // #ifndef IO_CIGAR_BUILDER_HPP
