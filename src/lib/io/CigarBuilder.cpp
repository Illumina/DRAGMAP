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

#include "io/CigarBuilder.hpp"

/* ============ end of #include directives ============ */
//#include "edico_memdebug.h"

namespace dragenos {
namespace io {

const uint8_t CigarBuilder::MATCH    = 0;
const uint8_t CigarBuilder::INS      = 1;
const uint8_t CigarBuilder::DEL      = 2;
const uint8_t CigarBuilder::REFSKIP  = 3;
const uint8_t CigarBuilder::SOFTCLIP = 4;
const uint8_t CigarBuilder::HARDCLIP = 5;
const uint8_t CigarBuilder::PAD      = 6;
const uint8_t CigarBuilder::SEQMATCH = 7;
const uint8_t CigarBuilder::MISMATCH = 8;

//--------------------------------------------------------------------------------adam
// consructor
//
CigarBuilder::CigarBuilder() : m_ReadStart(0), m_ReadEnd(0), m_RefStart(0), m_RefEnd(0) {}

//--------------------------------------------------------------------------------adam
// GetCurrentOp - get the operation (e.g. match, insert, etc) that we're currently
// tallying
//
uint8_t CigarBuilder::GetCurrentOp() const
{
  if (m_Records.empty()) return 255;

  // Force tallying a new record if we've filled one up.
  if (GetCurrentCount() == 63) return 255;

  const size_t i = m_Records.size() - 1;
  return m_Records[i] & 0xF;
}

//--------------------------------------------------------------------------------adam
// GetCurrentOp - get the number of bases that have been tallied up for the current
// op (e.g. match, insert, del, etc)
//
uint16_t CigarBuilder::GetCurrentCount() const
{
  if (m_Records.empty()) return 0;

  const size_t i = m_Records.size() - 1;
  return m_Records[i] >> 4;
}

//--------------------------------------------------------------------------------adam
// IncrementCurrentOp - add another match/insert/delete to the cigar.  If it's the
// same op as the previous one, increment the number.  Otherwise finish up the current
// cigar record and start a new one for the new op.
//
void CigarBuilder::IncrementCurrentOp(const uint8_t op)
{
  if (op != GetCurrentOp()) {
    m_Records.push_back(0);
  }

  const size_t i     = m_Records.size() - 1;
  uint16_t     count = GetCurrentCount() + 1;
  m_Records[i]       = (count << 4) | (static_cast<uint16_t>(op) & 0xF);
}

//--------------------------------------------------------------------------------adam
// ConsolidateRecords - combine neighboring records of the same type, so long
// as they won't overflow the maximum count per record
//
void CigarBuilder::ConsolidateRecords()
{
  if (!m_Records.size()) return;

  std::vector<uint16_t> new_records;
  new_records.push_back(m_Records[0] & 0xF);
  const uint16_t BUCKET_MAX = 0xFFF;

  for (std::vector<uint16_t>::const_iterator it = m_Records.begin(); it != m_Records.end(); ++it) {
    uint32_t num = (*it) >> 4;
    uint8_t  op  = (*it) & 0xF;

    const uint8_t  current_record_op    = new_records.back() & 0xF;
    const uint32_t current_record_count = new_records.back() >> 4;
    if ((op != current_record_op) || (current_record_count + num > BUCKET_MAX)) {
      new_records.push_back((num << 4) | op);
    } else {
      new_records.back() = ((current_record_count + num) << 4) | op;
    }
  }
  m_Records.swap(new_records);
}

/*
//--------------------------------------------------------------------------------adam
// Trim - cut off the record to ensure it covers only read_len bases
//
void CigarBuilder::Trim(const uint32_t read_len)
{
  int64_t count = 0;
  uint32_t i = 0;
  for ( ; i < m_Records.size() && (count < read_len); ++i ) {
    uint16_t num = m_Records[i] >> 4;
    uint8_t op = m_Records[i] & 0xF;
    if ( op == MATCH || op == DELETION )
      count += num;
  }

  if ( count > read_len ) {
    uint16_t num = m_Records[i] >> 4;

    uint8_t op = m_Records[i] & 0xF;
    uint16_t newnum = read_len - count;
    ASSERT(newnum < num);
    m_Records[i] = (newnum << 4) | op;
  }
  m_Records.resize(i+1);
}
*/

}  // namespace io
}  // namespace dragenos
