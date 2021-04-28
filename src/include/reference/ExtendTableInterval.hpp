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

#ifndef REFERENCE_EXTEND_TABLE_INTERVAL_HPP
#define REFERENCE_EXTEND_TABLE_INTERVAL_HPP

#include <vector>

#include "reference/HashRecord.hpp"

namespace dragenos {
namespace reference {

/**
 ** \brief an interval in the extend table
 **
 ** Initially represented as a set of 1, 2 or 3 HashRecords in the hash table,
 ** the interval is then converted int an instance of a class that encapdulates
 ** the definition of the interval (start, length and count of extra liftover
 ** matches) and additional characteristics to enable various sampling strategies
 ** in the interval, according to the outcome of the mapping.
 **
 ** Expected combinations are:
 ** - SL0
 ** - SL1, S
 ** - SLE(exlifts>0, MSB=0)
 ** - SLE(exlifts>0, MSB=1), S
 ** - SLE(exlifts=0, MSB=0), L
 ** - SLE(exlifts=0, MSB=1), S, L
 ** - S, L
 **
 ** Start position:
 ** INTERVAL_SL0:  The 15-bit Start Field contains bits 14:0 of the Start Value
 ** INTERVAL_SL1:  The 8-bit Start Field contains bits 31:24 of the Start Value
 ** INTERVAL_SLE:  The 8-bit Start Field contains bits 7:0 of the Start Value when flag Msb=0, or bits 31:24
 *of the Start Value when Msb=1
 ** INTERVAL_S:  The 24-bit Start Field contains bits 23:0 of the Start Value.  Furthermore, flag Cry is a
 *carry bit that adds into Start Value bits 31:24, which may already have been populated by INTERVAL_SLE with
 *Msb=1.
 **
 ** Length:
 ** INTERVAL_SL0:  The 9-bit Length Field contains bits 8:0 of the Length Value
 ** INTERVAL_SL1:  The 16-bit Length Field contains bits 15:0 of the Length Value
 ** INTERVAL_SLE:  The 8-bit Length Field contains bits 7:0 of the Length Value when Exlifts > 0, or bits
 *31:24 of the Length Value when Exlifts=0
 ** INTERVAL_L:  The 24-bit Length Field contains bits 23:0 of the Length Value
 **
 **/
class ExtendTableInterval {
public:
  template <typename I>
  ExtendTableInterval(I begin, I end);
  uint32_t getStart() const { return start_; }
  uint32_t getLength() const { return length_; }
  uint8_t  getExtraLiftoverMatchCount() const { return extraLiftoverMatchCount_; }

private:
  uint32_t start_;
  uint32_t length_;
  uint8_t  extraLiftoverMatchCount_;
  /// checks if the range [begin, end) contains a valid combination for an interval set
  template <typename I>
  static bool isValid(I begin, I end);
  /// validate the range [begin, end) end extract the start value
  template <typename I>
  static uint32_t extractStart(I begin, I end);
  /// extract the length (no validation)
  template <typename I>
  static uint32_t extractLength(I begin, I end);
  /// extract the exlift (no validation)
  template <typename I>
  static uint8_t extractExtraLiftoverMatchCount(I begin, I end);
  /// throws an exception reporting unexpected hash record
  static void unexpectedHashRecord(const HashRecord& hashRecord);
};  // class ExtendTableInterval

}  // namespace reference
}  // namespace dragenos
#endif  // #ifndef REFERENCE_EXTEND_TABLE_INTERVAL_HPP
