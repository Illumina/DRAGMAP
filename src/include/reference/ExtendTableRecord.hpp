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

#ifndef REFERENCE_EXTEND_TABLE_RECORD_HPP
#define REFERENCE_EXTEND_TABLE_RECORD_HPP

namespace dragenos {
namespace reference {

class ExtendTableRecord {
public:
  enum class LiftCode {
    NONE    = 0,
    ALT     = 1,
    PRI     = 2,
    DIF_PRI = 3,
  };
  ExtendTableRecord(uint64_t value) : value_(value) {}
  uint64_t getValue() const { return value_; }
  uint32_t getPosition() const { return value_ & 0xFFFFFFFF; }
  bool     isReverseComplement() const { return (value_ >> 32) & 1; }
  LiftCode getLiftCode() const { return LiftCode((value_ >> 33) & 3); }
  uint32_t getLiftGroup() const { return (value_ >> 35) & 0xFFFFFFF; }

private:
  uint64_t value_;
};

}  // namespace reference
}  // namespace dragenos

#endif  // #ifndef REFERENCE_EXTEND_TABLE_RECORD_HPP
