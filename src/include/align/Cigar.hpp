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

#ifndef ALIGN_CIGAR_HPP
#define ALIGN_CIGAR_HPP

#include <sstream>
#include <vector>

#include <boost/assert.hpp>
#include <boost/format.hpp>
#include <boost/range/iterator_range.hpp>

#include "common/Debug.hpp"

namespace dragenos {
namespace align {

/**
 ** \brief Component representing a BAM CIGAR string
 **/
class Cigar {
public:
  enum OperationCode {
    ALIGNMENT_MATCH   = 0,
    INSERT            = 1,
    DELETE            = 2,
    SKIP              = 3,
    SOFT_CLIP         = 4,
    HARD_CLIP         = 5,
    PAD               = 6,
    SEQUENCE_MATCH    = 7,
    SEQUENCE_MISMATCH = 8,
    INV               = 255
  };
  static const char          OPERATION_NAMES[Cigar::SEQUENCE_MISMATCH + 1];
  static const OperationCode OPERATION_CODES[256];
  static char getOperationName(const OperationCode operationCode) { return OPERATION_NAMES[operationCode]; }
  static OperationCode getOperationCode(
      const char operationName);  // {return OPERATION_CODES[operationName];};

  struct Operation : public std::pair<OperationCode, unsigned> {
    typedef std::pair<OperationCode, unsigned> BaseT;
    Operation(OperationCode f, unsigned s) : BaseT(f, s) {}
    bool operator!=(const Operation& that) const { return that.first != first || that.second != second; }
    bool operator==(const Operation& that) const { return that.first == first && that.second == second; }
  };
  typedef std::vector<Operation> Operations;

  const Operation* getOperations() const { return operations_.data(); }
  /// set the cigar operations from the individual operations in operations sequence string
  uint32_t setOperationSequence(const std::string& operationsSequence, int softClipStart = 0);
  unsigned getNumberOfOperations() const { return operations_.size(); }
  void     clear() { operations_.clear(); }
  void     emplace_back(const OperationCode operationCode, unsigned length)
  {
    operations_.emplace_back(operationCode, length);
  }
  void                 push_back(const Operation& operation) { operations_.push_back(operation); }
  friend std::ostream& operator<<(std::ostream& os, const Cigar& cigar)
  {
    for (const auto& operation : cigar.operations_) {
      os << operation.second << getOperationName(operation.first);
    }
    return os;
  }
  uint32_t getReferenceLength() const;
  int      countStartClips() const;
  int      countStartHardClips() const;

  int      countEndClips() const;
  int      countEndHardClips() const;
  uint32_t getReferenceLengthPlusEndClips() const;
  uint32_t getClippedLength() const;
  void     softClipsToHardClips();

  bool operator!=(const Cigar& that) const
  {
    return that.operations_.size() != operations_.size() ||
           that.operations_.end() !=
               std::mismatch(operations_.begin(), operations_.end(), that.operations_.begin()).second;
  }

  bool operator==(const Cigar& that) const { return !(*this != that); }

  bool empty() const { return operations_.empty(); }

private:
  Operations operations_;
};

class SerializedCigar {
  unsigned short   cigarLength_;
  Cigar::Operation operations_[];

public:
  void operator<<(const Cigar& cigar)
  {
    cigarLength_ = cigar.getNumberOfOperations();
    std::copy(cigar.getOperations(), cigar.getOperations() + cigarLength_, operations_);
  }

  const Cigar::Operation* begin() const { return operations_; }
  const Cigar::Operation* end() const { return operations_ + cigarLength_; }

  friend std::ostream& operator<<(std::ostream& os, const SerializedCigar& cig)
  {
    for (const auto& operation : boost::make_iterator_range(cig.begin(), cig.end())) {
      os << operation.second << Cigar::getOperationName(operation.first);
    }
    return os;
  }

  std::size_t getByteSize() const { return sizeof(SerializedCigar) + cigarLength_ * sizeof(*operations_); }

  static std::size_t getByteSize(const Cigar& cig)
  {
    return sizeof(SerializedCigar) + cig.getNumberOfOperations() * sizeof(*operations_);
  }

  uint32_t getReferenceLength() const;

  int countStartClips()  // includes hard and soft clips
      const;
  int countStartHardClips() const;

  int countEndClips()  // includes hard and soft clips
      const;
  int countEndHardClips() const;

  int getReferenceLengthPlusEndClips()  // includes hard and soft clips
      const;

  bool empty() const { return !cigarLength_; }
};

}  // namespace align
}  // namespace dragenos
#endif  // #ifndef ALIGN_CIGAR_HPP
