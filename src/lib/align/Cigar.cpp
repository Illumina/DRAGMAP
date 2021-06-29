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

#include "align/Cigar.hpp"

#include <algorithm>
#include <type_traits>

namespace dragenos {
namespace align {

const char                 Cigar::OPERATION_NAMES[]    = {'M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X'};
const Cigar::OperationCode Cigar::OPERATION_CODES[256] = {INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          SEQUENCE_MATCH,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          DELETE,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          HARD_CLIP,
                                                          INSERT,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          ALIGNMENT_MATCH,
                                                          SKIP,
                                                          INV,
                                                          PAD,
                                                          INV,
                                                          INV,
                                                          SOFT_CLIP,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          SEQUENCE_MISMATCH,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV,
                                                          INV};

Cigar::OperationCode Cigar::getOperationCode(const char operationName)
{
  return OPERATION_CODES[operationName];
}

/**
 * \return position adjustment for the number of database bases not present in the query query
 */
uint32_t Cigar::setOperationSequence(const std::string& operationsSequence, int softClipStart)
{
  operations_.clear();
  const std::string::const_iterator firstNotN =
      std::find_if(
          operationsSequence.cbegin(),
          operationsSequence.cend(),
          [](const char c) { return getOperationCode(c) != SKIP; }) +
      softClipStart;

  if (softClipStart) {
    emplace_back(SOFT_CLIP, softClipStart);
  }
  std::string::const_iterator current = firstNotN;
  while (operationsSequence.end() != current) {
    const OperationCode               operationCode = getOperationCode(*current);
    const std::string::const_iterator next =
        std::find_if(current, operationsSequence.cend(), [current](const char c) { return *current != c; });
    const unsigned length = next - current;
    emplace_back(operationCode, length);
    current = next;
  }

  return std::distance(operationsSequence.cbegin(), firstNotN);
}

uint32_t Cigar::getReferenceLength() const
{
  std::size_t ret = 0;
  for (const auto& op : operations_) {
    if (op.first == ALIGNMENT_MATCH || op.first == DELETE || op.first == SKIP
        /* || op.first == SEQUENCE_MATCH || op.first == SEQUENCE_MISMATCH*/) {
      ret += op.second;
    }
  }

  return ret;
}

//--------------------------------------------------------------------------------adam
// Count how many hard or soft clips are at the start of the cigar
//
int Cigar::countStartClips() const
{
  int rv = 0;

  for (const auto& op : operations_) {
    if (op.first != Cigar::SOFT_CLIP && op.first != Cigar::HARD_CLIP) {
      break;
    }
    rv += op.second;
  }

  return rv;
}

//--------------------------------------------------------------------------------adam
// Count how many hard or soft clips are at the start of the cigar
//
int Cigar::countStartHardClips() const
{
  int rv = 0;

  for (const auto& op : operations_) {
    if (op.first != Cigar::HARD_CLIP) {
      break;
    }
    rv += op.second;
  }

  return rv;
}

//--------------------------------------------------------------------------------adam
// Count how many hard or soft clips are at the end of the cigar
//
int Cigar::countEndClips() const
{
  int rv = 0;

  for (int32_t i = operations_.size() - 1; i >= 0; --i) {
    const auto& op = operations_[i];

    if (op.first != Cigar::SOFT_CLIP && op.first != Cigar::HARD_CLIP) {
      break;
    }
    rv += op.second;
  }

  return rv;
}

int Cigar::countEndHardClips() const
{
  int rv = 0;

  for (int32_t i = operations_.size() - 1; i >= 0; --i) {
    const auto& op = operations_[i];

    if (op.first != Cigar::HARD_CLIP) {
      break;
    }
    rv += op.second;
  }

  return rv;
}

//--------------------------------------------------------------------------------adam
// return the equivalent of GetReferenceLength() + CountEndClips()
//
uint32_t Cigar::getReferenceLengthPlusEndClips() const
{
  return getReferenceLength() + countEndClips();
}

uint32_t Cigar::getClippedLength() const
{
  std::size_t ret = 0;
  for (const auto& op : operations_) {
    if (op.first == Cigar::ALIGNMENT_MATCH || op.first == Cigar::INSERT
        /* || op.first == SEQUENCE_MATCH || op.first == SEQUENCE_MISMATCH*/) {
      ret += op.second;
    }
  }

  return ret;
}

void Cigar::softClipsToHardClips()
{
  for (auto& op : operations_) {
    if (Cigar::SOFT_CLIP == op.first) {
      op.first = Cigar::HARD_CLIP;
    }
  }
}

uint32_t SerializedCigar::getReferenceLength() const
{
  std::size_t ret = 0;
  for (const auto& op : boost::make_iterator_range(operations_, operations_ + cigarLength_)) {
    if (op.first == Cigar::ALIGNMENT_MATCH || op.first == Cigar::DELETE || op.first == Cigar::SKIP
        /* || op.first == SEQUENCE_MATCH || op.first == SEQUENCE_MISMATCH*/) {
      ret += op.second;
    }
  }

  return ret;
}

//--------------------------------------------------------------------------------adam
// Count how many hard or soft clips are at the start of the cigar
//
int SerializedCigar::countStartClips() const
{
  int rv = 0;

  for (const auto& op : boost::make_iterator_range(operations_, operations_ + cigarLength_)) {
    if (op.first != Cigar::SOFT_CLIP && op.first != Cigar::HARD_CLIP) {
      break;
    }
    rv += op.second;
  }

  return rv;
}

//--------------------------------------------------------------------------------adam
// Count how many hard or soft clips are at the start of the cigar
//
int SerializedCigar::countStartHardClips() const
{
  int rv = 0;

  for (const auto& op : boost::make_iterator_range(operations_, operations_ + cigarLength_)) {
    if (op.first != Cigar::HARD_CLIP) {
      break;
    }
    rv += op.second;
  }

  return rv;
}

//--------------------------------------------------------------------------------adam
// Count how many hard or soft clips are at the end of the cigar
//
int SerializedCigar::countEndClips() const
{
  int rv = 0;

  for (int32_t i = cigarLength_ - 1; i >= 0; --i) {
    const auto& op = operations_[i];

    if (op.first != Cigar::SOFT_CLIP && op.first != Cigar::HARD_CLIP) {
      break;
    }
    rv += op.second;
  }

  return rv;
}

int SerializedCigar::countEndHardClips() const
{
  int rv = 0;

  for (int32_t i = cigarLength_ - 1; i >= 0; --i) {
    const auto& op = operations_[i];

    if (op.first != Cigar::HARD_CLIP) {
      break;
    }
    rv += op.second;
  }

  return rv;
}

//--------------------------------------------------------------------------------adam
// return the equivalent of GetReferenceLength() + CountEndClips()
//
int SerializedCigar::getReferenceLengthPlusEndClips() const
{
  return getReferenceLength() + countEndClips();
}

}  // namespace align
}  // namespace dragenos
