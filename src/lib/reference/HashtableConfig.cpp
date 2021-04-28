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

#include "reference/HashtableConfig.hpp"
#include <algorithm>
#include <boost/format.hpp>
#include <numeric>

namespace dragenos {
namespace reference {

HashtableConfig::HashtableConfig(const char* const data, const size_t size)
  : header_(*header(data, size)),
    sequences_(sequences(data, size)),
    sequenceNames_(sequenceNames(data, size)),
    hostVersion_(stringAddress(data, &hostVersion_)),
    commandLine_(stringAddress(data, &commandLine_)),
    refInput_(stringAddress(data, &refInput_)),
    refOutput_(stringAddress(data, &refOutput_)),
    refIdxFname_(refIdxFname(data, size)),
    refFname_(stringAddress(data, &refFname_)),
    hashFname_(stringAddress(data, &hashFname_)),
    altLiftover_(stringAddress(data, &altLiftover_))
//  , referenceDictionary_(referenceDictionary())
{
}

unsigned long HashtableConfig::getReferenceLength() const
{
  return std::accumulate(sequences_.begin(), sequences_.end(), 0UL, [](unsigned long t, const Sequence& s) {
    return t + s.seqLen;
  });
}

std::vector<HashtableConfig::Region> HashtableConfig::getTrimmedRegions() const
{
  std::vector<Region> trimmedRegions;
  for (const auto& sequence : sequences_) {
    trimmedRegions.push_back({sequence.seqStart, sequence.seqStart + sequence.begTrim});
    trimmedRegions.push_back(
        {sequence.seqStart + sequence.seqLen - sequence.endTrim, sequence.seqStart + sequence.seqLen});
  }
  return trimmedRegions;
}

std::pair<size_t, int64_t> HashtableConfig::convertToReferenceCoordinates(int64_t position) const
{
  auto sequence = std::lower_bound(
      sequences_.begin(),
      sequences_.end(),
      position,
      [](const HashtableConfig::Sequence& sq, int64_t position) {
        const auto positionRange = getPositionRange(sq);
        return positionRange.first < position;
      });

  size_t sequenceIndex = std::distance(sequences_.begin(), sequence);
  if (sequences_.begin() != sequence) {
    --sequence;
    --sequenceIndex;
  }
  while (sequences_.end() != sequence) {
    if (sequence->seqStart > position) {
      //      boost::format message = boost::format("reference position outside of valid range: current sequence start=%i: queried position=%i") % sequence->seqStart % position;
      // TODO: verify that the position can actually be less than the start of the sequence (because of
      // padding)
      //      BOOST_THROW_EXCEPTION(common::InvalidParameterException(message.str()));
      const int64_t offset = sequence->begTrim + (position - sequence->seqStart);

      /////////////////
      //std::cerr << "position before beginning of sequence:" << sequence.begTrim << " + (" << position << " - " <<sequence.seqStart << ") = " << offset << std::endl;
      /////////////////

      return std::pair<size_t, int64_t>(sequenceIndex, offset);
    }

    const auto positionRange = getPositionRange(*sequence);
    if ((positionRange.first <= position) && (positionRange.second >= position)) {
      const int64_t offset = sequence->begTrim + (position - sequence->seqStart);
      return std::pair<size_t, int64_t>(sequenceIndex, offset);
    }

    ++sequenceIndex;
    ++sequence;
  }

  const auto lastRange = getPositionRange(sequences_.back());
  //  boost::format message = boost::format("reference position greater than last sequence position: last sequence end=%i: queried position=%i") % lastRange.second % position;
  // TODO: verify that the position can actually be greater than the end of the sequence (because of padding)
  //  BOOST_THROW_EXCEPTION(common::InvalidParameterException(message.str()));
  const auto&   last   = sequences_.back();
  const int64_t offset = last.begTrim + (position - last.seqStart);
  /////////////////
  //std::cerr << "position after beginning of last sequence:" << sequence.begTrim << " + (" << position << " - " <<sequence.seqStart << ") = " << offset << std::endl;
  /////////////////

  return std::pair<size_t, int64_t>(sequences_.size() - 1, offset);
}

#if 0
ReferenceDictionary HashtableConfig::referenceDictionary() const
{
  assert(sequences_.size() == header_.numRefSeqs);
  assert(sequenceNames_.size() == header_.numRefSeqs);
  ReferenceDictionary ret;
  uint64_t totalTrimmed = 0;
  uint64_t globalPos = 0;
  size_t   padLen     = REF_SEQ_END_PAD_BASES;
  uint32_t seqLenNoNs = 0;
  for (size_t i = 0; i < header_.numRefSeqs; ++i) {
    globalPos += padLen;
    ret.add_entry(
        i,
        sequenceNames_[i],
        sequences_[i].seqStart,
        sequences_[i].begTrim,
        sequences_[i].endTrim,
        sequences_[i].seqLen,
        padLen,
        globalPos);
    totalTrimmed += sequences_[i].begTrim + sequences_[i].endTrim;

    // seqLen (from hashTableConfig_t) is the number of bases in the sequence including
    // trimmed N's.  The sequence stored in reference.bin does not include trimmed N's
    // and is used to calculate padding
    seqLenNoNs = sequences_[i].seqLen - sequences_[i].begTrim - sequences_[i].endTrim;

    // Global position includes padding - calculate the same padding calculated by the
    // hash table builder (from hash_table.c)
    padLen = REF_SEQ_MIN_PAD_BASES;
    padLen += (REF_SEQ_ALIGN_BASES - ((seqLenNoNs + padLen) % REF_SEQ_ALIGN_BASES)) % REF_SEQ_ALIGN_BASES;
    globalPos += sequences_[i].seqLen;  // Includes N's
  }

  return ret;
}
#endif

}  // namespace reference
}  // namespace dragenos
