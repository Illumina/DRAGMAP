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

#ifndef ALIGN_SAM_HPP
#define ALIGN_SAM_HPP

#include <boost/range/adaptor/reversed.hpp>
#include <iostream>

#include "align/Mapq.hpp"
#include "reference/HashtableConfig.hpp"
#include "sequences/Read.hpp"

namespace dragenos {
namespace sam {

class SamGenerator {
  static const char                 Q0_ = 33;
  const reference::HashtableConfig& hashtableConfig_;

public:
  SamGenerator(const reference::HashtableConfig& hashtableConfig) : hashtableConfig_(hashtableConfig) {}

  template <typename ReadT>
  static std::string getReadName(const ReadT& read)
  {
    // make sure that only the actual name is used - skip everything after the first space
    const auto& fullName = read.getName();
    const auto  nameEnd  = std::find_if(std::begin(fullName), std::end(fullName), isspace);
    return std::string(fullName.begin(), nameEnd);
  }
  // generate record mapped as described by an alignment structure
  template <typename ReadT, typename AlignmenT>
  std::ostream& generateRecord(
      std::ostream& os, const ReadT& read, const AlignmenT& alignment, const std::string& rgid) const
  {
    os << getReadName(read) << '\t' << alignment.getFlags() << '\t';
    if (alignment.isUnmapped() && (!alignment.hasMultipleSegments() || alignment.isUnmappedNextSegment())) {
      os << "*\t0\t0\t*\t";
    } else {
      os << (-1 == alignment.getReference() ? std::string("=")
                                            : hashtableConfig_.getSequenceName(alignment.getReference()))
         << '\t' << alignment.getPosition() + 1 << '\t'
         << std::min<align::MapqType>(alignment.getMapq(), align::MAPQ_MAX) << '\t';
      if (alignment.getCigar().empty()) {
        os << "*\t";
      } else {
        os << alignment.getCigar() << '\t';
      }
    }

    if (!alignment.hasMultipleSegments() || (alignment.isUnmapped() && alignment.isUnmappedNextSegment())) {
      os << "*\t0\t";
    } else {
      os << (-1 == alignment.getNextReference() || alignment.getReference() == alignment.getNextReference()
                 ? std::string("=")
                 : hashtableConfig_.getSequenceName(alignment.getNextReference()))
         << '\t' << alignment.getNextPosition() + 1 << '\t';
    }
    os << (alignment.isUnmapped() ? 0 : alignment.getTemplateLength()) << '\t';
    generateSequence(os, read, alignment) << '\t';
    generateQualities(os, read, alignment) << '\t';
    os << "RG:Z:" << rgid;
    if (-1 != alignment.getScore()) {
      os << "\tAS:i:" << alignment.getScore();
    }
    if (align::INVALID_SCORE != alignment.getXs()) {
      os << "\tXS:i:" << alignment.getXs();
    }
    if (-1 != alignment.getMismatchCount()) {
      os << "\tNM:i:" << alignment.getMismatchCount();
    }
    if (align::MAPQ_MAX < alignment.getMapq()) {
      os << "\tXQ:i:" << std::min<align::MapqType>(alignment.getMapq(), align::HW_MAPQ_MAX);
    }

    if (alignment.getSa()) {
      const auto& sa = *alignment.getSa();
      os << "\tSA:Z:" << hashtableConfig_.getSequenceName(sa.getReference()) << ',' << (sa.getPosition() + 1)
         << ',' << (sa.reverse() ? "-," : "+,") << sa.getCigar()
         << ','
         //         << std::min<MapqType>(sa.getMapq(), MAPQ_MAX) <<
         << std::min<align::MapqType>(sa.getMapq(), align::HW_MAPQ_MAX) << ',' << sa.getNm() << ';';
    }
    return os;
  }
  // generate an unmapped record
  template <typename ReadT>
  static std::ostream& generateRecord(
      std::ostream& os, const ReadT& read, const std::string& rgid, unsigned flags = 0x4)
  {
    os << getReadName(read) << '\t' << flags << '\t' << '*' << '\t'  // reference name
       << '0' << '\t'                                                // reference position
       << 0 << '\t'                                                  // MAPQ unavailable
       << '*' << '\t'                                                // CIGAR
       << '*' << '\t'                                                // next reference name
       << '0' << '\t'                                                // next reference position
       << '0' << '\t';                                               // template length
    generateSequence(os, read, false) << '\t';
    generateQualities(os, read, false) << '\t';
    os << "RG:Z:" << rgid;
    return os;
  }

  template <typename ReadT, typename AlignmenT>
  static std::ostream& generateSequence(std::ostream& os, const ReadT& read, const AlignmenT& a)
  {
    const auto& bases = read.getBases();
    const auto& cigar = a.getCigar();
    if (cigar.countEndHardClips() + cigar.countStartHardClips() >= int(bases.size())) return os;
    if (a.isReverseComplement()) {
      auto range = boost::adaptors::reverse(bases);
      range.advance_begin(cigar.countStartHardClips());
      range.advance_end(-cigar.countEndHardClips());
      for (auto it = std::begin(range); std::end(range) != it; ++it) {
        os << ReadT::decodeRcBase(*it);
      }
    } else {
      auto range = boost::make_iterator_range(bases);
      range.advance_begin(cigar.countStartHardClips());
      range.advance_end(-cigar.countEndHardClips());
      for (auto it = std::begin(range); std::end(range) != it; ++it) {
        os << ReadT::decodeBase(*it);
      }
    }
    return os;
  }

  template <typename ReadT, typename AlignmenT>
  static std::ostream& generateQualities(std::ostream& os, const ReadT& read, const AlignmenT& a)
  {
    const auto& qualities = read.getQualities();
    const auto& cigar     = a.getCigar();
    if (cigar.countEndHardClips() + cigar.countStartHardClips() >= int(qualities.size())) return os;
    if (a.isReverseComplement()) {
      auto range = boost::adaptors::reverse(qualities);
      range.advance_begin(cigar.countStartHardClips());
      range.advance_end(-cigar.countEndHardClips());
      for (auto it = std::begin(range); std::end(range) != it; ++it) {
        os << static_cast<char>(*it + Q0_);
      }
    } else {
      auto range = boost::make_iterator_range(qualities);
      range.advance_begin(cigar.countStartHardClips());
      range.advance_end(-cigar.countEndHardClips());
      for (auto it = std::begin(range); std::end(range) != it; ++it) {
        os << static_cast<char>(*it + Q0_);
      }
    }
    return os;
  }

  static std::ostream& generateHeader(
      std::ostream&                     os,
      const reference::HashtableConfig& hashtableConfig,
      const std::string&                commandLine,
      const std::string&                rgid,
      const std::string                 rgsm)
  {
    os << "@HD\tVN:1.4\tSO:unsorted\n";
    os << "@PG\tID: DRAGEN-OS\tVN:" DRAGEN_OS_VERSION "\tCL:" << commandLine << "\n";
    os << "@RG\tID:" << rgid << "\tLB:LB0\tPL:PL0\tPU:PU0\tSM:" << rgsm << "\n";

    // sequences must be generated in the original order but they are internally sorted by increasing
    // positions
    auto                                         sequences = hashtableConfig.getSequences();
    typedef reference::HashtableConfig::Sequence Sequence;
    std::sort(sequences.begin(), sequences.end(), [](const Sequence& lhs, const Sequence& rhs) {
      return lhs.id_ < rhs.id_;
    });
    const auto& sequenceNames = hashtableConfig.getSequenceNames();
    for (std::size_t s = 0; s < sequences.size(); ++s) {
      const auto& sequence = sequences.at(s);
      assert(sequence.id_ == s);
      os << "@SQ\tSN:" << sequenceNames[s] << "\tLN:" << sequence.seqLen << "\n";
    }
    return os;
  }
};

}  // namespace sam
}  // namespace dragenos

#endif  // #ifndef ALIGN_SAM_HPP
