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

#ifndef ALIGN_ALIGNMENT_HPP
#define ALIGN_ALIGNMENT_HPP

#include <boost/iterator/reverse_iterator.hpp>
#include <string>

#include "align/Cigar.hpp"
#include "align/Mapq.hpp"
#include "align/Query.hpp"
#include "map/SeedChain.hpp"

namespace dragenos {
namespace align {

typedef uint32_t FlagType;

struct AlignmentHeader {
  enum Flag : FlagType {
    NONE                            = 0,
    MULTIPLE_SEGMENTS               = 0x1,
    ALL_PROPERLY_ALIGNED            = 0x2,
    UNMAPPED                        = 0x4,
    UNMAPPD_NEXT_SEGMENT            = 0x8,
    REVERSE_COMPLEMENT              = 0x10,
    REVERSE_COMPLEMENT_NEXT_SEGMENT = 0x20,
    FIRST_IN_TEMPLATE               = 0x40,
    LAST_IN_TEMPLATE                = 0x80,
    SECONDARY_ALIGNMENT             = 0x100,
    FAILED_FILTERS                  = 0x200,
    DUPLICATE                       = 0x400,
    SUPPLEMENTARY_ALIGNMENT         = 0x800
  };

  AlignmentHeader(FlagType flags = -1, ScoreType score = -1) : score_(score), flags_(flags) {}

  ScoreType score_             = INVALID_SCORE;
  int       mismatches_        = -1;
  ScoreType potentialScore_    = INVALID_SCORE;
  ScoreType subScore_          = INVALID_SCORE;
  bool      perfect_           = false;
  MapqType  mapq_              = -1;
  FlagType  flags_             = -1;
  int       position_          = -1;
  short     reference_         = -1;
  int       nextPosition_      = -1;
  short     nextReference_     = -1;
  int       templateLength_    = -1;
  int32_t   mateCoordinate_    = -1;
  bool      smithWatermanDone_ = false;
  bool      filtered_          = false;

  void setUnmapped()
  {
    flags_ =
        UNMAPPED | (flags_ & (MULTIPLE_SEGMENTS | UNMAPPD_NEXT_SEGMENT | REVERSE_COMPLEMENT_NEXT_SEGMENT |
                              FIRST_IN_TEMPLATE | LAST_IN_TEMPLATE | FAILED_FILTERS));
  }

  void     setFlags(FlagType flags) { flags_ |= flags; }
  void     resetFlags(FlagType flags = 0) { flags_ = flags; }
  void     unsetFlags(FlagType flags = 0) { flags_ &= ~flags; }
  FlagType getFlags() const { return flags_; }
  bool     hasMultipleSegments() const { return flags_ & MULTIPLE_SEGMENTS; }
  bool     areAllProperlyAligned() const { return flags_ & ALL_PROPERLY_ALIGNED; }
  bool     isUnmapped() const { return flags_ & UNMAPPED; }
  bool     isUnmappedNextSegment() const { return flags_ & UNMAPPD_NEXT_SEGMENT; }
  bool     isReverseComplement() const { return flags_ & REVERSE_COMPLEMENT; }
  bool     isReverseComplementNextSegment() const { return flags_ & REVERSE_COMPLEMENT_NEXT_SEGMENT; }
  bool     isFirstInTemplate() const { return flags_ & FIRST_IN_TEMPLATE; }
  bool     isLastInTemplate() const { return flags_ & LAST_IN_TEMPLATE; }
  bool     isSecondaryAlignment() const { return flags_ & SECONDARY_ALIGNMENT; }
  bool     hasFailedFilters() const { return flags_ & FAILED_FILTERS; }
  bool     isDuplicate() const { return flags_ & DUPLICATE; }
  bool     isSupplementaryAlignment() const { return flags_ & SUPPLEMENTARY_ALIGNMENT; }
  int      getPosition() const { return position_; }
  void     setPosition(const int position) { position_ = position; }
  short    getReference() const { return reference_; }
  void     setReference(unsigned short reference) { reference_ = reference; }
  int      getNextPosition() const { return nextPosition_; }
  void     setNextPosition(const int nextPosition) { nextPosition_ = nextPosition; }
  short    getNextReference() const { return nextReference_; }
  void     setNextReference(unsigned short nextReference) { nextReference_ = nextReference; }
  int      getTemplateLength() const { return templateLength_; }
  void     setTemplateLength(const int templateLength) { templateLength_ = templateLength; }
  int32_t  getMateCoordinate() const { return mateCoordinate_; }
  void     setMateCoordinate(const int32_t mateCoordinate) { mateCoordinate_ = mateCoordinate; }

  void      setScore(ScoreType score) { score_ = score; }
  ScoreType getScore() const { return score_; }
  void      setXs(ScoreType score) { subScore_ = score; }
  ScoreType getXs() const { return subScore_; }
  void      setMismatchCount(int mismatches) { mismatches_ = mismatches; }
  int       getMismatchCount() const { return mismatches_; }
  void      setPotentialScore(ScoreType score) { potentialScore_ = score; }
  ScoreType getPotentialScore() const { return potentialScore_; }
  void      setPerfect(bool perfect) { perfect_ = perfect; }
  bool      isPerfect() const { return perfect_; }
  void      setMapq(MapqType mapq)
  {
    mapq_ = std::max(0, mapq);
    //    if (!mapq_)
    //    {
    //      unsetFlags(ALL_PROPERLY_ALIGNED);
    //      setFlags(UNMAPPED);
    //      cigar_.clear();
    //    }
  }
  MapqType getMapq() const { return mapq_; }
  void     setSmithWatermanDone(bool flag) { smithWatermanDone_ = flag; }
  bool     isSmithWatermanDone() const { return smithWatermanDone_; }
  void     setFiltered(bool filtered) { filtered_ = filtered; }
  bool     isFiltered() const { return filtered_; }

  friend std::ostream& operator<<(std::ostream& os, const AlignmentHeader& aln)
  {
    return os << "AlignmentHeader(" << aln.flags_ << "f," << aln.score_ << "s," << aln.potentialScore_
              << "ps," << aln.mapq_ << "mq," << aln.getReference() << ":" << aln.getPosition()
              << (aln.isReverseComplement() ? "r," : "f,") << "-" << aln.getNextReference() << ":"
              << aln.getNextPosition() << (aln.isReverseComplementNextSegment() ? "r," : "f,") << ")";
  }
};

/**
 ** \brief encapsulates all the alignment information
 **
 **/
class Alignment : public AlignmentHeader {
public:
  Alignment(FlagType flags = 0, ScoreType score = -1) : AlignmentHeader(flags, score) {}

  const Cigar& getCigar() const { return cigar_; }
  Cigar&       cigar() { return cigar_; }
  template <typename QueryIt, typename DbIt>
  uint32_t setCigarOperations(
      const std::string& operations,
      DbIt               dbBegin,
      DbIt               dbEnd,
      QueryIt            queryBegin,
      QueryIt            queryEnd,
      bool               reverse,
      int                softClipStart = 0)
  {
    auto ret = cigar_.setOperationSequence(operations, softClipStart);
    setTemplateLength(cigar_.getReferenceLength());
    if (reverse) {
      setMismatchCount(countEdits(
          operations,
          dbBegin,
          dbEnd,
          boost::make_reverse_iterator(queryEnd),
          boost::make_reverse_iterator(queryBegin)));
    } else {
      setMismatchCount(countEdits(operations, dbBegin, dbEnd, queryBegin, queryEnd));
    }
    return ret;
  }

  const map::SeedChain& chain() const
  {
    assert(chain_);
    return *chain_;
  }
  void setChain(const map::SeedChain& c)
  {
    assert(!chain_ || &c == chain_);
    chain_ = &c;
  }
  bool hasOnlyRandomSamples() const { return isUnmapped() || chain().hasOnlyRandomSamples(); }
  bool isExtra() const { return isUnmapped() || chain().isExtra(); }

  uint32_t setCigarOperations(const std::string& operations, int softClipStart = 0)
  {
    auto ret = cigar_.setOperationSequence(operations, softClipStart);
    setTemplateLength(cigar_.getReferenceLength());
    return ret;
  }

  const AlignmentHeader& header() const { return *this; }

  int64_t getUnclippedAlignmentCoordinate() const
  {
    const int64_t offset = isUnmapped()
                               ? 0
                               : (isReverseComplement()) ? (getCigar().getReferenceLengthPlusEndClips() - 1)
                                                         : -int64_t(getCigar().countStartClips());
    int64_t result = static_cast<int64_t>(getPosition()) + offset;
    return result;
  }

  int64_t getUnclippedStartPosition() const { return position_ - getCigar().countStartClips(); }

  int64_t getUnclippedEndPosition() const { return position_ + cigar_.getReferenceLengthPlusEndClips() - 1; }

  uint32_t getRefLen() const { return cigar_.getReferenceLength(); }

  int64_t getEndPosition() const { return position_ + cigar_.getReferenceLength() - 1; }

  void setIneligibilityStatus(bool status) { ineligible_ = status; }

  bool getIneligibilityStatus() const { return ineligible_; }

  bool isDuplicate(const Alignment& other) const
  {
    if (isUnmapped() == other.isUnmapped() && isReverseComplement() == other.isReverseComplement() &&
        getReference() == other.getReference() &&
        ((getCigar().countStartClips() == other.getCigar().countStartClips() &&
          getPosition() == other.getPosition()) ||
         (getCigar().countEndClips() == other.getCigar().countEndClips() &&
          getEndPosition() == other.getEndPosition()))) {
      //            std::cerr << "duplicate:" << " our:" << *this << " other:" << other << std::endl;
      return true;
    }
    //    std::cerr << "no duplicate:" << " our:" << *this << " other:" << other << std::endl;
    return false;
  }

  bool isOverlap(const Alignment& other) const
  {
    assert(!other.getCigar().empty());
    const int ourStartClips   = getCigar().countStartClips();
    const int theirStartClips = other.isReverseComplement() != isReverseComplement()
                                    ? other.getCigar().countEndClips()
                                    : other.getCigar().countStartClips();

    const int ourClippedLength   = getCigar().getClippedLength();
    const int theirClippedLength = other.getCigar().getClippedLength();
    const int overlap = std::min(ourStartClips + ourClippedLength, theirStartClips + theirClippedLength) -
                        std::max(ourStartClips, theirStartClips);

    if (overlap * 2 >= int((std::min(ourClippedLength, theirClippedLength)))) {
      //            std::cerr << "overlap:" << overlap << " our:" << *this << " other:" << other << std::endl;
      return true;
    } else {
      //            std::cerr << "no overlap:" << overlap << " our:" << *this << " other:" << other << std::endl;
    }

    return false;
  }

  bool reverse() const { return isReverseComplement(); }

  void             setSa(const Alignment* sa) { sa_ = sa; }
  const Alignment* getSa() const { return sa_; }

  short getNm() const { return getMismatchCount(); }

private:
  const map::SeedChain* chain_ = 0;
  Cigar                 cigar_;
  const Alignment*      sa_ = 0;
  /// lengths clipped at the beginning and at the end of the sequence
  //int clip_[2];
  //std::vector<unsigned char> optional_;
  bool                 ineligible_ = false;
  friend std::ostream& operator<<(std::ostream& os, const Alignment& aln)
  {
    return os << "Alignment(" << aln.flags_ << "f," << aln.score_ << "s," << aln.potentialScore_ << "ps,"
              << aln.mapq_ << "mq," << aln.getReference() << ":" << aln.getPosition()
              << (aln.isReverseComplement() ? "r" : "f") << (aln.isUnmapped() ? "u," : ",") << "-"
              << aln.getNextReference() << ":" << aln.getNextPosition()
              << (aln.isReverseComplementNextSegment() ? "r," : "f,") << aln.cigar_
              << (aln.perfect_ ? "perf," : "") << (aln.ineligible_ ? "inel" : "") << ")";
  }

  template <typename QueryIt, typename DbIt>
  static uint32_t countEdits(
      const std::string& operations, DbIt dbIt, DbIt dbEnd, QueryIt queryIt, QueryIt queryEnd)
  {
    std::string::const_iterator opIt  = operations.begin();
    int                         edits = 0;
    while (queryEnd != queryIt) {
      assert(operations.end() != opIt);
      if ('I' == *opIt) {
        ++edits;
        ++queryIt;
      } else if ('S' == *opIt) {
        ++queryIt;
      } else if ('D' == *opIt) {
        assert(dbEnd != dbIt);
        ++edits;
        ++dbIt;
      } else if ('N' == *opIt) {
        assert(dbEnd != dbIt);
        ++dbIt;
      } else  // if ('M' == *opIt)
      {
        assert('M' == *opIt);
        assert(dbEnd != dbIt);
        edits += *queryIt != *dbIt;
        ++queryIt;
        ++dbIt;
      }

      ++opIt;
    }

    return edits;
  }
};

struct SerializedSaTag {
  struct Header {
    short    reference_ = 0;
    int      position_  = 0;
    bool     reverse_   = 0;
    MapqType mapq_      = 0;
    short    nm_        = 0;
  } header_;

  SerializedCigar cigar_;

  void operator<<(const Alignment& alignment)
  {
    header_.reference_ = alignment.getReference();
    header_.position_  = alignment.getPosition();
    header_.reverse_   = alignment.isReverseComplement();
    header_.mapq_      = alignment.getMapq();
    header_.nm_        = alignment.getMismatchCount();
    cigar_ << alignment.getCigar();
  }

  static std::size_t getByteSize(const Alignment& alignment)
  {
    return sizeof(SerializedSaTag) + SerializedCigar::getByteSize(alignment.getCigar()) -
           sizeof(SerializedCigar);
  }

  std::size_t getByteSize() const { return sizeof(*this) + cigar_.getByteSize() - sizeof(cigar_); }

  short                  getReference() const { return header_.reference_; }
  int                    getPosition() const { return header_.position_; }
  bool                   reverse() const { return header_.reverse_; }
  MapqType               getMapq() const { return header_.mapq_; }
  short                  getNm() const { return header_.nm_; }
  const SerializedCigar& getCigar() const { return cigar_; }
};

class SerializedAlignment : public AlignmentHeader {
  static const std::size_t MAX_SA_TAG_LENGTH_ = 1024;
  bool                     hasSa_             = false;
  SerializedCigar          cigar_;

public:
  const SerializedCigar& getCigar() const { return cigar_; }

  void operator<<(const Alignment& alignment)
  {
    header() = alignment.header();
    cigar_ << alignment.getCigar();
    if (alignment.getSa()) {
      hasSa_ = true;
      sa() << *alignment.getSa();
    } else {
      hasSa_ = false;
    }
  }

  std::size_t getByteSize() const
  {
    return sizeof(SerializedAlignment) + cigar_.getByteSize() - sizeof(SerializedCigar) +
           (getSa() ? getSa()->getByteSize() : 0);
  }

  static std::size_t getByteSize(const Alignment& alignment)
  {
    return sizeof(SerializedAlignment) + SerializedCigar::getByteSize(alignment.getCigar()) -
           sizeof(SerializedCigar) +
           (alignment.getSa() ? SerializedSaTag::getByteSize(*alignment.getSa()) : 0);
  }

  AlignmentHeader&       header() { return *this; }
  const AlignmentHeader& header() const { return *this; }

  const SerializedSaTag* getSa() const
  {
    return hasSa_ ? (SerializedSaTag*)((const char*)&cigar_ + cigar_.getByteSize()) : 0;
  }

  SerializedSaTag& sa() { return *(SerializedSaTag*)((const char*)&cigar_ + cigar_.getByteSize()); }

  friend std::ostream& operator<<(std::ostream& os, const SerializedAlignment& a)
  {
    return os << "SerializedAlignment(" << a.header() << ")";
  }
};

class AlignmentPair {
  typedef std::array<Alignment*, 2>            Alignments;
  typedef std::array<const map::SeedChain*, 2> SeedChains;
  Alignments                                   alignments_;
  SeedChains                                   seedChains_;
  ScoreType                                    score_;
  ScoreType                                    potentialScore_;
  bool                                         properPair_;

public:
  AlignmentPair(Alignment& a1, Alignment& a2)
    : alignments_{&a1, &a2}, seedChains_{nullptr, nullptr}, score_(0), potentialScore_(0), properPair_(false)
  {
  }

  AlignmentPair(
      Alignment&            unmapped1,
      Alignment&            a2,
      const map::SeedChain* s2,
      ScoreType             score,
      ScoreType             potentialScore)
    : alignments_{&unmapped1, &a2},
      seedChains_{0, s2},
      score_(score),
      potentialScore_(potentialScore),
      properPair_(false)
  {
  }

  AlignmentPair(
      Alignment&            a1,
      const map::SeedChain* s1,
      Alignment&            unmapped2,
      ScoreType             score,
      ScoreType             potentialScore)
    : alignments_{&a1, &unmapped2},
      seedChains_{s1, 0},
      score_(score),
      potentialScore_(potentialScore),
      properPair_(false)
  {
  }

  Alignment&       at(std::size_t i) { return *alignments_.at(i); }
  const Alignment& cat(std::size_t i) const { return *alignments_.at(i); }
  const Alignment& at(std::size_t i) const { return cat(i); }
  Alignment&       operator[](std::size_t i) { return *alignments_.at(i); }
  const Alignment& operator[](std::size_t i) const { return *alignments_.at(i); }
  void             setSeedChains(const map::SeedChain* s0, const map::SeedChain* s1)
  {
    seedChains_[0] = s0;
    seedChains_[1] = s1;
  }
  const SeedChains& getSeedChains() const { return seedChains_; }
  ScoreType         getScore() const { return score_; }
  void              setScore(ScoreType score) { score_ = score; }
  ScoreType         getPotentialScore() const { return potentialScore_; }
  void              setPotentialScore(ScoreType score) { potentialScore_ = score; }
  void              setProperPair(bool properPair) { properPair_ = properPair; }
  bool              isProperPair() const { return properPair_; }
  bool              isPerfect() const
  {
    return (nullptr != seedChains_[0]) && (nullptr != seedChains_[1]) && seedChains_[0]->isPerfect() &&
           seedChains_[1]->isPerfect();
  }
  bool hasOnlyRandomSamples() const
  {
    return (at(0).isUnmapped() || at(0).hasOnlyRandomSamples()) &&
           (at(1).isUnmapped() || at(1).hasOnlyRandomSamples());
  }
  bool isExtra() const
  {
    return (at(0).isUnmapped() || at(0).isExtra()) && (at(1).isUnmapped() || at(1).isExtra());
  }
  bool                 isFiltered() const { return at(0).isFiltered() || at(1).isFiltered(); }
  friend std::ostream& operator<<(std::ostream& os, const AlignmentPair& p)
  {
    return os << "AlignmentPair(" << (p.isProperPair() ? "p," : "np,")
              << (p.isSingleEnded() ? " sep, " : " pep,") << p.score_ << "s," << p.potentialScore_ << "ps,"
              << p[0] << "," << p[1] << ")";
  }

  bool isSingleEnded() const { return at(0).isUnmapped() || at(1).isUnmapped(); }
};

typedef std::vector<AlignmentPair> AlignmentPairs;

}  // namespace align
}  // namespace dragenos

#endif  // #ifndef ALIGN_ALIGNMENT_HPP
