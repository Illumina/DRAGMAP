// Copyright 2016 Edico Genome Corporation. All rights reserved.
//
// This file contains confidential and proprietary information of the Edico
// Genome Corporation and is protected under the U.S. and international
// copyright and other intellectual property laws.
//
// $Id$
// $Author$
// $Change$
// $DateTime$
//

#ifndef __OUTPUT_DBAM_HEADER_HPP__
#define __OUTPUT_DBAM_HEADER_HPP__

#include "align/Alignment.hpp"
#include "align/Cigar.hpp"

namespace BamTools {
namespace Constants {

const uint8_t BAM_SIZEOF_INT = 4;

// header magic number
const char *const BAM_HEADER_MAGIC = "BAM\1";
const uint8_t BAM_HEADER_MAGIC_LENGTH = 4;

// BAM alignment core size
const uint8_t BAM_CORE_SIZE = 32;
const uint8_t BAM_CORE_BUFFER_SIZE = 8;

// BAM alignment flags
const int BAM_ALIGNMENT_PAIRED = 0x0001;
const int BAM_ALIGNMENT_PROPER_PAIR = 0x0002;
const int BAM_ALIGNMENT_UNMAPPED = 0x0004;
const int BAM_ALIGNMENT_MATE_UNMAPPED = 0x0008;
const int BAM_ALIGNMENT_REVERSE_STRAND = 0x0010;
const int BAM_ALIGNMENT_MATE_REVERSE_STRAND = 0x0020;
const int BAM_ALIGNMENT_READ_1 = 0x0040;
const int BAM_ALIGNMENT_READ_2 = 0x0080;
const int BAM_ALIGNMENT_SECONDARY = 0x0100;
const int BAM_ALIGNMENT_QC_FAILED = 0x0200;
const int BAM_ALIGNMENT_DUPLICATE = 0x0400;

} // namespace Constants
} // namespace BamTools

class DbamHeader {
  const dragenos::align::SerializedAlignment &alignment_;
  uint8_t peStatsInterval_ = 0;

  const dragenos::sequences::SerializedRead &read_;

public:
  enum { ALIGNMENT_FLAG_SUPPLEMENTARY = 0x800 };
  enum { SUPPRESS_OUTPUT = 0x1000 };

  DbamHeader(const dragenos::align::SerializedAlignment &alignment,
             const dragenos::sequences::SerializedRead &read)
      : alignment_(alignment), read_(read) {}

  int32_t getTemplateLen() const { return alignment_.getTemplateLength(); }
  int64_t getUnclippedAlignmentCoordinate() const {
    const int64_t offset =
        (alignment_.isReverseComplement())
            ? (alignment_.getCigar().getReferenceLengthPlusEndClips() - 1)
            : -int64_t(alignment_.getCigar().countStartClips());
    int64_t result = static_cast<int64_t>(alignment_.getPosition()) + offset;
    return result;
  }

  const dragenos::align::SerializedCigar &getCigar() const {
    return alignment_.getCigar();
  }

  int32_t getMateCoordinate() const { return alignment_.getMateCoordinate(); }

  // GR  what is called mismatch in Alignment class seems to be the edit
  // distance
  int getEditDistance() const { return alignment_.getMismatchCount(); }

  bool isSecondary() const { return alignment_.isSecondaryAlignment(); }

  bool isSupplementary() const { return alignment_.isSupplementaryAlignment(); }

  bool isPrimary() const { return (not(isSupplementary() or isSecondary())); }

  bool isFirstInPair() const {
    return alignment_.hasMultipleSegments() && alignment_.isFirstInTemplate();
  }

  uint8_t getMapQuality() const { return alignment_.getMapq(); }

  bool hasMate() const { return alignment_.hasMultipleSegments(); }

  bool isMateUnmapped() const { return alignment_.isUnmappedNextSegment(); }

  bool isProperlyPaired() const { return alignment_.areAllProperlyAligned(); }

  uint16_t getFlag() const { return alignment_.getFlags(); }

  bool suppressOutput() const {
    return (alignment_.getFlags() & SUPPRESS_OUTPUT);
  }

  bool isPairMapped() const {
    int flag = alignment_.getFlags();
    return (flag & dragenos::align::AlignmentHeader::MULTIPLE_SEGMENTS) &&
           (0 ==
            (flag & dragenos::align::AlignmentHeader::UNMAPPD_NEXT_SEGMENT)) &&
           (0 == (flag & dragenos::align::AlignmentHeader::UNMAPPED));
  }

  short getSequenceLen() const { return read_.getReadLen(); }

  bool isMateOnSameContig() const {
    return alignment_.getReference() == alignment_.getNextReference();
  }

  bool isDisqualified() const { return (alignment_.hasFailedFilters()); }

  bool isUnmapped() const { return (alignment_.isUnmapped()); }

  bool isDuplicate() const { return (alignment_.isDuplicate()); }

  const uint8_t *getConstQualities() const {
    return reinterpret_cast<const uint8_t *>(read_.getQualities().first);
  }

  uint8_t getPeStatsInterval() const { return peStatsInterval_; }

  void setPeStatsInterval(const uint8_t interval) {
    peStatsInterval_ = interval;
  }
};

#endif
