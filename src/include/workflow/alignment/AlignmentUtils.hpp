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

#pragma once

#include "align/Aligner.hpp"
#include "align/Pairs.hpp"
#include "align/Tlen.hpp"

namespace dragenos {
namespace workflow {
namespace alignment {

template <typename StoreOp>
void storeSeSupplementary(
    const align::SinglePicker&            singlePicker,
    const sequences::Read&                read,
    dragenos::align::Aligner::Alignments& alignments,
    dragenos::align::Alignment*           primary,
    StoreOp                               store)
{
  align::Alignments::iterator sup = singlePicker.findSupplementary(read.getLength(), alignments, primary);

  if (alignments.end() != sup) {
    sup->setFlags(align::Alignment::SUPPLEMENTARY_ALIGNMENT);
    sup->cigar().softClipsToHardClips();
    sup->setTemplateLength(0);
    primary->setSa(&*sup);
    sup->setSa(primary);
    store(read, *sup);
  }
}

template <typename StoreOp>
bool storeSeSecondary(
    const align::SinglePicker&            singlePicker,
    const sequences::Read&                read,
    dragenos::align::Aligner::Alignments& alignments,
    dragenos::align::Alignment*           primary,
    StoreOp                               store)
{
  return singlePicker.findSecondary(read.getLength(), alignments, primary, [&](align::Alignment& a) {
    a.setFlags(align::Alignment::SECONDARY_ALIGNMENT);
    a.cigar().softClipsToHardClips();
    a.setTemplateLength(0);
    a.setXs(primary->getScore());
    a.setMapq(0);
    store(read, a);
  });
}

template <typename StoreOp>
void storePeSupplementary(
    const align::InsertSizeParameters&    insertSizeParameters,
    const align::SinglePicker&            singlePicker,
    const sequences::Read&                read,
    dragenos::align::Aligner::Alignments& alignments,
    dragenos::align::Alignment*           primary,
    const dragenos::align::Alignment&     mate,
    StoreOp                               store)
{
  align::Alignments::iterator sup = singlePicker.findSupplementary(read.getLength(), alignments, primary);

  if (alignments.end() != sup) {
    sup->setFlags(
        align::Alignment::SUPPLEMENTARY_ALIGNMENT | align::Alignment::MULTIPLE_SEGMENTS |
        (mate.isReverseComplement() ? align::Alignment::REVERSE_COMPLEMENT_NEXT_SEGMENT
                                    : align::Alignment::NONE));
    sup->cigar().softClipsToHardClips();
    sup->setNextPosition(mate.getPosition());
    if (mate.getReference() != sup->getReference()) {
      sup->setTemplateLength(0);
      sup->setNextReference(mate.getReference());
    } else {
      sup->setTemplateLength(getTlen(*sup, mate, insertSizeParameters.orientation_));
    }
    primary->setSa(&*sup);
    sup->setSa(primary);
    sup->unsetFlags(align::Alignment::ALL_PROPERLY_ALIGNED);
    store(read, *sup);
  }
}

/**
 * \brief           Find all secondary alignments and store them properly
 * \precondition    PairBuilder::findSecondary must not fail due to sec-aligns-hard
 */
template <typename StoreOp>
void storePeSecondary(
    const align::InsertSizeParameters&    insertSizeParameters,
    const align::PairBuilder&             pairBuilder,
    const sequences::ReadPair&            pair,
    align::AlignmentPairs&                alignmentPairs,
    const align::Alignments&              unpaired,
    align::AlignmentPairs::const_iterator primary,
    int                                   readIdx,
    StoreOp                               store)
{
  //  std::cerr << "p:" << *primary << std::endl;
  bool ok = pairBuilder.findSecondary(
      pair.getLength(), alignmentPairs, unpaired, primary, readIdx, [&](align::AlignmentPair& ap) {
        align::Alignment& a = ap.at(readIdx);
        align::Alignment& m = ap.at(!readIdx);
        //      std::cerr << "a:" << a << std::endl;
        //      std::cerr << "m:" << m << std::endl;
        a.setFlags(
            align::Alignment::SECONDARY_ALIGNMENT | align::Alignment::MULTIPLE_SEGMENTS |
            (m.isUnmapped() ? align::Alignment::UNMAPPD_NEXT_SEGMENT : align::Alignment::NONE) |
            (m.isReverseComplement() ? align::Alignment::REVERSE_COMPLEMENT_NEXT_SEGMENT
                                     : align::Alignment::NONE) |
            (ap.isProperPair() ? align::Alignment::ALL_PROPERLY_ALIGNED : align::Alignment::NONE));

        a.cigar().softClipsToHardClips();
        a.setXs(primary->at(readIdx).getScore());
        a.setMapq(0);
        a.setNextPosition(m.isUnmapped() ? a.getPosition() : m.getPosition());
        if (!m.isUnmapped() && m.getReference() != a.getReference()) {
          a.setNextReference(m.getReference());
          a.setTemplateLength(0);
        } else {
          a.setTemplateLength(getTlen(a, m, insertSizeParameters.orientation_));
        }
        store(pair.at(readIdx), a);
      });
  // expected to work
  assert(ok);
}

template <typename StoreOp>
void storePairedBest(
    const align::InsertSizeParameters& insertSizeParameters,
    const sequences::Read&             read,
    align::Alignment&                  a,
    const align::AlignmentPair&        pair,
    StoreOp                            store)
{
  const align::Alignment& m = pair.at(a.isFirstInTemplate());
  a.setNextPosition(m.getPosition());
  a.setFlags(
      align::Alignment::MULTIPLE_SEGMENTS |
      (m.isUnmapped() ? align::Alignment::UNMAPPD_NEXT_SEGMENT : align::Alignment::NONE) |
      (m.isReverseComplement() ? align::Alignment::REVERSE_COMPLEMENT_NEXT_SEGMENT : align::Alignment::NONE) |
      (pair.isProperPair() ? align::Alignment::ALL_PROPERLY_ALIGNED : align::Alignment::NONE));

  if (a.getReference() == m.getReference()) {
    const int tlen = getTlen(a, m, insertSizeParameters.orientation_);

    if (!align::pairMatch(insertSizeParameters, a, m) || !pair.isProperPair()) {
      a.unsetFlags(align::Alignment::ALL_PROPERLY_ALIGNED);
    }

    a.setTemplateLength(tlen);
    a.setMateCoordinate(m.getUnclippedAlignmentCoordinate());
  } else {
    a.unsetFlags(align::Alignment::ALL_PROPERLY_ALIGNED);
    a.setNextReference(m.getReference());
    a.setTemplateLength(0);
  }

  store(read, a);
}

template <typename StoreOp>
void storeUnmappedPair(const align::Aligner::ReadPair& pair, StoreOp store)
{
  static const align::Alignment unmappedR1(
      align::AlignmentHeader::MULTIPLE_SEGMENTS | align::AlignmentHeader::UNMAPPED |
      align::AlignmentHeader::UNMAPPD_NEXT_SEGMENT | align::AlignmentHeader::FIRST_IN_TEMPLATE);

  static const align::Alignment unmappedR2(
      align::AlignmentHeader::MULTIPLE_SEGMENTS | align::AlignmentHeader::UNMAPPED |
      align::AlignmentHeader::UNMAPPD_NEXT_SEGMENT | align::AlignmentHeader::LAST_IN_TEMPLATE);

  store(pair.at(0), unmappedR1);
  store(pair.at(1), unmappedR2);
}

template <typename StoreOp>
void alignAndStorePair(
    const align::InsertSizeParameters& insertSizeParameters,
    const sequences::ReadPair&         pair,
    align::Aligner&                    aligner,
    const align::SinglePicker&         singlePicker,
    const align::PairBuilder&          pairBuilder,
    align::AlignmentPairs&             alignmentPairs,
    StoreOp                            store)
{
  alignmentPairs.clear();
  const auto best = aligner.getAlignments(pair, alignmentPairs, insertSizeParameters, pairBuilder);
  if (alignmentPairs.end() != best) {
    // if we go over sec-aligns when sec-aligns-hard is set, make sure we don't store anything.
    if (!pairBuilder.findSecondary(
            pair.getLength(),
            alignmentPairs,
            aligner.unpaired(0),
            best,
            0,
            [&](align::AlignmentPair& /*ap*/) {}) ||
        !pairBuilder.findSecondary(
            pair.getLength(),
            alignmentPairs,
            aligner.unpaired(1),
            best,
            1,
            [&](align::AlignmentPair& /*ap*/) {})) {
      storeUnmappedPair(pair, store);
      return;
    }

    if (!best->at(0).isUnmapped()) {
      storePeSecondary(
          insertSizeParameters, pairBuilder, pair, alignmentPairs, aligner.unpaired(0), best, 0, store);
      storePeSupplementary(
          insertSizeParameters,
          singlePicker,
          pair.at(0),
          aligner.unpaired(0),
          &best->at(0),
          best->at(1),
          store);
    }
    storePairedBest(insertSizeParameters, pair.at(0), best->at(0), *best, store);

    if (!best->at(1).isUnmapped()) {
      storePeSecondary(
          insertSizeParameters, pairBuilder, pair, alignmentPairs, aligner.unpaired(1), best, 1, store);
      storePeSupplementary(
          insertSizeParameters,
          singlePicker,
          pair.at(1),
          aligner.unpaired(1),
          &best->at(1),
          best->at(0),
          store);
    }
    storePairedBest(insertSizeParameters, pair.at(1), best->at(1), *best, store);
  } else {
    storeUnmappedPair(pair, store);
  }
}

template <typename StoreOp>
void storeSingleEnded(const sequences::Read& read, align::Alignment& a, StoreOp store)
{
  a.unsetFlags(align::Alignment::FIRST_IN_TEMPLATE);
  store(read, a);
}

template <typename StoreOp>
void alignAndStoreSingle(
    const align::Aligner::Read& read,
    align::Aligner&             aligner,
    const align::SinglePicker&  singlePicker,
    align::Aligner::Alignments& alignments,
    StoreOp                     store)
{
  static const align::Alignment unmappedSE(align::AlignmentHeader::UNMAPPED);

  alignments.clear();
  aligner.getAlignments(read, alignments);
  const auto best = singlePicker.pickBest(read.getLength(), alignments);
  if (alignments.end() != best) {
    if (!storeSeSecondary(
            singlePicker, read, alignments, &*best, [&](const sequences::Read& read, align::Alignment& a) {
              storeSingleEnded(read, a, store);
            })) {
      // sec-aligns-hard is set and we have too many secondaries.
      store(read, unmappedSE);
      return;
    }

    storeSeSupplementary(
        singlePicker, read, alignments, &*best, [&](const sequences::Read& read, align::Alignment& a) {
          storeSingleEnded(read, a, store);
        });

    // for validation maintain the order same way as in dragen. Supplementary alignment is printed first
    best->setTemplateLength(0);
    storeSingleEnded(read, *best, store);
  } else {
    store(read, unmappedSE);
  }
}

}  // namespace alignment
}  // namespace workflow
}  // namespace dragenos
