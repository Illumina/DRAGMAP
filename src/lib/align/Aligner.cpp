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

#include <cassert>
#include <cerrno>
#include <cstring>
#include <queue>

#include <fcntl.h>
#include <sys/mman.h>

#include "common/DragenLogger.hpp"
#include "common/Exceptions.hpp"

#include "align/Aligner.hpp"
#include "align/CalculateRefStartEnd.hpp"
#include "align/Mapq.hpp"
#include "align/PairBuilder.hpp"
#include "align/Pairs.hpp"
#include "align/Tlen.hpp"
#include "sequences/Seed.hpp"

namespace dragenos {
namespace align {

Aligner::Aligner(
    const reference::ReferenceDir& referenceDir,
    const reference::Hashtable&    hashtable,
    const bool                     mapOnly,
    const int                      swAll,
    const SimilarityScores&        similarity,
    const int                      gapInit,
    const int                      gapExtend,
    const int                      unclipScore,
    const int                      alnMinScore,
    const uint32_t                 aln_cfg_unpaired_pen)
  : referenceDir_(referenceDir),
    mapOnly_(mapOnly),
    swAll_(swAll),
    mapper_(&hashtable),
    similarity_(similarity),
    gapInit_(gapInit),
    gapExtend_(gapExtend),
    unclipScore_(unclipScore),
    alnMinScore_(alnMinScore),
    aln_cfg_unpaired_pen_(aln_cfg_unpaired_pen),
    smithWaterman_(similarity, gapInit, gapExtend, unclipScore),
    alignmentGenerator_(referenceDir_, smithWaterman_)
{
}

// void Aligner::buildAlignments(map::ChainBuilder& chainBuilder, const Read &read, Alignments &alignments)
//{
//  // TODO: encapsulate these local buffers into the Alignment to enable reuse - possibly into the Alignment instance
//  SmithWaterman::Database database;
//  // sort the seed chains by decreasing length
//  // TODO: implement the exact same sorting mechanism as in DRAGEN
//  // const auto compare = [] (const map::SeedChain &lhs, const map::SeedChain &rhs) -> bool {return rhs.size() < lhs.size();};
//  // chainBuilder.sort(compare);
//  alignments.clear();
//  if (mapOnly_)
//  {
//    generateDummyAlignments(read, chainBuilder, alignments);
//  }
//  else
//  {
//    generateAlignments(read, chainBuilder, alignments);
//    // TODO order the alignments and calculate the single end alignment score and MAPQ
//  }
//}

void Aligner::buildUngappedAlignments(
    map::ChainBuilder& chainBuilder, const Read& read, Alignments& alignments)
{
  // sort the seed chains by decreasing length
  // TODO: implement the exact same sorting mechanism as in DRAGEN
  // const auto compare = [] (const map::SeedChain &lhs, const map::SeedChain &rhs) -> bool {return rhs.size() < lhs.size();};
  // chainBuilder.sort(compare);
  alignments.clear();
  if (mapOnly_) {
    generateDummyAlignments(read, chainBuilder, alignments);
  } else {
    generateUngappedAlignments(read, chainBuilder, alignments);
    // TODO order the alignments and calculate the single end alignment score and MAPQ
  }
}

void Aligner::generateUngappedAlignments(
    const Read& read, map::ChainBuilder& chainBuilder, Alignments& alignments)
{
  alignments.clear();
  for (auto& seedChain : chainBuilder) {
    //    if(seedChain.isFiltered()) continue;
    auto& alignment = alignments.addAlignment();
    generateUngappedAlignment(read, seedChain, alignment);
  }
}

void removeDuplicates(Alignments& alignments)
{
  std::sort(alignments.begin(), alignments.end(), [](const Alignment& left, const Alignment& right) {
    return left.getUnclippedEndPosition() < right.getUnclippedEndPosition() ||
           (left.getUnclippedEndPosition() == right.getUnclippedEndPosition() &&
            (left.isReverseComplement() < right.isReverseComplement() ||
             (left.isReverseComplement() == right.isReverseComplement() &&
              left.getScore() > right.getScore())));
  });

  alignments.erase(
      std::unique(
          alignments.begin(),
          alignments.end(),
          [](const Alignment& left, const Alignment& right) {
            return left.getUnclippedEndPosition() == right.getUnclippedEndPosition() &&
                   left.isReverseComplement() == right.isReverseComplement();
          }),
      alignments.end());

  std::sort(alignments.begin(), alignments.end(), [](const Alignment& left, const Alignment& right) {
    return left.getUnclippedStartPosition() < right.getUnclippedStartPosition() ||
           (left.getUnclippedStartPosition() == right.getUnclippedStartPosition() &&
            (left.isReverseComplement() < right.isReverseComplement() ||
             (left.isReverseComplement() == right.isReverseComplement() &&
              left.getScore() > right.getScore())));
  });

  alignments.erase(
      std::unique(
          alignments.begin(),
          alignments.end(),
          [](const Alignment& left, const Alignment& right) {
            return left.getUnclippedStartPosition() == right.getUnclippedStartPosition() &&
                   left.isReverseComplement() == right.isReverseComplement();
          }),
      alignments.end());
}

void Aligner::runSmithWatermanAll(
    const Read& read, const map::ChainBuilder& chainBuilder, Alignments& alignments)
{
  for (int i = 0; i < chainBuilder.size(); ++i) {
    alignmentGenerator_.generateAlignment(0, read, chainBuilder.at(i), alignments.at(i));
  }

  ////  if (DEBUG_FILES)
  //  {
  //    removeDuplicates(alignments);
  //  }
}

void Aligner::generateUngappedAlignment(const Read& read, map::SeedChain& seedChain, Alignment& alignment)
{
  FlagType flags = !read.getPosition() ? Alignment::FIRST_IN_TEMPLATE : Alignment::LAST_IN_TEMPLATE;
  flags |= seedChain.isReverseComplement() ? Alignment::REVERSE_COMPLEMENT : 0;
  alignment.resetFlags(flags);
  const size_t referenceOffset = seedChain.firstReferencePosition();

  int move =
      initializeUngappedAlignmentScores(read, seedChain.isReverseComplement(), referenceOffset, alignment);

  if (alnMinScore_ > alignment.getScore()) {
    alignment.setUnmapped();
  }
  // TODO: calculate score
  alignment.setPerfect(seedChain.isPerfect());
  //  alignment.setPerfect(alignment.getCigar().getClippedLength() == read.getLength());

#ifdef TRACE_ALIGNMENTS
  std::cerr << "[ALIGNMENT]\t"
            << "ungapped" << read.getPosition() << "\t" << alignment << std::endl;
#endif
}

unsigned Aligner::calculatePotentialScore(
    const Read& read, const unsigned score, const unsigned softClipTo, const unsigned softClipFrom) const
{
  const int SEED_LENGTH = mapper_.getHashtable()->getPrimarySeedBases();
  assert(21 == SEED_LENGTH);
  int ret = score;
  if (softClipTo > 0) {
    // assume single-base gap at the point of clipping
    ret -= gapInit_;
    // assume one mismatch per one whole seed
    const int mismatchCount = (softClipTo) / SEED_LENGTH;
    const int matchCount    = softClipTo - mismatchCount;
    ret += (mismatchCount * similarity_.mismatch_);
    ret += (matchCount * similarity_.match_);
  }
  if (softClipFrom < read.getLength()) {
    ret -= gapInit_;
    const int mismatchCount = (read.getLength() - softClipFrom) / SEED_LENGTH;
    const int matchCount    = read.getLength() - softClipFrom - mismatchCount;
    ret += (mismatchCount * similarity_.mismatch_);
    ret += (matchCount * similarity_.match_);
  }
  return std::max(0, ret);
}

int Aligner::initializeUngappedAlignmentScores(
    const Read& read, const bool rcFlag, const size_t referenceOffset, Alignment& alignment)
{
  std::string                operations;
  static const char          ALIGNMENT_MATCH      = Cigar::getOperationName(Cigar::ALIGNMENT_MATCH);
  static const char          SOFT_CLIP            = Cigar::getOperationName(Cigar::SOFT_CLIP);
  static const int           SOFT_CLIP_ADJUSTMENT = -5;
  static const unsigned char N                    = 0xF;
  operations.clear();
  const auto& readBases      = rcFlag ? read.getRcBases() : read.getBases();
  int         alignmentScore = -SOFT_CLIP_ADJUSTMENT;
  int         maxScore       = 0;
  int         lastPosition   = 0;
  // best stretch with strictly positive scores
  int bestFirst = 0;
  int bestLast  = 0;
  int bestScore = 0;
  // current stretch with strictly positive scores
  int currentFirst = 0;
  int currentLast  = 0;
  int currentScore = 0;

  const reference::HashtableConfig& hashtableConfig = referenceDir_.getHashtableConfig();
  auto refCoords = hashtableConfig.convertToReferenceCoordinates(referenceOffset);
  const reference::HashtableConfig::Sequence& seq      = hashtableConfig.getSequences().at(refCoords.first);
  const auto                                  posRange = hashtableConfig.getPositionRange(seq);
  const int seqLeft = std::min(readBases.size(), posRange.second - referenceOffset);

  for (unsigned i = 0; i < seqLeft; ++i) {
    const auto          readBase          = readBases[i];
    const unsigned char referenceBase     = referenceDir_.getReferenceSequence().getBase(referenceOffset + i);
    const bool          goodReadBase      = (readBase > 0) && (readBase != N);
    const bool          goodReferenceBase = (referenceBase > 0) && (referenceBase != N);
    const bool          goodBase          = goodReadBase && goodReferenceBase;
    const bool          match             = goodBase && (readBase == referenceBase);
    alignmentScore += similarity_(readBase, referenceBase);
    alignmentScore = std::max(0, alignmentScore);
    if (0 == alignmentScore) {
      if (currentScore > bestScore) {
        bestScore = currentScore;
        bestFirst = currentFirst;
        bestLast  = currentLast;
      }
      currentScore = 0;
      currentFirst = i;
      currentLast  = i;
    } else {
      if (0 == currentScore) {
        currentFirst = i;
      }

      if (alignmentScore > currentScore) {
        currentScore = alignmentScore;
        currentLast  = i;
      }
    }
  }
  if (currentScore > bestScore) {
    bestScore = currentScore;
    bestFirst = currentFirst;
    bestLast  = currentLast;
  }
  assert(bestScore > 0);
  assert(bestFirst <= bestLast);
  operations.clear();
  // check if initial soft clip is needed
  if (bestFirst > 0) {
    operations.insert(operations.size(), bestFirst, SOFT_CLIP);
  } else {
    // else, remove the initial score adjustment that was meant to give handycap to unclipped alignment
    bestScore += SOFT_CLIP_ADJUSTMENT;
  }
  // TODO: check if final soft clip is needed
  int malus = 0;
  for (unsigned i = bestLast + 1; i < seqLeft; ++i) {
    const auto          readBase      = readBases[i];
    const unsigned char referenceBase = referenceDir_.getReferenceSequence().getBase(referenceOffset + i);
    malus += similarity_(readBase, referenceBase);
  }
  if (bestLast + 1 < seqLeft && malus >= SOFT_CLIP_ADJUSTMENT) {
    bestLast = seqLeft - 1;
    bestScore += malus;
  }
  operations.insert(operations.size(), bestLast + 1 - bestFirst, ALIGNMENT_MATCH);
  if (read.getLength() > operations.size()) {
    const auto n = read.getLength() - operations.size();
    operations.insert(operations.size(), n, SOFT_CLIP);
  }
  assert(operations.size() == read.getLength());
  alignment.setScore(std::max(0, bestScore));
  alignment.setPotentialScore(calculatePotentialScore(read, bestScore, bestFirst, bestLast + 1));
  Database database;
  if (alignment.isReverseComplement()) {
    referenceDir_.getReferenceSequence().getRcBases(
        referenceOffset + bestFirst, referenceOffset + read.getLength(), std::back_inserter(database));
  } else {
    referenceDir_.getReferenceSequence().getBases(
        referenceOffset + bestFirst, referenceOffset + read.getLength(), std::back_inserter(database));
  }
  // deal with before reference starts
  const auto referenceCoordinates =
      referenceDir_.getHashtableConfig().convertToReferenceCoordinates(referenceOffset);
  const int refClip = std::max(-referenceCoordinates.second, int64_t(0));
  const int move    = alignment.setCigarOperations(
      operations,
      database.begin(),
      database.end(),
      read.getBases().begin(),
      read.getBases().end(),
      alignment.isReverseComplement(),
      std::max(refClip, bestFirst));

  alignment.setReference(referenceCoordinates.first);
  alignment.setPosition(referenceCoordinates.second + move);
  return bestFirst;
}

void Aligner::getAlignments(const Read& read, Alignments& alignments)
{
  alignments.clear();
  map::ChainBuilder& chainBuilder = chainBuilders_[0];
  chainBuilder.clear();
  mapper_.getPositionChains(read, chainBuilder);

  if (0 != chainBuilder.size()) {
    int bestScore = 0;
    buildUngappedAlignments(chainBuilder, read, alignments);
    if (swAll_) {
      runSmithWatermanAll(read, chainBuilder, alignments);
    } else {
      assert(chainBuilder.size() == alignments.size());
      const std::size_t toTry = alignments.size();
      for (unsigned i = 0; toTry > i; ++i) {
        auto&       alignment = alignments.at(i);
        const auto& seedChain = chainBuilder.at(i);
        if ((!alignment.isPerfect()) || (alignment.getPotentialScore() >= bestScore)) {
          if (alignment.isPerfect() && alignment.getScore() > bestScore) {
            bestScore = alignment.getScore();
          }

          if (alignment.getPotentialScore() > alignment.getScore()) {
            alignmentGenerator_.generateAlignment(0, read, seedChain, alignment);
          }
        }
      }
    }
  }
}

bool Aligner::rescueMate(
    const InsertSizeParameters& insertSizeParameters,
    const Read&                 anchoredRead,
    const Read&                 rescuedRead,
    const map::SeedChain&       anchoredSeedChain,
    Alignment&                  anchoredAlignment,
    const AlignmentRescue&      alignmentRescue,
    map::SeedChain&             rescuedSeedChain,
    Alignment&                  rescuedAlignment)
{
  if (alignmentRescue.scan(
          anchoredRead,
          rescuedRead,
          anchoredSeedChain,
          anchoredAlignment,
          referenceDir_.getReferenceSequence(),
          rescuedSeedChain)) {
    //          std::cerr << "rescued:" << rescuedSeedChain << std::endl;
    if (alignmentGenerator_.generateAlignment(
            alnMinScore_, rescuedRead, rescuedSeedChain, rescuedAlignment)) {
      //          bestUnpairedScore[rescuedReadPosition] = std::max(bestUnpairedScore[rescuedReadPosition], rescuedAlignment.getScore());
      // do Smith-Waterman on the anchored alignment if necessary
      if (!anchoredAlignment.isSmithWatermanDone() &&
          ((anchoredAlignment.getPotentialScore() > anchoredAlignment.getScore()) ||
           (!anchoredAlignment.isPerfect()))) {
        alignmentGenerator_.generateAlignment(
            alnMinScore_, anchoredRead, anchoredSeedChain, anchoredAlignment);
      }
#ifdef TRACE_ALIGNMENTS
      std::cerr << "[ALIGNMENT]\t"
                << "rescued"
                << "\t" << rescuedAlignment << std::endl;
#endif
      return true;
    }
  }

  return false;
}

void Aligner::makePair(
    const InsertSizeParameters& insertSizeParameters,
    const ReadPair&             readPair,
    const PairBuilder&          pairBuilder,
    const map::SeedChain*       s0,
    const map::SeedChain*       s1,
    Alignment&                  a0,
    Alignment&                  a1,
    AlignmentPairs&             alignmentPairs)
{
  alignmentPairs.insert(alignmentPairs.end(), AlignmentPair(a0, a1));
  const bool properPair =
      s0 && s1 ? pairMatch(insertSizeParameters, alignmentPairs.back()[0], alignmentPairs.back()[1]) : false;
  //  const bool properPair = s0 && s1 ? pairMatch(insertSizeParameters, readPair, *s0, *s1) : false;
  const int pair_pen = pairBuilder.computePairPenalty(insertSizeParameters, readPair, s0, s1, properPair);
  alignmentPairs.back().setScore(
      alignmentPairs.back()[0].getScore() + alignmentPairs.back()[1].getScore() - pair_pen);
  alignmentPairs.back().setPotentialScore(
      alignmentPairs.back()[0].getPotentialScore() + alignmentPairs.back()[1].getPotentialScore() - pair_pen);
  alignmentPairs.back().setSeedChains(s0, s1);
  alignmentPairs.back().setProperPair(properPair);
}

int findBest(const Alignments& alignments)
{
  const auto maxElement = std::max_element(
      alignments.begin(), alignments.end(), [](const Alignment& left, const Alignment& right) {
        return left.getScore() < right.getScore();
      });

  if (alignments.end() == maxElement || maxElement->isUnmapped()) {
    return -1;
  }
  return std::distance(alignments.begin(), maxElement);
}

bool Aligner::rescuePair(
    const InsertSizeParameters& insertSizeParameters,
    const PairBuilder&          pairBuilder,
    const ReadPair&             readPair,
    const int                   anchoredIdx,
    const bool                  any_pair_match,
    const map::SeedChain&       anchoredSeedChain,
    Alignment&                  anchoredAlignment,
    const AlignmentRescue       alignmentRescue,
    AlignmentPairs&             alignmentPairs)
{
  const int   rescuedIdx   = !anchoredIdx;
  const Read& anchoredRead = readPair.at(anchoredIdx);
  if (alignmentRescue.triggeredBy(anchoredRead, anchoredSeedChain, any_pair_match)) {
    const Read&    rescuedRead = readPair.at(rescuedIdx);
    map::SeedChain rescuedSeedChain;

    Alignment rescued;
    if (rescueMate(
            insertSizeParameters,
            anchoredRead,
            rescuedRead,
            anchoredSeedChain,
            anchoredAlignment,
            alignmentRescue,
            rescuedSeedChain,
            rescued)) {
      chainBuilders_[rescuedIdx].addSeedChain(rescuedSeedChain);
      unpairedAlignments_[rescuedIdx].append(rescued);
      makePair(
          insertSizeParameters,
          readPair,
          pairBuilder,
          anchoredIdx ? &chainBuilders_[rescuedIdx].back() : &anchoredSeedChain,
          anchoredIdx ? &anchoredSeedChain : &chainBuilders_[rescuedIdx].back(),
          anchoredIdx ? unpairedAlignments_[rescuedIdx].back() : anchoredAlignment,
          anchoredIdx ? anchoredAlignment : unpairedAlignments_[rescuedIdx].back(),
          alignmentPairs);
      return true;
    }
  }

  return false;
}

AlignmentPairs::iterator Aligner::getAlignments(
    const ReadPair&             readPair,
    AlignmentPairs&             alignmentPairs,
    const InsertSizeParameters& insertSizeParameters,
    const PairBuilder&          pairBuilder)
{
  alignmentPairs.clear();
  chainBuilders_[0].clear();
  chainBuilders_[1].clear();
  mapper_.getPositionChains(readPair[0], chainBuilders_[0]);
  mapper_.getPositionChains(readPair[1], chainBuilders_[1]);
  //  std::vector<std::array<map::ChainBuilder *, 2> > seedChainPairs; // keeping trace of the seed chains used for each
  // max number of chains is seed chains + rescued chains. Rescued is at most one per mate seed chain
  chainBuilders_[0].reserve(chainBuilders_[0].size() + chainBuilders_[1].size());
  chainBuilders_[1].reserve(chainBuilders_[0].size() + chainBuilders_[1].size());

  unpairedAlignments_[0].clear();
  unpairedAlignments_[1].clear();
  // max number is seed alignments + rescued alignments + 1. Rescued is at most one per mate seed chain. 1 for
  // unmapped mate
  unpairedAlignments_[0].reserve(chainBuilders_[0].size() + chainBuilders_[1].size() + 1);
  unpairedAlignments_[1].reserve(chainBuilders_[0].size() + chainBuilders_[1].size() + 1);

  buildUngappedAlignments(chainBuilders_[0], readPair[0], unpairedAlignments_[0]);
  buildUngappedAlignments(chainBuilders_[1], readPair[1], unpairedAlignments_[1]);

  if (swAll_) {
    runSmithWatermanAll(readPair[0], chainBuilders_[0], unpairedAlignments_[0]);
    runSmithWatermanAll(readPair[1], chainBuilders_[1], unpairedAlignments_[1]);
  }

  const auto bestOffset0 = findBest(unpairedAlignments_[0]);
  const auto bestOffset1 = findBest(unpairedAlignments_[1]);

  const AlignmentRescue alignmentRescue(
      insertSizeParameters.min_, insertSizeParameters.max_, insertSizeParameters.orientation_);

  ScoreType bestPairedScore = 0;

  // flag keeping track of the existence of any proper pair amongst the initial seed chains
  bool any_pair_match = false;
  // bitvectors recording for each seed chain for each end if a proper pair was found
  std::array<std::vector<bool>, 2> pairsFound;
  pairsFound[0].resize(unpairedAlignments_[0].size(), false);
  pairsFound[1].resize(unpairedAlignments_[1].size(), false);

  // find all combinations from seeds
  for (unsigned i0 = 0; pairsFound[0].size() > i0; ++i0) {
    for (unsigned i1 = 0; pairsFound[1].size() > i1; ++i1) {
      const bool hasBestUnpaired = ((bestOffset0 == i0) || (bestOffset1 == i1));
      const bool isPair =
          pairMatch(insertSizeParameters, readPair, chainBuilders_[0].at(i0), chainBuilders_[1].at(i1));
      if (isPair) {
        pairsFound[0][i0] = true;
        pairsFound[1][i1] = true;
        any_pair_match    = true;
      }
      if (hasBestUnpaired || isPair) {
        makePair(
            insertSizeParameters,
            readPair,
            pairBuilder,
            &chainBuilders_[0].at(i0),
            &chainBuilders_[1].at(i1),
            unpairedAlignments_[0].at(i0),
            unpairedAlignments_[1].at(i1),
            alignmentPairs);
        bestPairedScore = std::max(bestPairedScore, alignmentPairs.back().getScore());

#ifdef TRACE_ALIGNMENTS
        std::cerr << "[ALIGNMENT]\tcombined\t" << alignmentPairs.back() << std::endl;
#endif
      }
    }
  }

  // rescue for those that did not make any seed pairs
  for (unsigned i0 = 0; pairsFound[0].size() > i0; ++i0) {
    if (!pairsFound[0][i0]) {
      const auto& anchoredSeedChain = chainBuilders_[0].at(i0);
      auto&       anchoredAlignment = unpairedAlignments_[0].at(i0);

      if (rescuePair(
              insertSizeParameters,
              pairBuilder,
              readPair,
              0,
              any_pair_match,
              anchoredSeedChain,
              anchoredAlignment,
              alignmentRescue,
              alignmentPairs)) {
        pairsFound[0][i0] = true;
        bestPairedScore   = std::max(bestPairedScore, alignmentPairs.back().getScore());
#ifdef TRACE_ALIGNMENTS
        std::cerr << "[ALIGNMENT]\trescued0\t" << alignmentPairs.back() << std::endl;
#endif
      }
    }
  }

  for (unsigned i1 = 0; pairsFound[1].size() > i1; ++i1) {
    if (!pairsFound[1][i1]) {
      const auto& anchoredSeedChain = chainBuilders_[1].at(i1);
      auto&       anchoredAlignment = unpairedAlignments_[1].at(i1);
      if (rescuePair(
              insertSizeParameters,
              pairBuilder,
              readPair,
              1,
              any_pair_match,
              anchoredSeedChain,
              anchoredAlignment,
              alignmentRescue,
              alignmentPairs)) {
        pairsFound[1][i1] = true;
        bestPairedScore   = std::max(bestPairedScore, alignmentPairs.back().getScore());
#ifdef TRACE_ALIGNMENTS
        std::cerr << "[ALIGNMENT]\trescued1\t" << alignmentPairs.back() << std::endl;
#endif
      }
    }
  }

  // look into new combinations that became available
  for (unsigned i0 = 0; pairsFound[0].size() > i0; ++i0) {
    if (!pairsFound[0][i0]) {
      for (unsigned i1 = pairsFound[1].size(); unpairedAlignments_[1].size() > i1; ++i1) {
        makePair(
            insertSizeParameters,
            readPair,
            pairBuilder,
            0,
            0,
            unpairedAlignments_[0].at(i0),
            unpairedAlignments_[1].at(i1),
            alignmentPairs);
        bestPairedScore = std::max(bestPairedScore, alignmentPairs.back().getScore());
#ifdef TRACE_ALIGNMENTS
        std::cerr << "[ALIGNMENT]\tremains0\t" << alignmentPairs.back() << std::endl;
#endif
      }
    }
  }

  for (unsigned i1 = 0; pairsFound[1].size() > i1; ++i1) {
    if (!pairsFound[1][i1]) {
      for (unsigned i0 = pairsFound[0].size(); unpairedAlignments_[0].size() > i0; ++i0) {
        makePair(
            insertSizeParameters,
            readPair,
            pairBuilder,
            0,
            0,
            unpairedAlignments_[0].at(i0),
            unpairedAlignments_[1].at(i1),
            alignmentPairs);
        bestPairedScore = std::max(bestPairedScore, alignmentPairs.back().getScore());
#ifdef TRACE_ALIGNMENTS
        std::cerr << "[ALIGNMENT]\tremains1\t" << alignmentPairs.back() << std::endl;
#endif
      }
    }
  }

  // make single-ended pairs
  if (-1 != bestOffset0) {
    Alignment&       best0 = unpairedAlignments_[0].at(bestOffset0);
    align::Alignment unmappedR2(
        align::AlignmentHeader::MULTIPLE_SEGMENTS | align::AlignmentHeader::UNMAPPED |
        align::AlignmentHeader::LAST_IN_TEMPLATE);
    unmappedR2.setReference(best0.getReference());
    unmappedR2.setNextReference(best0.getReference());
    unmappedR2.setPosition(best0.getPosition());
    unmappedR2.setNextPosition(best0.getPosition());
    unpairedAlignments_[1].append(unmappedR2);

    const double m2a_scale = mapq2aln(similarity_.getSnpCost(), readPair.getLength());
    const int    pair_pen  = aln_cfg_unpaired_pen_ * m2a_scale;
    alignmentPairs.push_back(AlignmentPair(
        best0,
        chainBuilders_[0].at(bestOffset0),
        unpairedAlignments_[1].back(),
        best0.getScore() + alnMinScore_ - pair_pen,
        best0.getPotentialScore() + alnMinScore_ - pair_pen));
    bestPairedScore = std::max(bestPairedScore, alignmentPairs.back().getScore());
#ifdef TRACE_ALIGNMENTS
    std::cerr << "[ALIGNMENT]\tsep0\t" << alignmentPairs.back() << std::endl;
#endif
  }

  if (-1 != bestOffset1) {
    Alignment& best1 = unpairedAlignments_[1].at(bestOffset1);

    align::Alignment unmappedR1(
        align::AlignmentHeader::MULTIPLE_SEGMENTS | align::AlignmentHeader::UNMAPPED |
        align::AlignmentHeader::FIRST_IN_TEMPLATE);
    unmappedR1.setReference(best1.getReference());
    unmappedR1.setNextReference(best1.getReference());
    unmappedR1.setPosition(best1.getPosition());
    unmappedR1.setNextPosition(best1.getPosition());
    unpairedAlignments_[0].append(unmappedR1);

    const double m2a_scale = mapq2aln(similarity_.getSnpCost(), readPair.getLength());
    const int    pair_pen  = aln_cfg_unpaired_pen_ * m2a_scale;
    alignmentPairs.push_back(AlignmentPair(
        unpairedAlignments_[0].back(),
        best1,
        chainBuilders_[1].at(bestOffset1),
        best1.getScore() + alnMinScore_ - pair_pen,
        best1.getPotentialScore() + alnMinScore_ - pair_pen));
    bestPairedScore = std::max(bestPairedScore, alignmentPairs.back().getScore());
#ifdef TRACE_ALIGNMENTS
    std::cerr << "[ALIGNMENT]\tsep1\t" << alignmentPairs.back() << std::endl;
#endif
  }

  // when running without sw-all, apply smith-waterman to all that matters
  for (auto& alignmentPair : alignmentPairs) {
    const auto seedChains = alignmentPair.getSeedChains();
    if ((nullptr != seedChains[0]) && (nullptr != seedChains[1]) &&
        ((!alignmentPair.isPerfect()) || (alignmentPair.getPotentialScore() + 20 >= bestPairedScore))) {
      for (unsigned i = 0; 2 > i; ++i) {
        if ((nullptr != seedChains[i]) && !alignmentPair[i].isSmithWatermanDone() &&
            ((alignmentPair[i].getPotentialScore() + 20 > alignmentPair[i].getScore()) ||
             (!alignmentPair[i].isPerfect()))) {
          alignmentGenerator_.generateAlignment(alnMinScore_, readPair[i], *seedChains[i], alignmentPair[i]);
          // bestUnpairedScore[i] = std::max(bestUnpairedScore[i], alignmentPair[i].getScore());
        }
      }
      // const bool isPair   = pairMatch(insertSizeParameters, readPair, *seedChains[0], *seedChains[1]);
      const bool properPair = pairMatch(insertSizeParameters, alignmentPair[0], alignmentPair[1]);

      const int pair_pen = pairBuilder.computePairPenalty(
          insertSizeParameters, readPair, seedChains[0], seedChains[1], properPair);
      alignmentPair.setScore(alignmentPair[0].getScore() + alignmentPair[1].getScore() - pair_pen);
      alignmentPair.setPotentialScore(alignmentPair.getScore());
      alignmentPair.setProperPair(properPair);
      bestPairedScore = std::max(bestPairedScore, alignmentPair.getScore());
    }
  }

  if (!unpairedAlignments_[0].empty() && !unpairedAlignments_[1].empty()) {
    //    filter(alignmentPairs, unpairedAlignments_);
    return pairBuilder.pickBest(readPair, alignmentPairs, unpairedAlignments_);
  }

  return alignmentPairs.end();
}

void Aligner::generateDummyAlignments(
    const Read& read, const map::ChainBuilder& chainBuilder, Alignments& alignments) const
{
  alignments.clear();
  // TODO: implement
  for (const auto& seedChain : chainBuilder) {
    FlagType flags = !read.getPosition() ? Alignment::FIRST_IN_TEMPLATE : Alignment::LAST_IN_TEMPLATE;
    flags |= seedChain.isReverseComplement() ? Alignment::REVERSE_COMPLEMENT : 0;
    alignments.resize(alignments.size() + 1);
    auto& alignment = alignments.back();
    alignment.resetFlags(flags);
    // get the reference name from the ReferenceDir and the first diagonal
    const size_t                      referenceOffset = seedChain.firstReferencePosition();
    const reference::HashtableConfig& hashtableConfig = referenceDir_.getHashtableConfig();
    const auto referenceCoordinates = hashtableConfig.convertToReferenceCoordinates(referenceOffset);
    alignment.setReference(referenceCoordinates.first);
    const unsigned templateLength =
        seedChain.lastReferencePosition() - seedChain.firstReferencePosition() + 1;
    alignment.setTemplateLength(templateLength);
    auto& cigar = alignment.cigar();
    cigar.clear();
    const int beginClipLength = seedChain.firstReadBase();
    const int endClipLength   = read.getLength() - seedChain.lastReadBase() - 1;
    alignment.setPosition(referenceCoordinates.second + beginClipLength);
    const int totalClipLength = beginClipLength + endClipLength;
    alignment.setScore(read.getLength() - totalClipLength);
    if (beginClipLength > 0) {
      cigar.emplace_back(Cigar::SOFT_CLIP, beginClipLength);
    }
    if (read.getLength() == templateLength)  // straignt match
    {
      cigar.emplace_back(Cigar::ALIGNMENT_MATCH, templateLength - totalClipLength);
    } else if (read.getLength() < templateLength)  // insertion to the reference
    {
      const int totalMatchLength = read.getLength() - totalClipLength;
      const int firstMatchLength = totalMatchLength / 2;
      const int lastMatchLength  = totalMatchLength - firstMatchLength;
      cigar.emplace_back(Cigar::ALIGNMENT_MATCH, firstMatchLength);
      cigar.emplace_back(Cigar::DELETE, templateLength - read.getLength());
      cigar.emplace_back(Cigar::ALIGNMENT_MATCH, lastMatchLength);
    } else  // (read.getLength() > templateLength) deletion from the reference
    {
      const int totalMatchLength = templateLength - totalClipLength;
      const int firstMatchLength = totalMatchLength / 2;
      const int lastMatchLength  = totalMatchLength - firstMatchLength;
      cigar.emplace_back(Cigar::ALIGNMENT_MATCH, firstMatchLength);
      cigar.emplace_back(Cigar::INSERT, read.getLength() - templateLength);
      cigar.emplace_back(Cigar::ALIGNMENT_MATCH, lastMatchLength);
    }
    if (endClipLength > 0) {
      cigar.emplace_back(Cigar::SOFT_CLIP, endClipLength);
    }
  }
}

bool Aligner::isPerfectAlignment(
    const char* queryBegin,
    const char* queryEnd,
    const char* databaseBegin,
    const char* databaseEnd,
    Alignment&  alignment)
{
  assert(nullptr != queryBegin);
  assert(nullptr != queryEnd);
  assert(nullptr != databaseBegin);
  assert(nullptr != databaseEnd);
  assert(queryBegin <= queryEnd);
  assert(databaseBegin <= databaseEnd);
  // query and database must be exactly the same length and non-empty
  const auto count = queryEnd - queryBegin;
  if ((0 >= count) || (databaseEnd - databaseBegin != count)) {
    return false;
  }
  const char*                    query          = queryBegin;
  const char*                    database       = databaseBegin;
  int                            score          = 0;
  unsigned                       snpCount       = 0;  // number of SNPs in the current burst window
  constexpr int                  BURST_WINDOW   = 8;
  constexpr int                  BURST_MINIMUM  = 4;
  constexpr int                  MATCH_SCORE    = 1;
  constexpr int                  MATCH_N_SCORE  = -1;
  constexpr int                  MISMATCH_SCORE = -4;
  std::array<bool, BURST_WINDOW> burst;
  std::fill(burst.begin(), burst.end(), false);
  size_t front = 0;
  while (query != queryEnd) {
    // remove the SNP falling out of the window if any
    snpCount -= burst[front];
    // Ns have their own score, count as mismatches but not as burst
    constexpr char N = 0xF;
    if ((0 == *query) || (0 == *database) || (N == *query) || (N == *database)) {
      score += MATCH_N_SCORE;
      ++snpCount;
      burst[front] = false;
    }
    // TODO: add support for 2-3 base codes
    else if (*query == *database) {
      score += MATCH_SCORE;
      burst[front] = false;
    } else {
      score += MISMATCH_SCORE;
      ++snpCount;
      burst[front] = true;
    }
    if (BURST_MINIMUM <= snpCount) {
      return false;
    }
    front = ((front + 1) % burst.size());
    ++query;
    ++database;
  }
  if (score <= 0) {
    return false;
  }
  // We have a score, no burst - fill in the alignment
  // TODO: check for clipping
  // TODO: implement indel detection
  alignment.setScore(score);
  auto& cigar = alignment.cigar();
  cigar.clear();
  cigar.emplace_back(Cigar::ALIGNMENT_MATCH, count);
  return true;
}

void Aligner::filter(AlignmentPairs& alignmentPairs, std::array<Alignments, 2>& unpairedAlignments)
{
  for (unsigned i : {0, 1}) {
    filter(unpairedAlignments[i]);
  }
  for (auto& pair : alignmentPairs) {
    if (pair.isFiltered()) {
      pair.setScore(0);
    }
  }
}

void Aligner::filter(Alignments& alignments)
{
  // TODO: optimize (quadratic complexity)
  for (unsigned int i = 0; alignments.size() > i; ++i) {
    auto& a = alignments.at(i);
    if (a.isFiltered()) {
      continue;
    }
    for (unsigned int j = i + 1; alignments.size() > j; ++j) {
      auto& b = alignments.at(j);
      if (b.isFiltered()) {
        continue;
      }
      // TODO: add filtering for matching end positions
      if ((a.getReference() == b.getReference()) && (a.getPosition() == b.getPosition())) {
        if (a.getScore() < b.getScore()) {
          a.setFiltered(true);
          a.setScore(0);
          break;
        } else {
          b.setFiltered(true);
          b.setScore(0);
        }
      }
    }
  }
}

#if 0

static void checkDirectoryAndFile(const boost::filesystem::path &dir, const boost::filesystem::path &file)
{
  using namespace dragenos::common;
  if (!exists(dir))
    BOOST_THROW_EXCEPTION(IoException(ENOENT, std::string("ERROR: directory ") + dir.string() + " doesn't exist"));
  const auto filePath = dir / file;
  if (!exists(filePath))
    BOOST_THROW_EXCEPTION(IoException(ENOENT, std::string("ERROR: file not found: ") + filePath.string()));
  if (!is_regular_file(filePath))
    BOOST_THROW_EXCEPTION(IoException(ENOENT, std::string("ERROR: not a file: ") + filePath.string()));
}

static uintmax_t getFileSize(const boost::filesystem::path &filePath)
{
  using namespace dragenos::common;
  boost::system::error_code ec;
  const auto fileSize = file_size(filePath, ec);
  if (ec)
    BOOST_THROW_EXCEPTION(IoException(ec.value(), std::string("ERROR: failed to get stats for file: ") + filePath.string()));
  return fileSize;
}

std::vector<char> Aligner::getHashtableConfigData(const boost::filesystem::path referenceDir) const
{
  using namespace dragenos::common;
  checkDirectoryAndFile(referenceDir, "hash_table.cfg.bin");
  const auto hashtableConfigFile = referenceDir / "hash_table.cfg.bin";
  const auto fileSize = getFileSize(hashtableConfigFile);
  std::vector<char> data(fileSize);
  std::ifstream is(hashtableConfigFile.string());
  if ( (!is.read(data.data(), fileSize)) || (static_cast<long>(fileSize) != is.gcount()) )
  {
    boost::format message = boost::format("ERROR: failed to read %i bytes from binary hashtable config: %i bytes read from %s") % fileSize % is.gcount() % hashtableConfigFile.string();
    BOOST_THROW_EXCEPTION(IoException(errno, message.str()));
  }
  return data;
}

uint64_t *Aligner::mmapHashtableData(const boost::filesystem::path referenceDir) const
{
  using namespace dragenos::common;
  checkDirectoryAndFile(referenceDir, "hash_table.bin");
  const auto hashtableDataFile = referenceDir / "hash_table.bin";
  const auto fileSize = getFileSize(hashtableDataFile);
  if (fileSize != hashtableConfig_.getHashtableBytes())
  {
    boost::format message = boost::format("ERROR: hashtable size different from size in config file: expected %i bytes: actual %i bytes") % hashtableConfig_.getHashtableBytes() % fileSize;
    BOOST_THROW_EXCEPTION(IoException(errno, message.str()));
  }
  const auto hashtableFd = open(hashtableDataFile.c_str(), O_RDONLY, 0);
  if (-1 == hashtableFd)
  {
    BOOST_THROW_EXCEPTION(IoException(errno, std::string("ERROR: failed to open hashtable data file ") + hashtableDataFile.string()));
  }
  const int prot = PROT_READ;
  const int flags = MAP_PRIVATE | MAP_NORESERVE;
  const int offset = 0;
  auto table = static_cast<uint64_t*> (mmap(NULL, fileSize, prot, flags, hashtableFd, offset));
  if (MAP_FAILED == table)
  {
     BOOST_THROW_EXCEPTION(IoException(errno, std::string("ERROR: failed to map hashtable data file ") + hashtableDataFile.string()));
  }
  return reinterpret_cast<uint64_t *>(table);
}

#endif

}  // namespace align
}  // namespace dragenos
