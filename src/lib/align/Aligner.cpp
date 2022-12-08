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
    const reference::ReferenceSequence& refSeq,
    const reference::HashtableConfig&   htConfig,
    const reference::Hashtable&         hashtable,
    const bool                          mapOnly,
    const int                           swAll,
    const SimilarityScores&             similarity,
    const int                           gapInit,
    const int                           gapExtend,
    const int                           unclipScore,
    const int                           alnMinScore,
    const int                           aln_cfg_mapq_min_len,
    const uint32_t                      aln_cfg_unpaired_pen,
    const double                        aln_cfg_filter_len_ratio,
    const bool                          vectorizedSW)
  : refSeq_(refSeq),
    htConfig_(htConfig),
    mapOnly_(mapOnly),
    swAll_(swAll),
    vectorizedSW_(vectorizedSW),
    mapper_(&hashtable),
    similarity_(similarity),
    gapInit_(gapInit),
    gapExtend_(gapExtend),
    unclipScore_(unclipScore),
    alnMinScore_(alnMinScore),
    aln_cfg_mapq_min_len_(aln_cfg_mapq_min_len),
    aln_cfg_unpaired_pen_(aln_cfg_unpaired_pen),
    smithWaterman_(similarity, gapInit, gapExtend, unclipScore),
    vectorSmithWaterman_(similarity, gapInit, gapExtend, unclipScore),
    alignmentGenerator_(refSeq_, htConfig_, smithWaterman_, vectorSmithWaterman_, vectorizedSW_),
    chainBuilders_{map::ChainBuilder(aln_cfg_filter_len_ratio), map::ChainBuilder(aln_cfg_filter_len_ratio)}
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

void Aligner::runSmithWatermanAll(
    const Read& read, const map::ChainBuilder& chainBuilder, Alignments& alignments, const int readIdx)
{
  for (std::size_t i = 0; i < chainBuilder.size(); ++i) {
    alignmentGenerator_.generateAlignment(alnMinScore_, read, chainBuilder.at(i), alignments.at(i), readIdx);
  }

  ////  if (DEBUG_FILES)
  //  {
  //    removeDuplicates(alignments);
  //  }
}

void Aligner::updateIneligibility(
    const reference::HashtableConfig& hashtableConfig,
    const size_t                      referenceOffset,
    const Read&                       read,
    Alignment&                        alignment)
{
  auto        referenceCoordinates = hashtableConfig.convertToReferenceCoordinates(referenceOffset);
  const auto& sequences            = hashtableConfig.getSequences();
  const auto& sequence             = sequences.at(referenceCoordinates.first);
  if (referenceCoordinates.second + read.getLength() <= 0 || sequence.seqLen <= referenceCoordinates.second) {
    alignment.setIneligibilityStatus(true);
  }
}

void Aligner::generateUngappedAlignment(const Read& read, map::SeedChain& seedChain, Alignment& alignment)
{
  alignment.setChain(seedChain);
  FlagType flags = !read.getPosition() ? Alignment::FIRST_IN_TEMPLATE : Alignment::LAST_IN_TEMPLATE;
  flags |= seedChain.isReverseComplement() ? Alignment::REVERSE_COMPLEMENT : Alignment::NONE;
  flags |= seedChain.isFiltered() ? Alignment::UNMAPPED : Alignment::NONE;
  alignment.resetFlags(flags);

  const size_t referenceOffset = seedChain.firstReferencePosition();
  updateIneligibility(htConfig_, referenceOffset, read, alignment);

  if (seedChain.isFiltered()) {
    return;
  }

  //  int move =
  initializeUngappedAlignmentScores(read, seedChain.isReverseComplement(), referenceOffset, alignment);

  if (alnMinScore_ > alignment.getScore()) {
    alignment.setUnmapped();
  } else {
    alignment.setPerfect(seedChain.isPerfect());
  }

#ifdef TRACE_ALIGNMENTS
  std::cerr << "[ALIGNMENT]\t"
            << "ungapped" << read.getPosition() << "\t" << alignment << std::endl;
#endif
}

void Aligner::deFilterChain(
    const Read& read, map::SeedChain& seedChain, Alignment& alignment, const int readIdx)
{
  if (seedChain.isFiltered()) {
    seedChain.setFiltered(false);
    generateUngappedAlignment(read, seedChain, alignment);
    if (swAll_) {
      alignmentGenerator_.generateAlignment(alnMinScore_, read, seedChain, alignment, readIdx);
    }
  }
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
  std::string       operations;
  static const char ALIGNMENT_MATCH      = Cigar::getOperationName(Cigar::ALIGNMENT_MATCH);
  static const char SOFT_CLIP            = Cigar::getOperationName(Cigar::SOFT_CLIP);
  static const int  SOFT_CLIP_ADJUSTMENT = -5;
  //  static const unsigned char N                    = 0xF;
  operations.clear();
  const auto& readBases      = rcFlag ? read.getRcBases() : read.getBases();
  int         alignmentScore = -SOFT_CLIP_ADJUSTMENT;
  //  int         maxScore       = 0;
  //  int         lastPosition   = 0;
  // best stretch with strictly positive scores
  int bestFirst = 0;
  int bestLast  = 0;
  int bestScore = 0;
  // current stretch with strictly positive scores
  int currentFirst = 0;
  int currentLast  = 0;
  int currentScore = 0;

  auto refCoords                                  = htConfig_.convertToReferenceCoordinates(referenceOffset);
  const reference::HashtableConfig::Sequence& seq = htConfig_.getSequences().at(refCoords.first);
  const auto                                  posRange = htConfig_.getPositionRange(seq);
  const int seqLeft = std::min(readBases.size(), posRange.second - referenceOffset);

  Database databaseSeqLeft;
  databaseSeqLeft.reserve(seqLeft + 1);
  refSeq_.getBases(referenceOffset, referenceOffset + seqLeft, databaseSeqLeft);

  for (int i = 0; i < seqLeft; ++i) {
    const auto          readBase      = readBases[i];
    const unsigned char referenceBase = databaseSeqLeft[i];
    //    const bool          goodReadBase      = (readBase > 0) && (readBase != N);
    //    const bool          goodReferenceBase = (referenceBase > 0) && (referenceBase != N);
    //    const bool          goodBase          = goodReadBase && goodReferenceBase;
    //    const bool          match             = goodBase && (readBase == referenceBase);
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

  Database databaseBestLastToSeqLeft;
  size_t   bestToLastRefStart = referenceOffset + bestLast + 1;
  databaseBestLastToSeqLeft.reserve(seqLeft + 1);
  refSeq_.getBases(bestToLastRefStart, bestToLastRefStart + seqLeft, databaseBestLastToSeqLeft);

  /*
  databaseBestLastToSeqLeft.reserve(seqLeft - bestLast + 1);
  refSeq_.getBases(referenceOffset + bestLast + 1, referenceOffset + bestLast + 1 + seqLeft, databaseBestLastToSeqLeft);
  */

  for (int i = bestLast + 1; i < seqLeft; ++i) {
    const auto          readBase      = readBases[i];
    const unsigned char referenceBase = databaseBestLastToSeqLeft[i - (bestLast + 1)];
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
  database.reserve(read.getLength() - bestFirst + 1);
  if (alignment.isReverseComplement()) {
    refSeq_.getRcBases(referenceOffset + bestFirst, referenceOffset + read.getLength(), database);
  } else {
    refSeq_.getBases(referenceOffset + bestFirst, referenceOffset + read.getLength(), database);
  }
  // deal with before reference starts
  const auto referenceCoordinates = htConfig_.convertToReferenceCoordinates(referenceOffset);
  const int  refClip              = std::max(-referenceCoordinates.second, int64_t(0));
  const int  move                 = alignment.setCigarOperations(
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

void Aligner::runSmithWatermanWorthy(
    const Read& read, map::ChainBuilder& chainBuilder, Alignments& alignments)
{
  int               bestScore = 0;
  const std::size_t toTry     = alignments.size();
  for (unsigned i = 0; toTry > i; ++i) {
    auto&       alignment = alignments.at(i);
    const auto& seedChain = chainBuilder.at(i);
    if ((!alignment.isPerfect()) || (alignment.getPotentialScore() >= bestScore)) {
      if (alignment.isPerfect() && alignment.getScore() > bestScore) {
        bestScore = alignment.getScore();
      }
      if (alignment.getPotentialScore() > alignment.getScore()) {
        alignmentGenerator_.generateAlignment(alnMinScore_, read, seedChain, alignment, 0);
      }
    }
  }
}

void Aligner::getAlignments(const Read& read, Alignments& alignments)
{
  if (vectorizedSW_) {
    const auto& query0 = read.getBases();
    vectorSmithWaterman_.initReadContext(query0.data(), query0.data() + query0.size(), 0);
  }

  alignments.clear();
  map::ChainBuilder& chainBuilder = chainBuilders_[0];
  chainBuilder.clear();
  mapper_.getPositionChains(read, chainBuilder);

  if (0 != chainBuilder.size()) {
    buildUngappedAlignments(chainBuilder, read, alignments);
    if (swAll_) {
      runSmithWatermanAll(read, chainBuilder, alignments, 0);
    } else {
      runSmithWatermanWorthy(read, chainBuilder, alignments);
    }
  }

  if (vectorizedSW_) {
    vectorSmithWaterman_.destroyReadContext(0);
  }
}

bool Aligner::rescueMate(
    const InsertSizeParameters& /*insertSizeParameters*/,
    const Read&            anchoredRead,
    const Read&            rescuedRead,
    const map::SeedChain&  anchoredSeedChain,
    Alignment&             anchoredAlignment,
    const AlignmentRescue& alignmentRescue,
    map::SeedChain&        rescuedSeedChain,
    Alignment&             rescuedAlignment,
    const int              anchoredIdx)
{
  const int rescuedIdx = !anchoredIdx;

  if (alignmentRescue.scan(rescuedRead, anchoredSeedChain, refSeq_, rescuedSeedChain)) {
    //          std::cerr << "rescued:" << rescuedSeedChain << std::endl;
    if (alignmentGenerator_.generateAlignment(
            alnMinScore_, rescuedRead, rescuedSeedChain, rescuedAlignment, rescuedIdx)) {
      //          bestUnpairedScore[rescuedReadPosition] = std::max(bestUnpairedScore[rescuedReadPosition], rescuedAlignment.getScore());
      // do Smith-Waterman on the anchored alignment if necessary
      if (!anchoredAlignment.isSmithWatermanDone() &&
          ((anchoredAlignment.getPotentialScore() > anchoredAlignment.getScore()) ||
           (!anchoredAlignment.isPerfect()))) {
        alignmentGenerator_.generateAlignment(
            alnMinScore_, anchoredRead, anchoredSeedChain, anchoredAlignment, anchoredIdx);
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
  AlignmentPair& pair       = alignmentPairs.back();
  const bool     properPair = s0 && s1 ? pairMatch(insertSizeParameters, pair[0], pair[1]) : false;
  //  const bool properPair = s0 && s1 ? pairMatch(insertSizeParameters, readPair, *s0, *s1) : false;
  int       insert_len  = 0;
  int       insert_diff = 0;
  const int pair_pen    = pairBuilder.computePairPenalty(
      insertSizeParameters, readPair, &a0, &a1, properPair, insert_len, insert_diff);
  pair.setScore(pair[0].getScore() + pair[1].getScore() - pair_pen);
  pair.setPotentialScore(pair[0].getPotentialScore() + pair[1].getPotentialScore() - pair_pen);
  pair.setSeedChains(s0, s1);
  pair.setProperPair(properPair);

  DRAGEN_ALN_RESULT_LOG << "PAIR_RESULT:"
                        << " type=" << (!properPair ? "se_pair" : "paired") << ", paired=" << properPair
                        << ", insert_len=" << insert_len << ", insert_diff=" << insert_diff
                        << ", pair_pen=" << pair_pen << ", pair_score=" << pair.getScore() << std::endl;
}

int findBest(const Alignments& alignments)
{
  const auto maxElement = std::max_element(
      alignments.begin(), alignments.end(), [](const Alignment& left, const Alignment& right) {
        return left.getIneligibilityStatus() < right.getIneligibilityStatus() ||
               (left.getIneligibilityStatus() == right.getIneligibilityStatus() &&
                left.getScore() < right.getScore());
      });

  if (alignments.end() == maxElement || maxElement->isUnmapped() || maxElement->getIneligibilityStatus()) {
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
  if (alignmentRescue.triggeredBy(anchoredSeedChain, any_pair_match)) {
    map::SeedChain rescuedSeedChain;
    const int      rescuedIdx   = !anchoredIdx;
    const Read&    rescuedRead  = readPair.at(rescuedIdx);
    const Read&    anchoredRead = readPair.at(anchoredIdx);

    Alignment rescued;
    if (rescueMate(
            insertSizeParameters,
            anchoredRead,
            rescuedRead,
            anchoredSeedChain,
            anchoredAlignment,
            alignmentRescue,
            rescuedSeedChain,
            rescued,
            anchoredIdx)) {
      chainBuilders_[rescuedIdx].addSeedChain(rescuedSeedChain);
      rescued.setChain(chainBuilders_[rescuedIdx].back());

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
  if (vectorizedSW_) {
    const auto& query0 = readPair[0].getBases();
    vectorSmithWaterman_.initReadContext(query0.data(), query0.data() + query0.size(), 0);

    const auto& query1 = readPair[1].getBases();
    vectorSmithWaterman_.initReadContext(query1.data(), query1.data() + query1.size(), 1);
  }

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
    runSmithWatermanAll(readPair[0], chainBuilders_[0], unpairedAlignments_[0], 0);
    runSmithWatermanAll(readPair[1], chainBuilders_[1], unpairedAlignments_[1], 1);
  }

  int bestOffset0 = findBest(unpairedAlignments_[0]);
  int bestOffset1 = findBest(unpairedAlignments_[1]);
  assert(-1 == bestOffset0 || !chainBuilders_[0].at(bestOffset0).isFiltered());
  assert(-1 == bestOffset1 || !chainBuilders_[1].at(bestOffset1).isFiltered());

  const AlignmentRescue alignmentRescue(
      insertSizeParameters.rescueMin_, insertSizeParameters.rescueMax_, insertSizeParameters.orientation_);

  ScoreType bestPairedScore = 0;

  // flag keeping track of the existence of any proper pair amongst the initial seed chains
  bool any_pair_match = false;
  // bitvectors recording for each seed chain for each end if a proper pair was found
  std::array<std::vector<bool>, 2> pairsFound;
  pairsFound[0].resize(unpairedAlignments_[0].size(), false);
  pairsFound[1].resize(unpairedAlignments_[1].size(), false);

  //  std::cerr << "pairsFound[0].size():" << pairsFound[0].size() << std::endl;
  //  std::cerr << "pairsFound[1].size():" << pairsFound[1].size() << std::endl;
  // find all combinations from seeds
  for (unsigned i0 = 0; pairsFound[0].size() > i0; ++i0) {
    auto& chain0 = chainBuilders_[0].at(i0);
    if (unpairedAlignments_[0].at(i0).getIneligibilityStatus()) {
      continue;
    }

    for (unsigned i1 = 0; pairsFound[1].size() > i1; ++i1) {
      if (unpairedAlignments_[1].at(i1).getIneligibilityStatus()) {
        continue;
      }
      auto&      chain1          = chainBuilders_[1].at(i1);
      const bool hasBestUnpaired = !chain0.isFiltered() && !chain1.isFiltered() &&
                                   ((std::size_t(bestOffset0) == i0) || (std::size_t(bestOffset1) == i1));
      const bool isPair = pairMatch(insertSizeParameters, readPair, chain0, chain1);
      if (isPair) {
        pairsFound[0].at(i0) = true;
        pairsFound[1].at(i1) = true;
        any_pair_match |= !chain0.hasOnlyRandomSamples() && !chain1.hasOnlyRandomSamples();
        deFilterChain(readPair.at(0), chain0, unpairedAlignments_[0].at(i0), 0);
        deFilterChain(readPair.at(1), chain1, unpairedAlignments_[1].at(i1), 1);
      }
      if (hasBestUnpaired || isPair) {
        makePair(
            insertSizeParameters,
            readPair,
            pairBuilder,
            &chain0,
            &chain1,
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

  if (insertSizeParameters.isInitDone()) {
    // rescue for those that did not make any seed pairs
    for (unsigned i0 = 0; pairsFound[0].size() > i0; ++i0) {
      if (unpairedAlignments_[0].at(i0).getIneligibilityStatus()) {
        continue;
      }
      if (!pairsFound[0][i0]) {
        const auto& anchoredSeedChain = chainBuilders_[0].at(i0);
        auto&       anchoredAlignment = unpairedAlignments_[0].at(i0);
        if (anchoredAlignment.getIneligibilityStatus()) {
          continue;
        }

        if (!anchoredSeedChain.isFiltered() && rescuePair(
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
      if (unpairedAlignments_[1].at(i1).getIneligibilityStatus()) {
        continue;
      }
      if (!pairsFound[1][i1]) {
        const auto& anchoredSeedChain = chainBuilders_[1].at(i1);
        auto&       anchoredAlignment = unpairedAlignments_[1].at(i1);
        if (anchoredAlignment.getIneligibilityStatus()) {
          continue;
        }

        if (!anchoredSeedChain.isFiltered() && rescuePair(
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
  }

  // look into new combinations that became available
  if (-1 != bestOffset0) {
    const int i0 = bestOffset0;
    for (unsigned i1 = pairsFound[1].size(); unpairedAlignments_[1].size() > i1; ++i1) {
      if (unpairedAlignments_[1].at(i1).getIneligibilityStatus()) {
        continue;
      }
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

  if (-1 != bestOffset1) {
    const int i1 = bestOffset1;
    for (unsigned i0 = pairsFound[0].size(); unpairedAlignments_[0].size() > i0; ++i0) {
      if (unpairedAlignments_[0].at(i0).getIneligibilityStatus()) {
        continue;
      }
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

  // in case rescues have added something better
  bestOffset0 = findBest(unpairedAlignments_[0]);
  // make single-ended pairs
  if (-1 != bestOffset0) {
    Alignment& best0 = unpairedAlignments_[0].at(bestOffset0);
    if (!best0.getIneligibilityStatus()) {
      align::Alignment unmappedR1(
          align::AlignmentHeader::MULTIPLE_SEGMENTS | align::AlignmentHeader::UNMAPPED |
              align::AlignmentHeader::LAST_IN_TEMPLATE,
          ((-1 != bestOffset1 && unpairedAlignments_[1].at(bestOffset1).getIneligibilityStatus())
               ? std::max(alnMinScore_, unpairedAlignments_[1].at(bestOffset1).getScore())
               : alnMinScore_));
      unmappedR1.setReference(best0.getReference());
      unmappedR1.setNextReference(best0.getReference());
      unmappedR1.setPosition(best0.getPosition());
      unmappedR1.setNextPosition(best0.getPosition());
      unpairedAlignments_[1].append(unmappedR1);

      const int m2a_scale =
          mapq2aln(similarity_.getSnpCost(), std::max(aln_cfg_mapq_min_len_, readPair.getLength()));
      const int pair_pen = (aln_cfg_unpaired_pen_ * m2a_scale) >> 10;
      alignmentPairs.push_back(AlignmentPair(
          best0,
          chainBuilders_[0].size() > std::size_t(bestOffset0) ? &chainBuilders_[0].at(bestOffset0) : 0,
          unpairedAlignments_[1].back(),
          best0.getScore() + unmappedR1.getScore() - pair_pen,
          best0.getPotentialScore() + unmappedR1.getScore() - pair_pen));
      bestPairedScore = std::max(bestPairedScore, alignmentPairs.back().getScore());
#ifdef TRACE_ALIGNMENTS
      std::cerr << "[ALIGNMENT]\tsep0\t" << alignmentPairs.back() << std::endl;
#endif
      const AlignmentPair& pair = alignmentPairs.back();
      DRAGEN_ALN_RESULT_LOG << "PAIR_RESULT:"
                            << " type="
                            << "se_pair"
                            << ", paired=" << 0 << ", pair_pen=" << pair_pen
                            << ", pair_score=" << pair.getScore() << std::endl;
    }
  }

  // in case rescues have added something better
  bestOffset1 = findBest(unpairedAlignments_[1]);
  if (-1 != bestOffset1) {
    Alignment& best1 = unpairedAlignments_[1].at(bestOffset1);
    if (!best1.getIneligibilityStatus()) {
      align::Alignment unmappedR0(
          align::AlignmentHeader::MULTIPLE_SEGMENTS | align::AlignmentHeader::UNMAPPED |
              align::AlignmentHeader::FIRST_IN_TEMPLATE,
          ((-1 != bestOffset0 && unpairedAlignments_[0].at(bestOffset0).getIneligibilityStatus())
               ? std::max(alnMinScore_, unpairedAlignments_[0].at(bestOffset0).getScore())
               : alnMinScore_));
      unmappedR0.setReference(best1.getReference());
      unmappedR0.setNextReference(best1.getReference());
      unmappedR0.setPosition(best1.getPosition());
      unmappedR0.setNextPosition(best1.getPosition());
      unpairedAlignments_[0].append(unmappedR0);

      const int m2a_scale =
          mapq2aln(similarity_.getSnpCost(), std::max(aln_cfg_mapq_min_len_, readPair.getLength()));
      const int pair_pen = (aln_cfg_unpaired_pen_ * m2a_scale) >> 10;
      alignmentPairs.push_back(AlignmentPair(
          unpairedAlignments_[0].back(),
          best1,
          chainBuilders_[1].size() > std::size_t(bestOffset1) ? &chainBuilders_[1].at(bestOffset1) : 0,
          best1.getScore() + unmappedR0.getScore() - pair_pen,
          best1.getPotentialScore() + unmappedR0.getScore() - pair_pen));
      bestPairedScore = std::max(bestPairedScore, alignmentPairs.back().getScore());
#ifdef TRACE_ALIGNMENTS
      std::cerr << "[ALIGNMENT]\tsep1\t" << alignmentPairs.back() << std::endl;
#endif
      const AlignmentPair& pair = alignmentPairs.back();
      DRAGEN_ALN_RESULT_LOG << "PAIR_RESULT:"
                            << " type="
                            << "se_pair"
                            << ", paired=" << 0 << ", pair_pen=" << pair_pen
                            << ", pair_score=" << pair.getScore() << std::endl;
    }
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
          alignmentGenerator_.generateAlignment(
              alnMinScore_, readPair[i], *seedChains[i], alignmentPair[i], i);
          // bestUnpairedScore[i] = std::max(bestUnpairedScore[i], alignmentPair[i].getScore());
        }
      }
      // const bool isPair   = pairMatch(insertSizeParameters, readPair, *seedChains[0], *seedChains[1]);
      const bool properPair = pairMatch(insertSizeParameters, alignmentPair[0], alignmentPair[1]);

      int       insert_len  = 0;
      int       insert_diff = 0;
      const int pair_pen    = pairBuilder.computePairPenalty(
          insertSizeParameters,
          readPair,
          &alignmentPair[0],
          &alignmentPair[1],
          properPair,
          insert_len,
          insert_diff);
      alignmentPair.setScore(alignmentPair[0].getScore() + alignmentPair[1].getScore() - pair_pen);
      alignmentPair.setPotentialScore(alignmentPair.getScore());
      alignmentPair.setProperPair(properPair);
      bestPairedScore = std::max(bestPairedScore, alignmentPair.getScore());
      //      DRAGEN_ALN_RESULT_LOG << "PAIR_RESULT:"
      //                            << " type=" << (alignmentPair.isSingleEnded() ? "se_pair" : "paired")
      //                            << ", paired=" << properPair << ", insert_len=" << insert_len
      //                            << ", insert_diff=" << insert_diff << ", pair_pen=" << pair_pen
      //                            << ", pair_score=" << alignmentPair.getScore() << std::endl;
    }
  }

  if (vectorizedSW_) {
    vectorSmithWaterman_.destroyReadContext(0);
    vectorSmithWaterman_.destroyReadContext(1);
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
    flags |= seedChain.isReverseComplement() ? Alignment::REVERSE_COMPLEMENT : Alignment::NONE;
    alignments.resize(alignments.size() + 1);
    auto& alignment = alignments.back();
    alignment.resetFlags(flags);
    // get the reference name from the ReferenceDir and the first diagonal
    const size_t referenceOffset      = seedChain.firstReferencePosition();
    const auto   referenceCoordinates = htConfig_.convertToReferenceCoordinates(referenceOffset);
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

}  // namespace align
}  // namespace dragenos
