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

#ifndef ALIGN_ALIGNER_HPP
#define ALIGN_ALIGNER_HPP

#include "align/AlignmentGenerator.hpp"
#include "align/AlignmentRescue.hpp"
#include "align/PairBuilder.hpp"
#include "reference/Hashtable.hpp"
#include "reference/ReferenceDir.hpp"
#include "sequences/Read.hpp"
#include "sequences/ReadPair.hpp"

namespace dragenos {
namespace align {

class Aligner {
public:
  Aligner()               = delete;
  Aligner(const Aligner&) = delete;
  Aligner& operator=(const Aligner&) = delete;
  // Aligner(Aligner&&) = delete;
  //Aligner &operator=(Aligner&&) = delete;
  Aligner(
      const reference::ReferenceSequence& refSeq,
      const reference::HashtableConfig&   htConfig,
      const reference::Hashtable&         hashtable,
      bool                                mapOnly,
      const int                           swAll,
      const SimilarityScores&             similarity,
      const int                           gapInit,
      const int                           gapExtend,
      const int                           unclipScore,
      const int                           alnMinScore,
      const int                           aln_cfg_mapq_min_len,
      const uint32_t                      alignerUnpairedPen,
      const double                        aln_cfg_filter_len_ratio,
      const bool                          vectorizedSW);
  typedef sequences::Read     Read;
  typedef sequences::ReadPair ReadPair;
  typedef align::Alignment    Alignment;
  typedef align::Alignments   Alignments;
  //typedef std::vector<Alignment> Alignments;
  //typedef align::AlignmentPair AlignmentPair;
  //typedef std::vector<AlignmentPair> AlignmentPairs;
  void                     getAlignments(const Read& read, Alignments& alignments);
  AlignmentPairs::iterator getAlignments(
      const ReadPair&             readPair,
      AlignmentPairs&             alignmentPairs,
      const InsertSizeParameters& insertSizeParameters,
      const PairBuilder&          pairBuilder);
  Alignments& unpaired(std::size_t readPosition) { return unpairedAlignments_.at(readPosition); }
  /// generate ungapped alignments from the seed chains
  void generateUngappedAlignments(const Read& read, map::ChainBuilder& chainBuilder, Alignments& alignments);
  void runSmithWatermanAll(
      const Read& read, const map::ChainBuilder& chainBuilder, Alignments& alignments, const int readIdx);
  /// generate the ungapped alignment for the given seed chain
  void generateUngappedAlignment(const Read& read, map::SeedChain& seedChain, Alignment& alignment);
  /// generate dummy alignments consistent with the seed chains
  void generateDummyAlignments(
      const Read& read, const map::ChainBuilder& chainBuilder, Alignments& alignments) const;
  void deFilterChain(const Read& read, map::SeedChain& seedChain, Alignment& alignment, const int readIdx);
  /**
   ** \brief Accepts query and reference sequences of identical length and compares
   ** corresponding bases one by one.
   **
   ** Scores matches and mismatches and counts mismatches as usual for regular bases
   ** using aln_cfg.mismatch_score and aln_cfg.match_score.
   ** Uses a special score (aln_cfg.match_n_score) for Ns and counts Ns in either
   ** reference or query as mismatches.
   **
   ** If there is a mismatch at the beginning or at the end, one base may be clipped
   ** if the mismatch penalty minus the match score exceeds the unclipped bonus
   ** (aln_cfg.unclip_score). Only a single base can be clipped from either end.
   **
   ** Bursts of SNPs are detected when there are more than a given number of SNPs
   ** (burst_minimum_c) within a window of fixed size (burst_window_c).
   **
   ** To be a perfect alignment, the query and database must be the same length and there
   ** shouldn't be any burts of SNPs.
   **
   ** \param queryBegin
   ** \param queryEnd
   ** \param databaseBegin
   ** \param databaseEnd
   ** \param alignment output parameter where the appropriate alignment information is stored
   ** \return false if Smith-Waterman is required, true otherwise.
   **/
  static bool isPerfectAlignment(
      const char* queryBegin,
      const char* queryEnd,
      const char* databaseBegin,
      const char* databaseEnd,
      Alignment&  alignment);
  /// calculate the ungapped alignment score and potential score for the given read at the specified
  /// orientation and position
  int initializeUngappedAlignmentScores(
      const Read& read, const bool rcFlag, const size_t referenceOffset, Alignment& alignment);
  /// calculate potential score, assuming that the soft clip regions are improved with a gap
  unsigned calculatePotentialScore(
      const Read& read, const unsigned score, const unsigned softClipTo, const unsigned softClipFrom) const;
  /// filter overlapping alignments sharing a start or end position
  void filter(Alignments& alignments);
  /**
   ** filter overlapping pairs
   **
   ** first identified for each read if there are overlapping alignments ansd filters them out. Then
   *identifies
   ** the pairs that are built with single alignments that hav been filtered out.
   ** Single alignments are overlapping if either the beginning or the end of the reads map at the same
   *position.
   **/
  void filter(AlignmentPairs& alignmentPairs, std::array<Alignments, 2>& unpairedAlignments);

private:
  const reference::ReferenceSequence& refSeq_;
  const reference::HashtableConfig&   htConfig_;
  const bool                          mapOnly_;
  const int                           swAll_;
  const bool                          vectorizedSW_;
  /// read the hashtable config data and throw on error
  //std::vector<char> getHashtableConfigData(const boost::filesystem::path referenceDir) const;
  /// maps hashtable data and throw on error
  //uint64_t *mmapHashtableData(const boost::filesystem::path referenceDir) const;
  // raw binary content of the hashtable config file
  //const std::vector<char> hashtableConfigData_;
  // TODO: replace with a placement new
  //const reference::HashtableConfig hashtableConfig_;
  const map::Mapper         mapper_;
  const SimilarityScores    similarity_;
  const ScoreType           gapInit_;
  // const ScoreType           gapExtend_;
  // const ScoreType           unclipScore_;
  const ScoreType           alnMinScore_;
  const int                 aln_cfg_mapq_min_len_;
  const int                 aln_cfg_unpaired_pen_;
  SmithWaterman             smithWaterman_;
  VectorSmithWaterman       vectorSmithWaterman_;
  std::array<Alignments, 2> unpairedAlignments_;
  AlignmentGenerator        alignmentGenerator_;

  std::array<map::ChainBuilder, 2> chainBuilders_;

  /// generate all the ungapped allignments for the seed chains
  void buildUngappedAlignments(map::ChainBuilder& chainBuilder, const Read& read, Alignments& alignments);

  bool rescueMate(
      const InsertSizeParameters& insertSizeParameters,
      const Read&                 anchoredRead,
      const Read&                 rescuedRead,
      const map::SeedChain&       anchoredSeedChain,
      Alignment&                  anchoredAlignment,
      const AlignmentRescue&      alignmentRescue,
      map::SeedChain&             rescuedSeedChain,
      Alignment&                  rescuedAlignment,
      const int                   anchoredIdx);

  void makePair(
      const InsertSizeParameters& insertSizeParameters,
      const ReadPair&             readPair,
      const PairBuilder&          pairBuilder,
      const map::SeedChain*       s0,
      const map::SeedChain*       s1,
      Alignment&                  a0,
      Alignment&                  a1,
      AlignmentPairs&             alignmentPairs);

  bool rescuePair(
      const InsertSizeParameters& insertSizeParameters,
      const PairBuilder&          pairBuilder,
      const ReadPair&             readPair,
      const int                   anchoredIdx,
      const bool                  any_pair_match,
      const map::SeedChain&       anchoredSeedChain,
      Alignment&                  anchoredAlignment,
      const AlignmentRescue       alignmentRescue,
      AlignmentPairs&             alignmentPairs);
  void runSmithWatermanWorthy(const Read& read, map::ChainBuilder& chainBuilder, Alignments& alignments);
  void updateIneligibility(
      const reference::HashtableConfig& hashtableConfig,
      const size_t                      referenceOffset,
      const Read&                       read,
      Alignment&                        alignment);
};

}  // namespace align
}  // namespace dragenos

#endif  // #ifndef ALIGN_ALIGNER_HPP
