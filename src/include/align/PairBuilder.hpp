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

#ifndef ALIGN_PAIR_BUILDER_HPP
#define ALIGN_PAIR_BUILDER_HPP

#include <sys/mman.h>
#include <functional>
#include <vector>

#include <boost/filesystem.hpp>

#include "align/AlignmentGenerator.hpp"
#include "align/Alignments.hpp"
#include "align/InsertSizeParameters.hpp"
#include "map/Mapper.hpp"
#include "reference/Hashtable.hpp"
#include "reference/HashtableConfig.hpp"
#include "reference/ReferenceDir.hpp"
#include "sequences/ReadPair.hpp"

namespace dragenos {
namespace align {

typedef std::array<Alignments, 2> UnpairedAlignments;

class PairBuilder {
  const SimilarityScores& similarity_;
  const int               alnMinScore_;
  const int               aln_cfg_unpaired_pen_;
  const int               aln_cfg_xs_pair_penalty_;
  const int               aln_cfg_sec_aligns_;
  const int               aln_cfg_sec_score_delta_;
  const int               aln_cfg_sec_phred_delta_;
  const bool              aln_cfg_sec_aligns_hard_;
  const int               aln_cfg_mapq_min_len_;
  const int               aln_cfg_sample_mapq0_;

  mutable std::vector<int> reported_;

public:
  typedef sequences::Read     Read;
  typedef sequences::ReadPair ReadPair;
  typedef align::Alignment    Alignment;
  typedef align::Alignments   Alignments;

  PairBuilder(
      const SimilarityScores& similarity,
      const int               alnMinScore,
      const uint32_t          aln_cfg_unpaired_pen,
      const int               aln_cfg_xs_pair_penalty,
      const int               aln_cfg_sec_aligns,
      const int               aln_cfg_sec_score_delta,
      const int               aln_cfg_sec_phred_delta,
      const bool              aln_cfg_sec_aligns_hard,
      const int               aln_cfg_mapq_min_len,
      const int               aln_cfg_sample_mapq0)
    : similarity_(similarity),
      alnMinScore_(alnMinScore),
      aln_cfg_unpaired_pen_(aln_cfg_unpaired_pen),
      aln_cfg_xs_pair_penalty_(aln_cfg_xs_pair_penalty),
      aln_cfg_sec_aligns_(aln_cfg_sec_aligns),
      aln_cfg_sec_score_delta_(aln_cfg_sec_score_delta),
      aln_cfg_sec_phred_delta_(aln_cfg_sec_phred_delta),
      aln_cfg_sec_aligns_hard_(aln_cfg_sec_aligns_hard),
      aln_cfg_mapq_min_len_(aln_cfg_mapq_min_len),
      aln_cfg_sample_mapq0_(aln_cfg_sample_mapq0)
  {
  }

  AlignmentPairs::iterator pickBest(
      const ReadPair&           readPair,
      AlignmentPairs&           alignmentPairs,
      const UnpairedAlignments& unpairedAlignments) const;

  void updateMapq(
      const int                 readLength,
      AlignmentPairs&           alignments,
      const UnpairedAlignments& unpairedAlignments,
      AlignmentPairs::iterator  best) const;

  int computePairPenalty(
      const InsertSizeParameters& insertSizeParameters,
      const ReadPair&             readPair,
      const Alignment*            a1,
      const Alignment*            a2,
      const bool                  properPair,
      int&                        insert_len,
      int&                        insert_diff) const;

  /**
   * \brief     Find all secondary alignments within config parameters. Return false iff sec-aligns is
   * exceeded and sec-aligns-hard is set
   */
  template <typename StoreOp>
  bool findSecondary(
      const int                            averageReadLength,
      AlignmentPairs&                      pairs,
      const align::Alignments&             unpaired,
      const AlignmentPairs::const_iterator best,
      int                                  readIdx,
      StoreOp                              store) const
  {
    const int m2a_scale =
        mapq2aln(similarity_.getSnpCost(), std::max(aln_cfg_mapq_min_len_, averageReadLength));
    const ScoreType scaled_max_pen = (m2a_scale * aln_cfg_sec_phred_delta_) >> 10;  //27;
    const ScoreType sec_aln_delta  = std::max(scaled_max_pen, aln_cfg_sec_score_delta_);

    reported_.clear();
    reported_.resize(unpaired.size(), false);
    int secAligns = aln_cfg_sec_aligns_;

    for (auto& p : pairs) {
      if (/*(!best->isProperPair() || p.isProperPair()) && */ !p.cat(readIdx).isUnmapped() &&
          alnMinScore_ <= p.at(readIdx).getScore() && sec_aln_delta >= (best->getScore() - p.getScore()) &&
          !best->at(readIdx).isDuplicate(p.at(readIdx)) &&
          // D0004:230:H08B1ADXX:1:1105:18508:35343          (p.at(readIdx).isUnmapped() ||
          // best->at(readIdx).isOverlap(p.at(readIdx))) &&
          !reported_[std::distance(&unpaired.front(), &p.cat(readIdx))]) {
        if (!secAligns) {
          return !aln_cfg_sec_aligns_hard_;
        }
        store(p);
        --secAligns;
        reported_[std::distance(&unpaired.front(), &p.cat(readIdx))] = true;
      }
    }

    return true;
  }

private:
  AlignmentPairs::const_iterator findSecondBest(
      const int                            averageReadLength,
      const AlignmentPairs&                alignments,
      const UnpairedAlignments&            unpairedAlignments,
      const AlignmentPairs::const_iterator best,
      int                                  readIdx,
      int&                                 sub_count) const;

  AlignmentPairs::const_iterator findSecondBestScore(
      const int                            averageReadLength,
      const AlignmentPairs&                pairs,
      const UnpairedAlignments&            unpairedAlignments,
      const AlignmentPairs::const_iterator best,
      int                                  readIdx,
      int&                                 sub_count,
      ScoreType&                           secondBestPeScore,
      ScoreType&                           secondBestSeScore) const;

  void updateEndMapq(
      const int                averageReadLength,
      AlignmentPairs::iterator best,
      const int                readIdx,
      int                      sub_count,
      ScoreType                sub_pair_score_v,
      const ScoreType          secondBestSeScore[],
      const ScoreType          xs[]) const;
};

}  // namespace align
}  // namespace dragenos

#endif  // #ifndef ALIGN_PAIR_BUILDER_HPP
