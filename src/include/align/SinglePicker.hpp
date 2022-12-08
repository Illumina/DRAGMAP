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

#ifndef ALIGN_SINGLE_PICKER_HPP
#define ALIGN_SINGLE_PICKER_HPP

#include <sys/mman.h>
#include <functional>
#include <vector>

#include <boost/filesystem.hpp>

#include "align/Alignments.hpp"
#include "align/InsertSizeParameters.hpp"
#include "map/Mapper.hpp"
#include "reference/Hashtable.hpp"
#include "reference/HashtableConfig.hpp"
#include "reference/ReferenceDir.hpp"
#include "sequences/Read.hpp"
#include "sequences/ReadPair.hpp"

namespace dragenos {
namespace align {

class SinglePicker {
  const SimilarityScores& similarity_;
  const int               alnMinScore_;
  const int               suppMinScoreAdj_;
  const int               aln_cfg_sec_aligns_;
  const int               aln_cfg_sec_score_delta_;
  const int               aln_cfg_sec_phred_delta_;
  const bool              aln_cfg_sec_aligns_hard_;
  const int               aln_cfg_mapq_min_len_;
  const int               aln_cfg_sample_mapq0_;

public:
  typedef sequences::Read   Read;
  typedef align::Alignment  Alignment;
  typedef align::Alignments Alignments;

  SinglePicker(
      const SimilarityScores& similarity,
      const int               alnMinScore,
      const int               suppMinScoreAdj,
      const int               aln_cfg_sec_aligns,
      const int               aln_cfg_sec_score_delta,
      const int               aln_cfg_sec_phred_delta,
      const bool              aln_cfg_sec_aligns_hard,
      const int               aln_cfg_mapq_min_len,
      const int               aln_cfg_sample_mapq0)
    : similarity_(similarity),
      alnMinScore_(alnMinScore),
      suppMinScoreAdj_(suppMinScoreAdj),
      aln_cfg_sec_aligns_(aln_cfg_sec_aligns),
      aln_cfg_sec_score_delta_(aln_cfg_sec_score_delta),
      aln_cfg_sec_phred_delta_(aln_cfg_sec_phred_delta),
      aln_cfg_sec_aligns_hard_(aln_cfg_sec_aligns_hard),
      aln_cfg_mapq_min_len_(aln_cfg_mapq_min_len),
      aln_cfg_sample_mapq0_(aln_cfg_sample_mapq0)
  {
  }

  Alignments::iterator pickBest(const int readLength, Alignments& alignments) const;
  Alignments::iterator findSupplementary(
      const double readLength, Alignments& alignments, Alignment* best) const;

  /**
   * \brief find secondary alignments
   *
   * \return false if aln_cfg_sec_aligns_hard_ is set and number of secondary alignments exceeds max allowed
   */
  template <typename StoreOp>
  bool findSecondary(const int readLength, Alignments& alignments, const Alignment* best, StoreOp store) const
  {
    const int m2a_scale = mapq2aln(similarity_.getSnpCost(), std::max(aln_cfg_mapq_min_len_, readLength));
    const ScoreType scaled_max_pen = (m2a_scale * aln_cfg_sec_phred_delta_) >> 10;  //27;
    const ScoreType sec_aln_delta  = std::max(scaled_max_pen, aln_cfg_sec_score_delta_);
    //    std::cerr << "scaled_max_pen:" << scaled_max_pen << std::endl;

    if (aln_cfg_sec_aligns_hard_) {
      int secAligns = 0;
      for (Alignments::iterator it = alignments.begin();
           secAligns < aln_cfg_sec_aligns_ + 1 && alignments.end() != it;
           ++it) {
        if (!it->isUnmapped() && alnMinScore_ <= it->getScore() && !best->isDuplicate(*it) &&
            best->isOverlap(*it) && sec_aln_delta >= (best->getScore() - it->getScore())) {
          ++secAligns;
        }
      }
      if (secAligns > aln_cfg_sec_aligns_) {
        return false;
      }
    }

    int secAligns = aln_cfg_sec_aligns_;
    for (Alignments::iterator it = alignments.begin(); secAligns && alignments.end() != it; ++it) {
      if (!it->isUnmapped() && alnMinScore_ <= it->getScore() && !best->isDuplicate(*it) &&
          best->isOverlap(*it) && sec_aln_delta >= (best->getScore() - it->getScore())) {
        store(*it);
        --secAligns;
      }
    }

    return true;
  }

private:
  void updateMapq(const int readLength, Alignments& alignments, Alignments::iterator best) const;
};

ScoreType findSecondBestScore(
    ScoreType snpCost, const Alignments& alignments, const Alignment* best, int& sub_count);

}  // namespace align
}  // namespace dragenos

#endif  // #ifndef ALIGN_SINGLE_PICKER_HPP
