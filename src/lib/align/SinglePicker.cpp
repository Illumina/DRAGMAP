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

#include "align/SinglePicker.hpp"

namespace dragenos {
namespace align {

/**
 * \brief find second best alignment
 */
Alignments::const_iterator findSecondBest(
    ScoreType snpCost, const Alignments& alignments, const Alignment* best, int& sub_count)
{
  Alignments::const_iterator ret = alignments.end();
  for (Alignments::const_iterator it = alignments.begin(); alignments.end() != it; ++it) {
    if (it->getIneligibilityStatus()) {
      // Single Picker is being used or paired data which may contain ineligibles
      continue;
    }
    if (!it->isUnmapped() && (alignments.end() == ret || ret->getScore() < it->getScore()) &&
        !best->isDuplicate(*it) && best->isOverlap(*it)) {
      ret = it;
    }
  }

  if (alignments.end() != ret) {
    //    std::cerr << "second best: " << *ret << std::endl;
    const ScoreType list_se_max = ret->getScore();
    const ScoreType list_se_min = list_se_max - snpCost;
    //  std::cerr << "list_se_max:" << list_se_max << std::endl;
    //  std::cerr << "list_se_min:" << list_se_min << std::endl;

    for (Alignments::const_iterator it = alignments.begin(); alignments.end() != it; ++it) {
      if (it->getIneligibilityStatus()) {
        // Single Picker is being used or paired data which may contain ineligibles
        continue;
      }
      if (!it->isUnmapped() && best != &*it &&
          (it->getScore() > list_se_min && it->getScore() <= list_se_max)) {
        //        std::cerr << "near sub:" << *it << std::endl;
        ++sub_count;
      }
    }
    assert(sub_count);
  }
  return ret;
}

ScoreType findSecondBestScore(
    ScoreType snpCost, const Alignments& alignments, const Alignment* best, int& sub_count)
{
  const Alignments::const_iterator secondBest = findSecondBest(snpCost, alignments, best, sub_count);

  return alignments.end() == secondBest ? INVALID_SCORE : secondBest->getScore();
}

Alignments::iterator SinglePicker::findSupplementary(
    const double readLength, Alignments& alignments, Alignment* primary) const
{
  Alignments::iterator ret = alignments.end();
  for (Alignments::iterator it = alignments.begin(); alignments.end() != it; ++it) {
    if (it->getIneligibilityStatus()) {
      // Single Picker is being used or paired data which may contain ineligibles
      continue;
    }
    if ((alignments.end() == ret || ret->getScore() < it->getScore()) && !primary->isDuplicate(*it) &&
        (!it->isUnmapped() && !primary->isOverlap(*it))) {
      ret = it;
    }
  }

  if (alignments.end() != ret &&
      suppMinScoreAdj_ <= (ret->getScore() - alnMinScore_))  // std::max(alnMinScore_, minScore)))
  {
    updateMapq(readLength, alignments, ret);
    return ret;
  } else {
    return alignments.end();
  }
}

/*
 * \return iterator of best scored alignment
 */
void SinglePicker::updateMapq(const int readLength, Alignments& alignments, Alignments::iterator best) const
{
  assert(!alignments.empty());

  int       sub_count       = 0;
  ScoreType secondBestScore = findSecondBestScore(similarity_.getSnpCost(), alignments, &*best, sub_count);
  //    std::cerr << "r0 sub_count:" << sub_count << std::endl;
  const MapqType a2m_mapq = computeMapq(
      similarity_.getSnpCost(),
      best->getScore(),
      std::max(alnMinScore_, secondBestScore),
      std::max(aln_cfg_mapq_min_len_, readLength));
  const int      sub_count_log2 = log2_simple(sub_count);
  const MapqType sub_mapq_pen_v = sub_count ? (3 * sub_count_log2) : 0;
#ifdef TRACE_SCORING
  std::cerr << "[SCORING]\t"
            << "r" << 0 << " a2m_mapq=" << a2m_mapq << " sub_count=" << sub_count
            << " sub_count_log2=" << sub_count_log2 << " sub_mapq_pen_v=" << sub_mapq_pen_v << std::endl;
#endif  // TRACE_SCORING
  const bool mapq0 = (aln_cfg_sample_mapq0_ >= 1 && best->hasOnlyRandomSamples()) ||
                     (aln_cfg_sample_mapq0_ >= 2 && best->isExtra());
  const MapqType mapq = mapq0 ? 0 : a2m_mapq - (sub_mapq_pen_v >> 7);

  best->setMapq(mapq);
  const ScoreType xs = secondBestScore >= alnMinScore_ ? secondBestScore : INVALID_SCORE;
  best->setXs(xs);
}

Alignments::iterator SinglePicker::pickBest(const int readLength, Alignments& alignments) const
{
  Alignments::iterator best = std::max_element(
      alignments.begin(), alignments.end(), [](const Alignment& left, const Alignment& right) {
        return left.getScore() < right.getScore();
      });

  // code below is more consistent with sim if the following is done to vhd:
  // changing ../../mapper/src/aln_result_proc.vhd around line 3938 to:
  //      inp_result.score.rand <= aln_res.ref_beg(score_rand_bits_c-1 downto 0);--lfsr(score_rand_bits_c-1 downto 0);

  //  assert(!alignments.empty());
  //
  //  Alignments::iterator best = alignments.begin();
  //  std::cerr << *best << " is best " << std::endl;
  //  for (Alignments::iterator it = best + 1; alignments.end() != it; ++it)
  //  {
  //    if (best->getScore() < it->getScore() || (best->getScore() == it->getScore() && (best->getPosition() &
  //    0xFFFF) < (it->getPosition()&0xFFFF)))
  //    {
  //      std::cerr << *it << " better than " << *best << std::endl;
  //      best = it;
  //    }
  //  }

  if (alignments.end() != best) {
    if (true == best->getIneligibilityStatus()) return alignments.end();

    updateMapq(readLength, alignments, best);
  }

  return best;
}

}  // namespace align
}  // namespace dragenos
