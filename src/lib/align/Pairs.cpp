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

#include "align/AlignmentRescue.hpp"
#include "align/PairBuilder.hpp"
#include "align/Tlen.hpp"

namespace dragenos {
namespace align {

std::pair<int, int> calculateEffBegEnd(const sequences::Read& read, const map::SeedChain& chain)
{
  //  std::cerr << chain << std::endl;
  const int beg_gap_v = chain.firstReadBase();
  const int end_gap_v = read.getLength() - chain.lastReadBase() - 1;

  if (chain.isReverseComplement()) {
    const int pipe3_eff_beg = chain.lastSeedReferencePosition() + beg_gap_v;
    const int pipe3_eff_end = chain.firstSeedReferencePosition() - end_gap_v;
    return std::make_pair(pipe3_eff_beg, pipe3_eff_end);
  } else {
    const int pipe3_eff_beg = chain.firstSeedReferencePosition() - beg_gap_v;
    const int pipe3_eff_end = chain.lastSeedReferencePosition() + end_gap_v;
    return std::make_pair(pipe3_eff_beg, pipe3_eff_end);
  }
}

bool pairMatch(
    const InsertSizeParameters& insertSizeParameters,
    const sequences::ReadPair&  readPair,
    const map::SeedChain&       first,
    const map::SeedChain&       second)
{
  using Orientation = InsertSizeParameters::Orientation;
  bool end1_pair_rc = false;
  if (insertSizeParameters.orientation_ == Orientation::pe_orient_fr_c ||
      insertSizeParameters.orientation_ == Orientation::pe_orient_rf_c) {
    end1_pair_rc = !first.isReverseComplement();
  } else {
    end1_pair_rc = first.isReverseComplement();
  }

  const std::pair<int, int> r1_eff_beg_end = calculateEffBegEnd(readPair[0], first);

  int beg_pm_pair_min_v = 0;
  int beg_pm_pair_max_v = 0;
  int end_pm_pair_min_v = 0;
  int end_pm_pair_max_v = 0;

  if (first.isReverseComplement()) {
    beg_pm_pair_min_v = r1_eff_beg_end.first - insertSizeParameters.max_;
    beg_pm_pair_max_v = r1_eff_beg_end.first - insertSizeParameters.min_;
    end_pm_pair_min_v = r1_eff_beg_end.second + insertSizeParameters.min_;
    end_pm_pair_max_v = r1_eff_beg_end.second + insertSizeParameters.max_;
  } else {
    beg_pm_pair_min_v = r1_eff_beg_end.first + insertSizeParameters.min_;
    beg_pm_pair_max_v = r1_eff_beg_end.first + insertSizeParameters.max_;
    end_pm_pair_min_v = r1_eff_beg_end.second - insertSizeParameters.max_;
    end_pm_pair_max_v = r1_eff_beg_end.second - insertSizeParameters.min_;
  }

  int end1_pair_min_beg = 0;
  int end1_pair_max_beg = 0;
  int end1_pair_min_end = 0;
  int end1_pair_max_end = 0;

  switch (insertSizeParameters.orientation_) {
    //  -- Forward-reverse: measure from beginning to beginning only
  case Orientation::pe_orient_fr_c:
    end1_pair_min_beg = beg_pm_pair_min_v;
    end1_pair_max_beg = beg_pm_pair_max_v;
    end1_pair_min_end = -1;  //unsigned(-1) >> 1;//(ref_pos_bits_c => '0', others => '1');
    end1_pair_max_end =
        0;  //unsigned(-1) << (sizeof(end1_pair_max_end) * 8 - 1);//(ref_pos_bits_c => '1', others => '0');
    break;
    //  -- Reverse-forward: measure from end to end only
  case Orientation::pe_orient_rf_c:
    end1_pair_min_beg = -1;  //-1 >> 1;//(ref_pos_bits_c => '0', others => '1');
    end1_pair_max_beg =
        0;  //-1 << (sizeof(end1_pair_max_end) * 8 - 1);//(ref_pos_bits_c => '1', others => '0');
    end1_pair_min_end = end_pm_pair_min_v;
    end1_pair_max_end = end_pm_pair_max_v;
    break;
    //  -- Forward-forward (or equivalently, reverse-reverse): measure beginning to end OR end to beginning
  default:
    end1_pair_min_beg = end_pm_pair_min_v;
    end1_pair_max_beg = end_pm_pair_max_v;
    end1_pair_min_end = beg_pm_pair_min_v;
    end1_pair_max_end = beg_pm_pair_max_v;
    break;
  }

  const std::pair<int, int> r2_eff_beg_end = calculateEffBegEnd(readPair[0], second);

  //  -- Pair match if orientation as expected, and either the effective beginning or end is in the expected
  //  inteval
  const bool pair_match =
      !(end1_pair_rc ^ second.isReverseComplement()) &&
      (((r2_eff_beg_end.first >= end1_pair_min_beg) && (r2_eff_beg_end.first <= end1_pair_max_beg)) ||
       ((r2_eff_beg_end.second >= end1_pair_min_end) && (r2_eff_beg_end.first <= end1_pair_max_end)));

  //    std::cerr << "end1_pair_rc:" << end1_pair_rc << " second.isReverseComplement():" << second.isReverseComplement() <<
  //      " pipe3_eff_beg:" << r2_eff_beg_end.first << " pipe3_eff_end:" << r2_eff_beg_end.second <<
  //      " end1_pair_min_beg:" << end1_pair_min_beg << " end1_pair_max_beg:" << end1_pair_max_beg <<
  //      " end1_pair_min_end:" << end1_pair_min_end << " end1_pair_max_end:" << end1_pair_max_end <<
  //      " pair_match:" << pair_match <<
  //      (((r2_eff_beg_end.first >= end1_pair_min_beg) && (r2_eff_beg_end.first <= end1_pair_max_beg)) ? " yes" : " no") <<
  //      (((r2_eff_beg_end.second >= end1_pair_min_end) && (r2_eff_beg_end.first <= end1_pair_max_end)) ? " yes" : " no") << std::endl;

  //  if (pair_match)
  //  {
  //    std::cerr << "pairMatch\n" << first << "\n" << second << std::endl;
  //  }

  return pair_match;
}

/**
 * \brief State machine states: S_EFF_EDGES S_PICK_EDGES S_INSERT_LEN S_INSERT_DIFF
 */
bool pairMatch(const InsertSizeParameters& insert_stats, const Alignment& a1, const Alignment& a2)
{
  if (a1.isUnmapped() || a2.isUnmapped() || a1.getReference() != a2.getReference()) {
    return false;
  }

  const Alignment& result_rrec = a1.isReverseComplement() ? a2 : a1;
  const Alignment& inp_result  = a1.isReverseComplement() ? a1 : a2;

  //  -- Compute effective begin and end reference positions for the input and recent-other-end alignments,
  //  -- extrapolating the alignment begin/end positions by the read begin or end gaps, per the orientation.
  const int inp_eff_beg =
      (!inp_result.isReverseComplement() ? inp_result.getUnclippedEndPosition()
                                         : inp_result.getUnclippedStartPosition());
  const int inp_eff_end =
      (inp_result.isReverseComplement() ? inp_result.getUnclippedEndPosition()
                                        : inp_result.getUnclippedStartPosition());
  const int rct_eff_beg =
      (result_rrec.isReverseComplement() ? result_rrec.getUnclippedEndPosition()
                                         : result_rrec.getUnclippedStartPosition());
  const int rct_eff_end =
      (!result_rrec.isReverseComplement() ? result_rrec.getUnclippedEndPosition()
                                          : result_rrec.getUnclippedStartPosition());
  //
  //  -- Pick the appropriate alignment edges to represent the insert begin and end,
  //  -- per the expected paired end orientation (RF, FR, FF)
  //    -- Forward-forward orientation (same as reverse-reverse): pick the outermost begin and end
  //    if insert_stats.pe_orientation = pe_orient_ff_c or insert_stats.pe_orientation = pe_orient_rr_c then
  int insert_beg  = 0;
  int insert_end  = 0;
  int wrong_beg_v = 0;
  int wrong_end_v = 0;
  if (InsertSizeParameters::Orientation::pe_orient_ff_c == insert_stats.orientation_ ||
      InsertSizeParameters::Orientation::pe_orient_rr_c == insert_stats.orientation_) {
    insert_beg  = std::min(inp_eff_beg, rct_eff_beg);
    insert_end  = std::max(inp_eff_end, rct_eff_end);
    wrong_beg_v = inp_result.getPosition();
    wrong_end_v = inp_result.getEndPosition();
  }
  //    -- Forward-reverse: beginning from the forward one, end from the reverse one
  //    -- Reverse-forward: beginning from the reverse one, end from the forward one
  else if (
      (InsertSizeParameters::Orientation::pe_orient_fr_c == insert_stats.orientation_) ^
      inp_result.isReverseComplement()) {
    insert_beg  = inp_eff_beg;
    insert_end  = rct_eff_end;
    wrong_beg_v = result_rrec.getPosition();
    wrong_end_v = inp_result.getEndPosition();
  } else {
    insert_beg  = rct_eff_beg;
    insert_end  = inp_eff_end;
    wrong_beg_v = inp_result.getPosition();
    wrong_end_v = result_rrec.getEndPosition();
  }
  //    -- In RNA mode, it's possible for an alignment to nest inside an intron of the mate alignment.
  //    -- This shouldn't be considered a proper pair, but we won't notice by considering only
  //    -- the chosen insert_beg and insert_end; we need to detect if the other begin or end
  //    -- extends significantly further outward.
  static const int template_overshoot_lim_c = 6;
  const int        tem_max_beg              = wrong_beg_v + template_overshoot_lim_c;
  const int        tem_min_end              = wrong_end_v - template_overshoot_lim_c;
  const int        insert_len               = insert_end - insert_beg + 1;

  //    -- Detect suspicious 3' overshoot of mate's 5' end
  //    template_bad  <= (test_gt_signed(insert_beg,tem_max_beg) or test_lt_signed(insert_end,tem_min_end))
  //                     and test_ne(insert_stats.pe_orientation,pe_orient_ff_c)
  //                     and test_ne(insert_stats.pe_orientation,pe_orient_rr_c)
  //                     and not aln_cfg.global;
  const bool template_bad = (insert_beg > tem_max_beg || insert_end < tem_min_end) &&
                            InsertSizeParameters::Orientation::pe_orient_ff_c != insert_stats.orientation_ &&
                            InsertSizeParameters::Orientation::pe_orient_rr_c != insert_stats.orientation_;

  //
  //    -- For RF or FR the two orientations must be opposite; for FF they must be the same.
  //    -- Furthermore the calculated insert length (which may be negative) must be within the configured
  //    limits,
  //    -- and both ends must be in the same reference sequence.
  //    if ((inp_result.rev_comp xor result_rrec.rev_comp xor
  //    test_eq(insert_stats.pe_orientation,pe_orient_ff_c))
  //            and test_gte_signed(insert_len,EXT(insert_stats.pe_min_insert,seq_pos_bits_c+1))
  //            and test_lte_signed(insert_len,EXT(insert_stats.pe_max_insert,seq_pos_bits_c+1))
  //            and test_eq(inp_result.ref_id,result_rrec.ref_id)
  //            and not template_bad) = '0' then
  //      -- No pair score valid
  //      state <= S_SE_COMP;
  //    end if;

  //  std::cerr << "S_INSERT_DIFF:, insert_len=" << insert_len << ", template_bad=" << template_bad <<
  //    ", tem_max_beg=" << tem_max_beg << ", tem_min_end=" << tem_min_end <<
  //    ", insert_beg=" << insert_beg << ", insert_end=" << insert_end <<
  //    ", inp_eff_beg=" << inp_eff_beg << ", inpt_eff_end=" << inp_eff_end <<
  //    ", rct_eff_beg=" << rct_eff_beg << ", rct_eff_end=" << rct_eff_end << std::endl;

  return (inp_result.isReverseComplement() xor result_rrec.isReverseComplement() ^
          (InsertSizeParameters::Orientation::pe_orient_ff_c == insert_stats.orientation_)) &&
         insert_len >= insert_stats.min_ && insert_len <= insert_stats.max_ &&
         inp_result.reference_ == result_rrec.reference_ && !template_bad;
}

}  // namespace align
}  // namespace dragenos
