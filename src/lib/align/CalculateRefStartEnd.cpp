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
//#include <queue>

//#include <sys/mman.h>
//#include <fcntl.h>

//#include "common/Exceptions.hpp"
#include "align/CalculateRefStartEnd.hpp"
#include "align/SmithWaterman.hpp"
//#include "align/PairBuilder.hpp"
//#include "align/Mapq.hpp"
//#include "sequences/Seed.hpp"

namespace dragenos {
namespace align {

static const int            ALN_CFG_SW_EARLY_MAX = 256;
std::pair<int64_t, int64_t> calculateRefStartEnd(const sequences::Read& read, const map::SeedChain& chain)
{
  //    std::cerr << "calculateRefStartEnd:" << chain << std::endl;
  //  const auto beginPosition = chain.firstReferencePosition();// - 36;// - 12;
  //  const auto endPosition =  chain.lastReferencePosition() + 1;// + 36;// + 12;
  //  return std::make_pair(beginPosition, endPosition);

  const int64_t read_len = read.getLength();
  //  const int64_t beg_gap =
  //      chain.isPerfect()
  //          ? 0
  //          //: chain.front().getSeed().getReadPosition();  //chain.firstReadBase(); //chain.first_seed_pos
  //          : chain.firstReadBase(); //chain.first_seed_pos
  const int64_t beg_gap = chain.firstReadBase();  // chain.first_seed_pos
  //  const int64_t end_gap = chain.isPerfect() ? 0
  ////                                            : read_len - chain.back().getLastBaseReadPosition() -
  //                                            : read_len - chain.lastReadBase()/*last_seed_pos*/ - 1;

  const int64_t end_gap = read_len - chain.lastReadBase() /*last_seed_pos*/ - 1;

  const int64_t beg_gap_minus = beg_gap - ALN_CFG_SW_EARLY_MAX;
  const int64_t end_gap_minus = end_gap - ALN_CFG_SW_EARLY_MAX;

  const bool perfect_align = false;

  const bool beg_gap_small = perfect_align || beg_gap <= ALN_CFG_SW_EARLY_MAX;
  const bool end_gap_small = perfect_align || end_gap <= ALN_CFG_SW_EARLY_MAX;

  const int skip_gap =
      chain.isReverseComplement() ? (end_gap_small ? 0 : end_gap_minus) : (beg_gap_small ? 0 : beg_gap_minus);

  const int64_t head_gap = chain.isReverseComplement() ? end_gap : beg_gap;
  const int64_t tail_gap = chain.isReverseComplement() ? beg_gap : end_gap;

  //  std::cerr << "chain.lastReferencePosition()=" << chain.lastReferencePosition() <<
  //    " chain.firstReferencePosition()=" << chain.firstReferencePosition() << std::endl;

  //  -- Compute the chain's length in the reference, which corresponds to the chain's length in the read,
  //  -- but may be different due to indels.  Add to this the head and tail gaps between the chain and the
  //  -- beginning and end of the read, yielding the expected length of the reference segment corresponding
  //  -- to the full read.
  const int64_t tmp_ref_len =
      //      (chain.isPerfect() ? (chain.lastReferencePosition() - chain.firstReferencePosition()) :
      //      (chain.lastSeedReferencePosition() - chain.firstSeedReferencePosition())) +
      (chain.lastSeedReferencePosition() - chain.firstSeedReferencePosition()) + (head_gap + tail_gap + 1);

  static const int64_t sw_cells_x2p1_c = SmithWaterman::SW_CELLS * 2 + 1;
  //  -- Subtract from the 'full' reference segment the length we are skipping in the read, yielding
  //  -- a corresponding alignment length.  But then add 2*SW_CELLS for the reference window extension
  //  -- at each end, plus 1 dummy reference position to initialize S-W scoring.  Also add 1/8 of the
  //  -- tail gap, the potentially long region we plan to align beyond the chain, in case that region
  //  -- has additional deletions.
  const int64_t ref_length = tmp_ref_len - skip_gap + (tail_gap >> 3) + sw_cells_x2p1_c;

  static const int64_t sw_cells_p1_c = SmithWaterman::SW_CELLS + 1;
  //  -- We may be aligning extra bases before hitting the chain, namely head gap minus skipped bases;
  //  -- the reference alignent start needs to be similarly adjusted.  Also, we always include SW_CELLS
  //  -- extra bases for the reference, plus 1 more dummy reference position to initialize S-W scoring.
  const int64_t ref_start_adj = head_gap + sw_cells_p1_c;
  //  const int64_t ref_start     = /*ch_cmd.ref_reverse ?
  //        (chain.lastSeedReferencePosition() + ref_start_adj) :*/
  //      ((chain.isPerfect() ? chain.firstReferencePosition() : chain.firstSeedReferencePosition()) -
  //       ref_start_adj);
  //  //     +
  //  //    chain.isReverseComplement();

  const int64_t ref_start = /*ch_cmd.ref_reverse ?
    (chain.lastSeedReferencePosition() + ref_start_adj) :*/
      chain.firstSeedReferencePosition() - ref_start_adj;

  //  -- Compute the last reference base to fetch
  const int64_t ref_m1_v = ref_length - 1;

  const int64_t ref_end = /*ch_cmd.ref_reverse ? (ref_start - ref_m1_v) :*/ (ref_start + ref_m1_v);

  //    std::cerr << "beg_gap=" << beg_gap << ", end_gap=" << end_gap << ", skip_gap=" << skip_gap << ", tail_gap=" << tail_gap <<
  //      ", head_gap=" << head_gap << std::uppercase <<
  //      ", ref_start=0x" << std::hex << std::setw(10) << std::setfill('0') << ref_start <<
  //      ", ref_end=0x" << std::hex << std::setw(10) << std::setfill('0') << ref_end <<
  //      std::dec <<
  //      ", qry_length=" << read_len <<
  //      ", ref_length=" << ref_length << ", tmp_ref_len=" << tmp_ref_len <<
  //      ", ref_start_adj=" << ref_start_adj << chain << std::endl;

  return std::make_pair(ref_start, ref_end);
}

}  // namespace align
}  // namespace dragenos
