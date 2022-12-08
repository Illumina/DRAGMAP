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

#include <emmintrin.h>
#include <boost/assert.hpp>
#include <iomanip>

#include "align/AlignmentRescue.hpp"
#include "common/DragenLogger.hpp"
#include "map/SeedPosition.hpp"

// Created and modified based on adam's legacy code in align/remove/AlignmentPlan
namespace dragenos {
namespace align {

bool AlignmentRescue::triggeredBy(const SeedChain& anchoredChain, const bool any_pair_match) const
{
  const int chain_len = anchoredChain.getReadSpanLength();
  return anchoredChain.isExtra() ||
         ((any_pair_match) ? (chain_len >= resc_ifpair_len_) : (chain_len >= resc_nopair_len_));

  /*
  if rescue_en and (
    rescuable(bit_to_int(not second))
      and (
        (
          (chn_rrec.chain.extra or (chn_rrec.chain.sample and aln_cfg.rescue_hifreq)) and not chn_rrec.paired)
  or (enable_rescue and (not chn_rrec.chain.filtered_out) and ((chn_rrec.rescue(0) and not any_pair_match) or
  (chn_rrec.rescue(1) and any_pair_match)))) ) = '1' then
 */
}

// clears referenceBases before data is added
uint32_t AlignmentRescue::getReferenceInterval(
    const SeedChain&                    anchoredChain,
    const int                           rescuedReadLength,
    const reference::ReferenceSequence& reference,
    std::vector<unsigned char>&         referenceBases,
    const bool /*second*/) const
{
  typedef map::SeedPosition::ReferencePosition ReferencePosition;
  const bool                                   rescue_end = anchoredChain.isReverseComplement();
  const ReferencePosition                      pos_v =
      rescue_end ? anchoredChain.lastReferencePosition() : anchoredChain.firstReferencePosition();

  ReferencePosition rescue_min = 0;
  ReferencePosition rescue_max = 0;
  if (anchoredChain.isReverseComplement()) {
    rescue_min = pos_v - pe_min_insert_;
    rescue_max = pos_v - pe_max_insert_;
  } else {
    rescue_min = pos_v + pe_min_insert_;
    rescue_max = pos_v + pe_max_insert_;
  }

  // dragen actually does -1, but then traced SCAN_CMD mismatches
  const int other_len_m1 = rescuedReadLength + 1;
  // NOTE: only Illumina pairs are supported: they look like so: >>>>>>>>>>>>>    <<<<<<<<<<<<<<<<
  assert(Orientation::pe_orient_fr_c == pe_orientation_);
  ReferencePosition first_ref_pos = 0;
  ReferencePosition last_ref_pos  = 0;
  if (!rescue_end) {
    last_ref_pos  = rescue_max;
    first_ref_pos = rescue_min - other_len_m1;
  } else {
    first_ref_pos = rescue_max;
    last_ref_pos  = rescue_min + other_len_m1;
  }

  const int ref_length = ((((last_ref_pos + 1 - first_ref_pos) + 1) >> 2) << 2);
  //  const int ref_length = ((((last_ref_pos + 1 - first_ref_pos) + 3) >> 2) << 2);
  //  const int ref_length = (last_ref_pos + 1 - first_ref_pos);
  assert(first_ref_pos <= last_ref_pos);
  referenceBases.clear();
  referenceBases.reserve(ref_length + 1);
  if (isReversedRescue(anchoredChain)) {
    reference.getRcBases(last_ref_pos - ref_length + 1, last_ref_pos + 1, referenceBases);
    std::reverse(referenceBases.begin(), referenceBases.end());
  } else {
    reference.getBases(first_ref_pos, first_ref_pos + ref_length, referenceBases);
  }

  // std::cerr << "after=" << after << ", anchorPosition=" << anchorPosition<< ", length=" << length << ", pe_min_insert_=" << pe_min_insert_ << ", rescuedReadLength=" << rescuedReadLength << ", minPosition= " << minPosition << ", referenceBases.size()=" << referenceBases.size() << ", isReversedRescue(anchoredChain): " << isReversedRescue(anchoredChain) << std::endl;
  {
    boost::io::ios_flags_saver ifs(std::cerr);

    DRAGEN_RESCUE_LOG << "SCAN_CMD:  qry_start=0, ref_start="
                      << "0x" << std::uppercase << std::hex << std::setw(10) << std::setfill('0')
                      << (!isReversedRescue(anchoredChain) ? first_ref_pos : (last_ref_pos))
                      << ", qry_end=" << std::dec << (((rescuedReadLength >> 2) << 2) - 1) << ", ref_end="
                      << "0x" << std::uppercase << std::hex << std::setw(10) << std::setfill('0')
                      << (isReversedRescue(anchoredChain) ? (last_ref_pos - ref_length + 1)
                                                          : (first_ref_pos + ref_length - 1))
                      << ", qry_length=" << std::dec << rescuedReadLength
                      << ", ref_length=" << ref_length  // referenceBases.size()
                      << ", qry_comp=" << isReversedRescue(anchoredChain)
                      << ", ref_reverse=" << isReversedRescue(anchoredChain) << std::endl;
    //    std::cerr << anchoredChain << " pos_v=" << std::hex << pos_v <<" rescue_min=" << std::hex << rescue_min << " rescue_max=" << std::hex << rescue_max << std::endl;
  }

  return first_ref_pos;
}

__m128i mm_bitshift_left4(__m128i x)
{
  __m128i carry =
      _mm_bslli_si128(x, 8);          // old compilers only have the confusingly named _mm_slli_si128 synonym
  carry = _mm_srli_epi64(carry, 60);  // After bslli shifted left by 64b

  x = _mm_slli_epi64(x, 4);
  return _mm_or_si128(x, carry);
}

static inline int popcnt128(__m128i n)
{
  const __m128i n_hi = _mm_unpackhi_epi64(n, n);
  return __builtin_popcountll(_mm_cvtsi128_si64(n)) + __builtin_popcountll(_mm_cvtsi128_si64(n_hi));
}

bool AlignmentRescue::scan(
    const Read&                         rescuedRead,
    const SeedChain&                    anchoredChain,
    const reference::ReferenceSequence& reference,
    SeedChain&                          rescuedChain) const
{
  rescuedChain.clear();

  std::vector<unsigned char> referenceBases;
  if (((pe_orientation_ == Orientation::pe_orient_fr_c) ||
       (pe_orientation_ == Orientation::pe_orient_rf_c))) {
    // get the rescue reference interval
    const auto rescueIntervalStart =
        getReferenceInterval(anchoredChain, rescuedRead.getLength(), reference, referenceBases);

    // build the two rescue 32-mers
    unsigned   modOffset   = rescuedRead.getLength() % BASES_PER_CYCLE;
    const auto rescueKmers = getRescueKmers(rescuedRead, modOffset);

    // if (log)
    //{
    // for (int i= 0; i < 2; ++i)
    //{
    //std::cerr << "kmer " << i << ":";
    //for(auto b: rescueKmers[i]) std::cerr << " " << (int)b;
    //std::cerr << std::endl;
    //}
    //}

    // scan the rescue reference interval for each of the rescue 32-mers
    assert(pe_max_insert_ >= pe_min_insert_);
    const int scanLength = std::min(
        referenceBases.size() - rescuedRead.getLength(),
        referenceBases.size() - rescueKmers[1].size() - modOffset);
    // the offsets that give the least number of mismatches for each kmer
    int bestOffsets[] = {scanLength, scanLength};
    int bestCounts[]  = {rescueKmers[0].size(), rescueKmers[1].size()};
    const std::vector<unsigned char>::const_iterator startIterators[] = {
        referenceBases.begin(), referenceBases.end() - scanLength - rescueKmers[1].size() - modOffset};
    bool conflict = false;

    // if (log)
    //   std::cerr << "referenceBases.size(): " << referenceBases.size() << ", scanLength=" << scanLength
    //           << ", rescueKmers[1].size()=" << rescueKmers[1].size() << std::endl;

    __m128i kmer1_m = loadRescueKmer(rescueKmers[0]);
    __m128i kmer2_m = loadRescueKmer(rescueKmers[1]);

    __m128i refkmer1_m = loadRefKmer(startIterators[0], rescueKmers[0].size());
    __m128i refkmer2_m = loadRefKmer(startIterators[1], rescueKmers[0].size());

    for (int i = 0; scanLength > i; ++i) {
      // get the number of mismatches for both kmer
      //    const int mismatchCounts[] = {
      //        countMismatches(rescueKmers[0], startIterators[0] + i),
      //        countMismatches(rescueKmers[1], startIterators[1] + i)};

      // count matches = compare with &, then popcount
      __m128i comparek1 = _mm_and_si128(refkmer1_m, kmer1_m);
      __m128i comparek2 = _mm_and_si128(refkmer2_m, kmer2_m);
      int     nbmis0    = rescueKmers[0].size() - popcnt128(comparek1);
      int     nbmis1    = rescueKmers[1].size() - popcnt128(comparek2);

      const int mismatchCounts[] = {nbmis0, nbmis1};

      // update best counts and best offsets
      for (int j = 0; j != 2; ++j) {
        if (bestCounts[j] > mismatchCounts[j] or
            (bestCounts[j] == mismatchCounts[j] and i <= scanLength / 2)) {
          // flag a conflic if the kmer maps at multiple locations
          conflict |= (bestCounts[j] <= RESCUE_MAX_SNPS);
          bestCounts[j]  = mismatchCounts[j];
          bestOffsets[j] = i;
        }
      }
      // if (log) std::cerr << i << '\t' << mismatchCounts[0] << ":" << mismatchCounts[1] << " - " << bestCounts[0] << ":" << bestCounts[1] << " - " << bestOffsets[0] << ":" << bestOffsets[1] << std::endl;

      unsigned char newBase;
      // shift kmer window on reference by one base : shift left vector 4 bits and insert new base
      refkmer1_m = mm_bitshift_left4(refkmer1_m);
      newBase    = *(startIterators[0] + rescueKmers[0].size() - 1 + i + 1);
      if (newBase == 0xF) newBase = 0;
      __m128i newBase_m = _mm_setr_epi8(newBase, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
      refkmer1_m        = _mm_or_si128(refkmer1_m, newBase_m);

      // idem for second kmer
      refkmer2_m = mm_bitshift_left4(refkmer2_m);
      newBase    = *(startIterators[1] + rescueKmers[0].size() - 1 + i + 1);
      if (newBase == 0xF) newBase = 0;
      newBase_m  = _mm_setr_epi8(newBase, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
      refkmer2_m = _mm_or_si128(refkmer2_m, newBase_m);
    }

    // flag a conflict and adjust the best offset if they are on different diagonals
    if (bestOffsets[0] != bestOffsets[1]) {
      conflict = true;
      if (bestCounts[1] < bestCounts[0]) {
        bestOffsets[0] = bestOffsets[1];
      } else {
        bestOffsets[1] = bestOffsets[0];
      }
    }

    // generate the new seed chain and alignment pair if any
    if ((bestCounts[0] <= RESCUE_MAX_SNPS) || (bestCounts[1] <= RESCUE_MAX_SNPS)) {
      const bool reversedRescue = isReversedRescue(anchoredChain);
      rescuedChain.setReverseComplement(reversedRescue);
      if (bestCounts[0] <= RESCUE_MAX_SNPS) {
        const sequences::Seed seed(
            &rescuedRead, 0, std::min((size_t)RESCUE_SEED_LENGTH, rescuedRead.getLength()));
        const int64_t referencePosition =
            rescueIntervalStart +
            (reversedRescue ? (referenceBases.size() - bestOffsets[0] - RESCUE_SEED_LENGTH) : bestOffsets[0]);
        rescuedChain.addSeedPosition(map::SeedPosition(seed, referencePosition, 0), false);
      }
      if (bestCounts[1] <= RESCUE_MAX_SNPS) {
        const sequences::Seed seed(
            &rescuedRead,
            std::max(0, (int)rescuedRead.getLength() - RESCUE_SEED_LENGTH),
            std::min((size_t)RESCUE_SEED_LENGTH, rescuedRead.getLength()));

        const int64_t referencePosition =
            rescueIntervalStart +
            (reversedRescue
                 ? (referenceBases.size() - rescuedRead.getLength() - bestOffsets[1] + modOffset - 1)
                 : rescuedRead.getLength() + bestOffsets[1] - RESCUE_SEED_LENGTH + modOffset - 1);
        rescuedChain.addSeedPosition(map::SeedPosition(seed, referencePosition, 0), false);
      }
      rescuedChain.setPerfect(not conflict);
      rescuedChain.setExtra(anchoredChain.isExtra());
      rescuedChain.setRandomSamplesOnly(anchoredChain.hasOnlyRandomSamples());
      //if (log) std::cerr << "scanning... succeeded" << std::endl;
      DRAGEN_RESCUE_LOG << "RESCUE: "
                        << "first_seed_pos=" << rescuedChain.firstReadBase()
                        << ", last_seed_pos=" << rescuedChain.lastReadBase() << ", first_ref_pos=0x"
                        << std::hex << std::uppercase << std::setw(9) << std::setfill('0')
                        << rescuedChain.firstSeedReferencePosition() - 1 << ", last_ref_pos=0x"
                        << std::setw(9) << rescuedChain.lastSeedReferencePosition() - 1 << std::dec
                        << ", rev_comp=" << rescuedChain.isReverseComplement()
                        << ", perfect_align=" << rescuedChain.isPerfect() << std::endl;
      return true;
    }
  }

  DRAGEN_RESCUE_LOG << "NO_RESCUE" << std::endl;
  return false;
}

int AlignmentRescue::countMismatches(
    const RescueKmer& rescueKmer, const std::vector<unsigned char>::const_iterator position) const
{
  //for (int i = 0; i < rescueKmer.size() ; ++i) std::cerr << " " << (int)rescueKmer[i] << ":" << (int)(*(position + i));
  //std::cerr << std::endl;

  static constexpr unsigned char N     = 0xF;
  int                            count = 0;
  for (size_t i = 0; rescueKmer.size() != i; ++i) {
    const unsigned char b[] = {rescueKmer[i], *(position + i)};
    // Ns and paddings count as mismatches
    if ((b[0] == 0) || (b[1] == 0) || (b[0] == N) || (b[1] == N)) {
      ++count;
    } else {
      bool mismatch = !(b[0] & b[1]);
      count += mismatch;
      if (count > RESCUE_MAX_SNPS) break;
    }
  }
  return count;
}

__m128i AlignmentRescue::loadRescueKmer(AlignmentRescue::RescueKmer rescueKmer) const
{
  // left most base on kmer =  highest address in vector
  std::array<unsigned char, 16> kmer_packed;
  int                           jj     = 0;
  unsigned char                 packed = 0;
  for (size_t ii = 0; rescueKmer.size() != ii; ++ii) {
    packed = packed << 4;
    packed |= rescueKmer[ii];

    if (ii & 1) {
      kmer_packed[15 - jj] = packed;
      jj++;
      packed = 0;
    }
  }
  __m128i kmer_mm = _mm_loadu_si128((__m128i*)kmer_packed.data());

  return kmer_mm;
}

// load kmer from reference into mm128 vector, and set all non ATCG to 0
__m128i AlignmentRescue::loadRefKmer(std::vector<unsigned char>::const_iterator refIter, size_t size) const
{
  std::array<unsigned char, 16> kmer_packed;

  unsigned char packed = 0;
  int           jj     = 0;
  for (size_t ii = 0; size != ii; ++ii) {
    packed             = packed << 4;
    unsigned char pval = *(refIter + ii);
    if (pval == 0xF) pval = 0;
    packed |= pval;
    if (ii & 1) {
      kmer_packed[15 - jj] = packed;
      jj++;
      packed = 0;
    }
  }

  __m128i refkmer_mm = _mm_loadu_si128((__m128i*)kmer_packed.data());

  return refkmer_mm;
}

std::array<AlignmentRescue::RescueKmer, 2> AlignmentRescue::getRescueKmers(
    const Read& rescuedRead, signed modOffset) const
{
  std::array<RescueKmer, 2> rescueKmers;
  std::fill(rescueKmers[0].begin(), rescueKmers[0].end(), 0);
  std::fill(rescueKmers[1].begin(), rescueKmers[1].end(), 0);
  const unsigned length = std::min(rescuedRead.getLength(), rescueKmers[0].size());
  for (unsigned i = 0; i != length; ++i) {
    rescueKmers[0][i] = rescuedRead.getBase4bpb(i);
  }

  const unsigned length2 = std::min(rescuedRead.getLength() - modOffset, rescueKmers[0].size());
  for (unsigned i = 0; i != length2; ++i) {
    rescueKmers[1][i] = rescuedRead.getBase4bpb(i + rescuedRead.getLength() - length2 - modOffset);
  }
  return rescueKmers;
}

bool AlignmentRescue::findRescueChain(
    const Read& /*rescuedRead*/,
    const std::vector<unsigned char>& /*referenceBases*/,
    SeedChain& /*rescuedChain*/) const
{
  return false;
}

#if 0
// Constructor
AlignmentRescue::AlignmentRescue( std::vector<Read*>& reads, std::vector<ChainBuilder*>& chainBuilders)
  : m_pe_orientation(PE_ORIENTATION_FR) 
{
  //  std::cerr << "AlignmentRescue constructor" << std::endl;
  BOOST_ASSERT(reads.size() >= 1);
  if ( reads.size() == 1 || reads[1] == 0) {
    // Single-ended.  Return a list of all of the unfiltered chains from the read.
    InitSingleEnded(reads[0], chainBuilders[0]);
  } else {
    InitPairedEnd(reads, chainBuilders);
  }
}

// AddChain - internal helper to store the read/chain pair in order
//
void AlignmentRescue::AddChain(Read* theread, const SeedChain& chain) 
{
  // chain.set_paired(true);
  m_Chains.push_back(std::make_pair(theread, chain));
}

// InitSingleEnded - initialize our list of chains to align for a single read
//
void AlignmentRescue::InitSingleEnded(Read* theread, ChainBuilder* chains) {
  //  std::cerr << "AlignmentRescue InitSingleEnded" << std::endl;
  for ( auto& chain : *chains) {
    if ( !chain.isFiltered() ) {
      AddChain(theread, chain);
    }
  }
}

#endif

#if 0

// InitPairedEnd - initialize our list of chains with interleaved chains from 
// each of the two reads
//
void AlignmentRescue::InitPairedEnd(
    std::vector<Read*>& reads,
    std::vector<ChainBuilder*> chainBuilders)
{
  //  std::cerr << "AlignmentRescue InitPairedEnd" << std::endl;
  BOOST_ASSERT(reads.size() >= 1);
  BOOST_ASSERT(reads[0] && reads[1]);

  // Do an nxm join of all of the unfiltered chains from the two reads, emitting any
  // naturally-occurring pairs.  
  bool found_a_rescue_proof_pair = false;
  for ( auto& chain1 : *chainBuilders[0]) {
    // if ( it->isFiltered() )
    //   continue;

    bool emitted_chain1 = false;
    for ( auto& chain2 : *chainBuilders[1]) {
      // if ( jt->isFiltered() )
      //   continue;

      // Check if they are within the insert-size limits, and oriented correctly
      if ( AreAPair(chain1, chain2, PE_MIN_INSERT, PE_MAX_INSERT) ) {
        chain1.setFiltered(false);
        chain2.setFiltered(false);

        // avoid duplicate
        if ( !emitted_chain1 ) { 
          emitted_chain1 = true;
          AddChain(reads[0], chain1);
        }
        AddChain(reads[1], chain2);

        // Check if they are actually close enough that we wouldn't need any rescue alignments
        // for the pair
        if ( AreAPair(chain1, chain2, RESCUE_MIN_INSERT, RESCUE_MAX_INSERT) ) {
          found_a_rescue_proof_pair = true;
          chain1.setNeedRescue(false);
          chain2.setNeedRescue(false);
        }
      }
    }
  }

  /*
  if ( found_a_rescue_proof_pair ) 
    std::cerr << " Found a rescue proof pair" << std::endl;
  else
    std::cerr << " Didn't find any rescue proof pairs" << std::endl;
  */
  // Now make a second pass through each chain, generating rescue pairs for any chains that
  // didn't have a really good match from the earlier loop.  In this pass we are ignoring the
  // filtered tag, but we are testing for chain quality based on seed length
  int32_t rescue_seed_len = (found_a_rescue_proof_pair) ? 
    static_cast<int32_t>(RESCUE_CHAIN_LENGTH_WITH_PAIRS) : 
    static_cast<int32_t>(RESCUE_CHAIN_LENGTH_NO_PAIRS);

  for ( int i = 0; i < 2; ++i ) {
    for ( auto chain : *chainBuilders[i]) {
      //TODO: RNA mode and maxRescue
      if ( !chain.needRescue() or !chain.isFiltered() ) {
        continue;
      }
      //TODO: chain flag extra
      else if(chain.getLength() >= rescue_seed_len or chain.hasOnlyRandomSamples()) {
        GenerateRescueAlignments(reads, i, chain);
      }
    }
  }
}


// AreAPair - test whether the two chains are a pair, based on their orientation
// and the distance between their two distal limits.
//
bool AlignmentRescue::AreAPair(
      const SeedChain& chain1,
      const SeedChain& chain2,
      const uint32_t min_insert_size,
      const uint32_t max_insert_size) const
{
  int64_t insert_size = 0;
  if ( m_pe_orientation == PE_ORIENTATION_FF ) {
    if ( chain1.isReverseComplement() != chain2.isReverseComplement() ) 
      return false;
  } else {
    if ( chain1.isReverseComplement() == chain2.isReverseComplement() )
      return false;

    /*  OR is it really this, i.e. first read has to be forward, second RC?
    if ( m_pe_orientation == PE_ORIENTATION_FR ) {
      if ( chain1.is_reverse_complement() || !chain2.is_reverse_complement() )
        return false;
    } else if ( m_pe_orientation == PE_ORIENTATION_RF ) {
      if ( !chain1.is_reverse_complement() || chain2.is_reverse_complement() ) 
        return false;
    } 
    */

    /*
    std::cerr << " are a pair? chain1_start=" << chain1.get_ref_start() 
              << " rc=" << (chain1.is_reverse_complement()?std::string("true"):std::string("false"))
              << ", chain2_end="
              << chain2.get_ref_end() 
              << " rc=" << (chain2.is_reverse_complement()?std::string("true"):std::string("false"))
              << std::endl;
    */
  }
   
  int64_t ref_end = std::max(chain1.lastReferencePosition(), chain2.lastReferencePosition());
  int64_t ref_start = std::min(chain1.firstReferencePosition(), chain2.firstReferencePosition());
  insert_size = ref_end - ref_start;
  return (insert_size >= min_insert_size) && (insert_size <= max_insert_size);
}

// GetRescueOffset - given an #exiting_chain# needing rescue (from the 
// #idx_read_with#'th read), return a pair that includes: the reference
// offset to to the target center of the set of rescue swaths, and 
// whether the rescue should be reverse-complement
//
std::pair<int32_t, bool> AlignmentRescue::GetRescueOffset(
    const SeedChain& existing_chain,
    const int idx_read_with)
{
  // Calculate the direction of the offset to the newly-generated chain
  // assuming that the left read has a chain, and we are generating one for the
  // second (right) one of the pair.  If that's not the case, we'll just negate
  int32_t offset = PE_MEAN_INSERT;
  bool reverse_comp = !existing_chain.isReverseComplement();
  if ( m_pe_orientation == PE_ORIENTATION_FF )
    reverse_comp = !reverse_comp;

  if ( idx_read_with == 1 )
    offset *= -1;

  return std::make_pair(offset, reverse_comp);
}

// GenerateRescue - generate a rescue chain for the specified read (#theread#), 
// with distal terminus at reference position #pos#.
//
void AlignmentRescue::GenerateRescue(
    sequences::Seed seed,
    uint64_t rescue_center,
    const bool reverse_comp,
    ChainBuilder& chainBuilder)
{
  //  std::cerr << "AlignmentRescue GenerateRescue at pos " << pos << std::endl;

  const bool isRandomSample = false;

  uint64_t kmer_pos = rescue_center + 0;//? TODO:need to figure out kmer_pos on reference

  chainBuilder.addSeedPosition(map::SeedPosition(seed, kmer_pos, 0), reverse_comp, isRandomSample);
}

#endif

#if 0
// RescueScan - TODO: scan the read to extract kmers with no more than RESCUE_MAX_SNPS 
//
void AlignmentRescue::RescueScan(
    Read* theread,
    const uint64_t rescue_center,
    const bool reverse_comp)
{
  ChainBuilder chainBuilder;
  //static const uint32_t DEFAULT_PERIOD = 2;
  //static const uint32_t DEFAULT_PATTERN = 0x01;
  //static const uint8_t DEFAULT_FORCE_LAST_N = 1;
  //const auto seedOffsets = sequences::Seed::getSeedOffsets(theread->getLength(), RESCUE_SEED_LENGTH, SEED_PERIOD, SEED_PATTERN, FORCE_LAST_N_SEEDS);

  for (const auto &seedOffset: seedOffsets)
  {
    const sequences::Seed seed(theread, seedOffset, RESCUE_SEED_LENGTH);
    assert(seed.isValid(0)); // getSeedOffset is supposed to produce offsets only for valid non-extended seeds
    GenerateRescue(seed, rescue_center, reverse_comp, chainBuilder);
  }
  
  // should be very few
  for(auto chain: chainBuilder) {
    AddChain(theread, chain);
  }
}
#endif

#if 0
// GenerateRescueAlignments - there is a #chain#, belonging to #read_with#, that
// is missing its pair.  Generate a set of rescue alignments in the #read_without#,
// i.e. in the read that doesn't have the corresponding chain.  
//
void AlignmentRescue::GenerateRescueAlignments(
    std::vector<Read*>& reads,
    const int idx_read_with_chain,
    SeedChain& existing_chain)
{
  /*
  std::cerr << "Generate Rescue Alignments for chain at " << existing_chain.get_ref_start() 
            << ", seed length=" << existing_chain.get_seed_len() << std::endl;
  */
    
  Read* read_with = reads[idx_read_with_chain];
  Read* read_without = reads[(idx_read_with_chain+1)%2];
  std::pair<int32_t, bool> offset_and_rc = GetRescueOffset(existing_chain, idx_read_with_chain);
  const uint64_t rescue_offset = offset_and_rc.first;
  uint64_t rescue_center = existing_chain.firstReferencePosition() + rescue_offset;
  const bool reverse_comp = offset_and_rc.second;

  AddChain(read_with, existing_chain);
  RescueScan(read_without, rescue_center, reverse_comp);
  if ( m_pe_orientation == PE_ORIENTATION_FF ) {
    rescue_center = existing_chain.firstReferencePosition() - rescue_offset;
    RescueScan(read_without, rescue_center, reverse_comp);
  }
}
#endif

}  // namespace align
}  // namespace dragenos
