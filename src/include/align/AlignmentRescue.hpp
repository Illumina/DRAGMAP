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

#ifndef ALIGN_ALIGNMENT_RESCUE_HPP
#define ALIGN_ALIGNMENT_RESCUE_HPP

#include <array>
#include <deque>
#include "align/Alignment.hpp"
#include "align/InsertSizeParameters.hpp"
#include "map/ChainBuilder.hpp"
#include "reference/ReferenceSequence.hpp"
#include "sequences/Read.hpp"

namespace dragenos {
namespace align {

/**
 ** \brief component encapsulating the method to rescue mapping positions missed by the standard search in the
 *hash table
 **
 ** Triggering the rescue:
 ** - rescue scans are never triggered if
 **   -- the seed chain is properly paired
 **   -- the seed chain is both low-coverage and flagged as filtered_out
 **   -- if the total number of seed chains not flagged filtered out exceeds max_rescues
 **   -- RNA mode
 ** - otherwise, for a seed chain spanning chain_len bases over the anchored read trigger rescue if
 **   -- any_pair_match = 0 and chain_len >= resc_nopair_len
 **   -- any_pair_match = 1 and chain_len >= resc_ifpair_len
 **   -- seed chain flagged as "extra"
 **   -- chain is flagged as "sample" and rescue_hifreq=1
 **
 ** Rescue interval to scan:
 ** - basically the reference interval implied by the minimum and maximum insert sizes configured for rescue
 *scans
 **   -- Let P be the outer edge of the anchored mate
 **   -- rescue_min = P +/- rescue_min_ins (sign depends on direction in template)
 **   -- rescue_max = P +/- rescue_max_ins (sign depends on direction in template)
 ** - orientation infered from expected paired end orientation and the orientation of anchored mate
 ** - rescued read is scanned forward - reference is reversed if needed
 **
 ** Method:
 ** - extract kmers of length rescue_kmer_len from the beginning and end of the mate to rescue
 ** - scan the reference rescue interval for matches to each of these K-mers with no more than rescue_max_snps
 *mismatches
 **   -- N on reference or read counts as a mismatch
 **   -- Non N Multi-base match as appropriate
 **   -- for reads shorter than read length, missing bases count as mismatches
 ** - select the best match position for each K-mer, breaking ties by smaller reference position and setting
 *the conflict to true if either K-mer match to more than one position
 ** - if both K-mers have matching positions but on different diagonals, keep the best-scoring and set the
 *conflict flag to true
 ** - fabricate a seed chain that represents the succesful rescue if any
 **   -- if both K-mer match to same diagonal, the seedchain from tha start of the start K-mer to the end of
 *the end K-mer
 **   -- else the seedchain spans from the begin to end of the matching K-mer
 ** - flag seed chain as perfect align if conflict is false (trigger Smith-Waterman if conflict is true)
 ** - inject the seed chain into the source seed chain for the other mate, flagged as paired with it
 **/
class AlignmentRescue {
public:
  typedef sequences::Read                       Read;
  typedef map::SeedChain                        SeedChain;
  typedef map::ChainBuilder                     ChainBuilder;
  typedef std::pair<Read*, SeedChain>           Read_and_chain;
  typedef std::deque<Read_and_chain>            Read_and_chain_c;
  typedef Read_and_chain_c::iterator            iterator;
  typedef align::InsertSizeParameters           InsertSizeParameters;
  typedef InsertSizeParameters::Orientation     Orientation;
  typedef decltype(Orientation::pe_orient_fr_c) PairedEndOrientation;

private:
  /// TODO: get this from command line or configuration options
  static constexpr int PE_MEAN_INSERT = 318;
  //static constexpr int PE_MIN_INSERT = 41;
  //static constexpr int PE_MAX_INSERT = 596;
  //static constexpr int RESCUE_MIN_INSERT = 41;
  //static constexpr int RESCUE_MAX_INSERT = 596;
  static constexpr int RESCUE_CHAIN_LENGTH_WITH_PAIRS = 48;
  static constexpr int RESCUE_CHAIN_LENGTH_NO_PAIRS   = 0;
  static constexpr int RESCUE_SEED_LENGTH             = 32;
  static constexpr int RESCUE_MAX_SNPS                = 7;

  static constexpr auto pe_orient_fr_c = Orientation::pe_orient_fr_c;
  static constexpr auto pe_orient_rf_c = Orientation::pe_orient_rf_c;
  static constexpr auto pe_orient_ff_c = Orientation::pe_orient_ff_c;
  static constexpr auto pe_orient_rr_c = Orientation::pe_orient_rr_c;

  //enum PairedEndOrientation {
  //  PE_ORIENTATION_FR = 0,
  //  PE_ORIENTATION_RF = 1,
  //  PE_ORIENTATION_FF = 2
  //};

public:
  /// initialize alignment rescue with all the relevant control parameters
  AlignmentRescue(
      int         pe_min_insert,
      int         pe_max_insert,
      Orientation pe_orientation,
      int         resc_nopair_len = RESCUE_CHAIN_LENGTH_NO_PAIRS,
      int         resc_ifpair_len = RESCUE_CHAIN_LENGTH_WITH_PAIRS)
    : pe_min_insert_(pe_min_insert),
      pe_max_insert_(pe_max_insert),
      pe_orientation_(pe_orientation),
      resc_nopair_len_(resc_nopair_len),
      resc_ifpair_len_(resc_ifpair_len)
  {
  }

  /**
   ** \brief check if the given read and seed chain trigger the rescue
   **
   ** This must be called only if the seed chain qualifies for rescue. These are the conditions for which
   *rescue should never be triggered:
   ** - the seed chain is properly paired
   ** - the seed chain is both low-coverage and flagged as filtered_out
   ** - if the total number of seed chains not flagged filtered out exceeds max_rescues
   ** - RNA mode
   **
   ** for a seed chain spanning over chin_len bases that qualifies for potential rescue, the triggering
   *conditions are:
   **   -- any_pair_match = 0 and chain_len >= resc_nopair_len
   **   -- any_pair_match = 1 and chain_len >= resc_ifpair_len
   **   -- seed chain flagged as "extra" (TODO - not implemented)
   **   -- chain is flagged as "sample" and rescue_hifreq=1 (legacy - not implemented)
   **
   ** \param read the anchored read
   ** \param anchoredChain a seed chain for the anchor read that qualifies for rescue
   ** \param any_pair_match flag indicating if there is a correctly matching pair for at least one other seed
   *chain
   ** \return true if the rescue scan should be executed
   **/
  bool triggeredBy(const SeedChain& anchoredChain, const bool any_pair_match) const;

  /**
   ** \brief scan the reference interval for a  suitable seed chain if any
   **/
  bool scan(
      const Read&                         rescuedRead,
      const SeedChain&                    anchoredChain,
      const reference::ReferenceSequence& reference,
      map::SeedChain&                     rescuedChain) const;

  // private:
  /**
   ** \brief retrieve the appropriate sequence of bases from the reference
   **
   ** \param anchoredChain the seed chain for the anchored read
   ** \param rescuedReadLength the length of the read to rescue
   ** \param reference the reference to get the data from
   ** \param referenceBases the buffer to store the reference (possibly reverse complemented
   ** \param second optional parameter for FF and RR orientations where two rescue regions must be scanned
   ** \return the reference offset of the leftmost base of the interval
   **/
  uint32_t getReferenceInterval(
      const SeedChain&                    anchoredChain,
      const int                           rescuedReadLength,
      const reference::ReferenceSequence& reference,
      std::vector<unsigned char>&         referenceBases,
      const bool                          second = false) const;

  /**
   ** \brief find a rescue chain for the given rean in the given reference seauence if any
   **/
  bool findRescueChain(
      const Read&                       rescuedRead,
      const std::vector<unsigned char>& referenceBases,
      SeedChain&                        rescuedChain) const;

  typedef std::array<unsigned char, 32> RescueKmer;
  /**
   ** \brief generate the two rescue kmers for the read to rescue
   **/
  std::array<RescueKmer, 2> getRescueKmers(const Read& rescuedRead) const;

  /**
   ** \brief count mismatches between the kmer and the reference
   **/
  int countMismatches(
      const RescueKmer& rescueKmer, const std::vector<unsigned char>::const_iterator position) const;

  /**
   **
   **/
  bool isReversedRescue(const SeedChain& anchoredChain) const
  {
    return (pe_orientation_ == pe_orient_fr_c)
               ? (!anchoredChain.isReverseComplement())
               : (pe_orientation_ == pe_orient_rf_c)
                     ? (!anchoredChain.isReverseComplement())
                     : (pe_orientation_ == pe_orient_ff_c) ? (false)
                                                           : /* (pe_orientation_==pe_orient_rr_c) ? */ (true);
  }

#if 0

  // input interface
  //AlignmentRescue( 
  //    std::vector<Read*>& reads, 
  //    std::vector<ChainBuilder*>& chainBuilders)
  //  ;

  // output interface
  template <typename OutputIter>
  OutputIter GetAllChains(OutputIter o) {
    for ( Read_and_chain_c::iterator it = m_Chains.begin(); it != m_Chains.end(); ++it ) {
      o++ = (*it);
    }
  }

public:
  iterator begin() {
    return m_Chains.begin();
  }

  iterator end() {
    return m_Chains.end();
  }

private:
  void InitSingleEnded(Read* read,
      ChainBuilder* chains)
    ;

  void InitPairedEnd(
      std::vector<Read*>& reads,
      std::vector<ChainBuilder*> chainBuilders)
    ;

  bool AreAPair(
      const SeedChain& chain1,
      const SeedChain& chain2,
      const uint32_t min_insert_size,
      const uint32_t max_insert_size) const
    ;

  void AddChain(
      Read* read, 
      const SeedChain& chain)
    ;

  std::pair<int32_t, bool> GetRescueOffset(
      const SeedChain& existing_chain,
      const int idx_read_with)
    ;

  void GenerateRescue(
    sequences::Seed seed,
    uint64_t rescue_center,
    const bool reverse_comp,
    ChainBuilder& chainBuilder)
    ;

  //void RescueScan(
  //    Read* theread,
  //    const uint64_t rescue_center,
  //    const bool reverse_comp)
  //  ;

  //void GenerateRescueAlignments(
  //    std::vector<Read*>& reads,
  //    const int idx_read_with_chain,
  //    SeedChain& existing_chain)
  //  ;

#endif

private:
  Read_and_chain_c m_Chains;  // the seed-chains, along with the corresponding read, in order
  int              pe_min_insert_;
  int              pe_max_insert_;
  Orientation      pe_orientation_;
  int              resc_nopair_len_;
  int              resc_ifpair_len_;
};

}  // namespace align
}  // namespace dragenos

#endif  // #ifndef ALIGN_ALIGNMENT_RESCUE_HPP
