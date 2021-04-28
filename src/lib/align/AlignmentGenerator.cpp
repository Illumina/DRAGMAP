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

#include "align/AlignmentGenerator.hpp"
#include "align/CalculateRefStartEnd.hpp"
#include "common/DragenLogger.hpp"

namespace dragenos {
namespace align {

void AlignmentGenerator::generateAlignments(
    const Read& read, const map::ChainBuilder& chainBuilder, Alignments& alignments)
{
  for (const auto& seedChain : chainBuilder) {
    auto& alignment = alignments.addAlignment();
    generateAlignment(0, read, seedChain, alignment);
  }
}

bool AlignmentGenerator::generateAlignment(
    const ScoreType alnMinScore, const Read& read, const map::SeedChain& seedChain, Alignment& alignment)
{
  Database database;
  // TODO: create the database as a vector of unsigned char, 1 base per unsigned char, encoded on 2 bits
  //const auto databaseBegin = referenceDir_.getReference() + seedChain.firstReferencePosition();
  //const auto databaseEnd = referenceDir_.getReference() + seedChain.lastReferencePosition() + 1;
  //    const auto beginPosition = seedChain.firstReferencePosition();
  //    const auto endPosition =  seedChain.lastReferencePosition() + 1;
  const auto refStartEnd   = calculateRefStartEnd(read, seedChain);
  auto       beginPosition = refStartEnd.first;
  auto       endPosition   = refStartEnd.second + 1;

  const reference::HashtableConfig& hashtableConfig = referenceDir_.getHashtableConfig();
  // endPosition should be bound by sequence end
  if (not hashtableConfig.beyondLastCfgSequence(refStartEnd.first)) {
    auto refCoords = hashtableConfig.convertToReferenceCoordinates(refStartEnd.first);
    const reference::HashtableConfig::Sequence& seq      = hashtableConfig.getSequences().at(refCoords.first);
    const auto                                  posRange = hashtableConfig.getPositionRange(seq);
    endPosition                                          = std::min<int64_t>(endPosition, posRange.second);
  } else {
    // RP: from FZ:
    //    The purpose is to flag it as ineligible to be as final alignment even it is aligned but to a
    //    implicit decoy sequence
    alignment.setIneligibilityStatus(true);
    // RP: I was tempted to simply return false from here but according to FZ we still need to perform the s-w
    // alignment.
  }
  //  std::cerr << "beginPosition=" << beginPosition << " endPosition=" << endPosition << std::endl;

  //  assert(endPosition > beginPosition);
  database.clear();
  if (seedChain.isReverseComplement()) {
    referenceDir_.getReferenceSequence().getRcBases(beginPosition, endPosition, std::back_inserter(database));
  } else {
    referenceDir_.getReferenceSequence().getBases(beginPosition, endPosition, std::back_inserter(database));
  }

  // align the read for the current seedChain
  // TODO: put this configuration parameter in the right location
  //  std::cerr << seedChain << std::endl;
  const size_t forcedDiagonalMotion = std::max<int>(
      1,
      (seedChain.isReverseComplement() ? (read.getLength() - seedChain.lastReadBase() - 1)
                                       : seedChain.firstReadBase()) *
          2);  //10;//1;
  static constexpr size_t forcedHorizontalMotion = smithWaterman_.width;
  // initialize the query from the base and the orientation of the seedChain
  const auto& query = read.getBases();
  int         move  = 0;
  std::string operations;
  FlagType    flags = !read.getPosition() ? Alignment::FIRST_IN_TEMPLATE : Alignment::LAST_IN_TEMPLATE;
  {
    const ScoreType score = smithWaterman_.align(
        query.data(),
        query.data() + query.size(),
        database.data(),
        database.data() + database.size(),
        forcedDiagonalMotion,
        forcedHorizontalMotion,
        // dragen right-shifts indels for reverse-complement alignments
        seedChain.isReverseComplement(),
        operations);

#ifdef TRACE_SMITH_WATERMAN
    std::cerr << "[SMITH-WATERMAN]\tdb\t" << database << "\n[SMITH-WATERMAN]\tops\t" << operations
              << "\n[SMITH-WATERMAN]\tquery\t" << Query(query.begin(), query.end())
              << "\trc:" << seedChain.isReverseComplement() << "\tscore:" << score << "("
              << smithWaterman_.getMaxScore() << ")"
              << "\talnMinScore:" << alnMinScore << std::endl;
#endif

    if (alnMinScore > score) {
      return false;
    }

    alignment.setScore(score);
    //    alignment.setPotentialScore(score); // no more improvement to expect
    alignment.setSmithWatermanDone(true);

    move = int(alignment.setCigarOperations(
               operations,
               database.begin(),
               database.end(),
               query.begin(),
               query.end(),
               seedChain.isReverseComplement())) +
           beginPosition - seedChain.firstReferencePosition();  // - SmithWaterman::SW_CELLS;

    //Alignment::FlagType flags = !read.getPosition() ? Alignment::FIRST_IN_TEMPLATE : Alignment::LAST_IN_TEMPLATE;
    flags |= seedChain.isReverseComplement() ? Alignment::REVERSE_COMPLEMENT : 0;
  }
  // get the reference name from the ReferenceDir and the first diagonal
  const size_t referenceOffset      = seedChain.firstReferencePosition();
  auto         referenceCoordinates = hashtableConfig.convertToReferenceCoordinates(referenceOffset + move);

  if (0 > referenceCoordinates.second) {
    if (std::abs(referenceCoordinates.second) < read.getLength()) {
      // oops! the sequence starts before reference and the read is long enough to reach the next contig. Redo
      // CIGAR.
      move = int(alignment.setCigarOperations(
                 operations,
                 database.begin(),
                 database.end(),
                 query.begin(),
                 query.end(),
                 seedChain.isReverseComplement(),
                 -referenceCoordinates.second)) +
             beginPosition - seedChain.firstReferencePosition();  // - SmithWaterman::SW_CELLS;
      // get the reference name from the ReferenceDir and the first diagonal
      const size_t referenceOffset = seedChain.firstReferencePosition();
      referenceCoordinates         = hashtableConfig.convertToReferenceCoordinates(referenceOffset + move);
    } else {
      // the whole read falls into holes of reference sequences listed in hash_table.cfg)
      // TODO:consider if there are second layer of holes(padding bytes) within reference.bin and also need to
      // redo CIGAR.
      referenceCoordinates.second = 0;
      alignment.setIneligibilityStatus(true);
    }
  }

  alignment.resetFlags(flags);
  //    alignment.setReferenceName(std::string(hashtableConfig.getSequenceNames()[referenceCoordinates.first]));
  alignment.setReference(referenceCoordinates.first);
  alignment.setPosition(referenceCoordinates.second);

  // debug info
#ifdef TRACE_ALIGNMENTS
  std::cerr << "[ALIGNMENT]\tsmith-waterman\t" << alignment << "\t" << seedChain
            << "\tcontig:" << std::string(hashtableConfig.getSequenceNames()[referenceCoordinates.first])
            << "(" << referenceCoordinates.first << "):" << referenceCoordinates.second << std::endl;
#endif

  return true;
}

}  // namespace align
}  // namespace dragenos
