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
    const Read& read, const map::ChainBuilder& chainBuilder, Alignments& alignments, const int readIdx)
{
  for (const auto& seedChain : chainBuilder) {
    auto& alignment = alignments.addAlignment();
    generateAlignment(0, read, seedChain, alignment, readIdx);
  }
}

void AlignmentGenerator::updateFetchChain(const Read& read, map::SeedChain& seedChain, Alignment& alignment)
{
  // only update chain if the alignment is valid
  const size_t referenceOffset = seedChain.firstReferencePosition();
  auto         refCoords       = htConfig_.convertToReferenceCoordinates(referenceOffset);
  int          move            = alignment.getPosition() - refCoords.second;
  if (alignment.isPerfect()) {
    if (seedChain.isReverseComplement()) {
      seedChain.clear();
      seedChain.addSeedPosition(
          map::SeedPosition(
              sequences::Seed(&read, read.getLength() - alignment.getCigar().countStartClips() - 1, 1),
              referenceOffset + move,
              0),
          false);
      //        std::cerr << "between chain update:" << seedChain << " " << alignment << std::endl;
      seedChain.addSeedPosition(
          map::SeedPosition(
              sequences::Seed(&read, alignment.getCigar().countEndClips(), 1),
              referenceOffset + move + alignment.getRefLen() - 1,
              0),
          false);
      seedChain.setReverseComplement(true);
    } else {
      seedChain.clear();
      seedChain.addSeedPosition(
          map::SeedPosition(
              sequences::Seed(&read, alignment.getCigar().countStartClips(), 1), referenceOffset + move, 0),
          false);
      seedChain.addSeedPosition(
          map::SeedPosition(
              sequences::Seed(&read, read.getLength() - alignment.getCigar().countEndClips() - 1, 1),
              referenceOffset + move + alignment.getRefLen() - 1,
              0),
          false);
    }
    seedChain.setPerfect(alignment.isPerfect());
  }
}

void updateIneligibility(const reference::HashtableConfig& hashtableConfig, Alignment& alignment)
{
  const auto& sequences = hashtableConfig.getSequences();
  const auto& sequence  = sequences.at(alignment.getReference());
  if (alignment.getPosition() + alignment.getCigar().getReferenceLength() <= 0 ||
      int(sequence.seqLen) <= alignment.getPosition()) {
    alignment.setIneligibilityStatus(true);
  }
}

bool AlignmentGenerator::generateAlignment(
    const ScoreType alnMinScore,
    const Read&     read,
    map::SeedChain  seedChain,
    Alignment&      alignment,
    const int       readIdx)
{
  if (seedChain.isFiltered()) {
    return false;
  }

  updateFetchChain(read, seedChain, alignment);

  DRAGEN_S_W_FETCH_LOG << seedChain << std::endl;

  Database database;
  // TODO: create the database as a vector of unsigned char, 1 base per unsigned char, encoded on 2 bits
  //const auto databaseBegin = referenceDir_.getReference() + seedChain.firstReferencePosition();
  //const auto databaseEnd = referenceDir_.getReference() + seedChain.lastReferencePosition() + 1;
  //    const auto beginPosition = seedChain.firstReferencePosition();
  //    const auto endPosition =  seedChain.lastReferencePosition() + 1;
  const auto refStartEnd   = calculateRefStartEnd(read, seedChain);
  auto       beginPosition = refStartEnd.first;
  auto       endPosition   = refStartEnd.second + 1;

  // endPosition should be bound by sequence end
  if (not htConfig_.beyondLastCfgSequence(refStartEnd.first)) {
    auto refCoords = htConfig_.convertToReferenceCoordinates(refStartEnd.first);
    const reference::HashtableConfig::Sequence& seq      = htConfig_.getSequences().at(refCoords.first);
    const auto                                  posRange = htConfig_.getPositionRange(seq);
    endPosition                                          = std::min<int64_t>(endPosition, posRange.second);
    if (endPosition <= beginPosition) return false;
  } else {
    // RP: from FZ:
    //    The purpose is to flag it as ineligible to be as final alignment even it is aligned but to a
    //    implicit decoy sequence I believe completely remove the else section should be the same logic
    //    because the seed would still be the same as in the “holes” on the reference and hence it will be
    //    flagged in the 0>referenceCordinate.second logic branch To explicitly flag it here would be easier
    //    to read the code
    //
    // Note: this is not the best place to update eligibility. It is done before further down.
    // alignment.setIneligibilityStatus(true);
    // RP: I was tempted to simply return false from here but according to FZ we still need to perform the s-w
    // alignment.
  }
  //  std::cerr << "beginPosition=" << beginPosition << " endPosition=" << endPosition << std::endl;

  //  assert(endPosition > beginPosition);
  database.clear();
  database.reserve(endPosition - beginPosition + 1);
  if (seedChain.isReverseComplement()) {
    refSeq_.getRcBases(beginPosition, endPosition, database);
  } else {
    refSeq_.getBases(beginPosition, endPosition, database);
  }

  // align the read for the current seedChain
  // TODO: put this configuration parameter in the right location
  //  std::cerr << seedChain << std::endl;
  const size_t forcedDiagonalMotion = std::max<int>(
      1,
      (seedChain.isReverseComplement() ? (read.getLength() - seedChain.lastReadBase() - 1)
                                       : seedChain.firstReadBase()));  //10;//1;
  static const size_t forcedHorizontalMotion = smithWaterman_.width;
  // initialize the query from the base and the orientation of the seedChain
  const auto& query = read.getBases();
  int         move  = 0;
  std::string operations;
  FlagType    flags = !read.getPosition() ? Alignment::FIRST_IN_TEMPLATE : Alignment::LAST_IN_TEMPLATE;
  {
    ScoreType scoreSW;

    if (vectorizedSW_ && query.size() > 30) {
      scoreSW = vectorSmithWaterman_.align(
          query.data(),
          query.data() + query.size(),
          database.data(),
          database.data() + database.size(),
          seedChain.isReverseComplement(),
          operations,
          readIdx);

    } else {
      scoreSW = smithWaterman_.align(
          query.data(),
          query.data() + query.size(),
          database.data(),
          database.data() + database.size(),
          forcedDiagonalMotion,
          forcedHorizontalMotion,
          // dragen right-shifts indels for reverse-complement alignments
          seedChain.isReverseComplement(),
          operations);
    }
    const ScoreType score = scoreSW;

#ifdef TRACE_SMITH_WATERMAN
    std::cerr << "[SMITH-WATERMAN]\tdb\t" << database << "\n[SMITH-WATERMAN]\tops\t" << operations
              << "\n[SMITH-WATERMAN]\tquery\t"
              << (seedChain.isReverseComplement() ? Query(query.rbegin(), query.rend())
                                                  : Query(query.begin(), query.end()))
              << "\trc:" << seedChain.isReverseComplement() << "\tscore:" << score << "("
              << smithWaterman_.getMaxScore() << ")"
              << "\talnMinScore:" << alnMinScore << "\tfdm:" << forcedDiagonalMotion
              << "\tfhm:" << forcedHorizontalMotion << std::endl;
#endif

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
    flags |= seedChain.isReverseComplement() ? Alignment::REVERSE_COMPLEMENT : Alignment::NONE;
  }
  // get the reference name from the ReferenceDir and the first diagonal
  const size_t referenceOffset      = seedChain.firstReferencePosition();
  auto         referenceCoordinates = htConfig_.convertToReferenceCoordinates(referenceOffset + move);

  if (0 > referenceCoordinates.second) {
    if (std::size_t(std::abs(referenceCoordinates.second)) < read.getLength()) {
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
      referenceCoordinates         = htConfig_.convertToReferenceCoordinates(referenceOffset + move);
    } else {
      // the whole read falls into holes of reference sequences listed in hash_table.cfg)
      // TODO:consider if there are second layer of holes(padding bytes) within reference.bin and also need to
      // redo CIGAR.
      referenceCoordinates.second = 0;
    }
  }

  alignment.resetFlags(flags);
  //    alignment.setReferenceName(std::string(hashtableConfig.getSequenceNames()[referenceCoordinates.first]));
  alignment.setReference(referenceCoordinates.first);
  alignment.setPosition(referenceCoordinates.second);

  updateIneligibility(htConfig_, alignment);

  if (alnMinScore > alignment.getScore()) {
    alignment.setUnmapped();
  }

  // debug info
#ifdef TRACE_ALIGNMENTS
  std::cerr << "[ALIGNMENT]\tsmith-waterman\t" << alignment << "\t" << seedChain
            << "\tcontig:" << std::string(htConfig_.getSequenceNames()[referenceCoordinates.first]) << "("
            << referenceCoordinates.first << "):" << referenceCoordinates.second << std::endl;
#endif

  return true;
}

}  // namespace align
}  // namespace dragenos
