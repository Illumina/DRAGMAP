
#include <boost/filesystem.hpp>
#include <sstream>

#include "gtest/gtest.h"

#include "align/PairBuilder.hpp"
#include "align/Pairs.hpp"
#include "reference/ReferenceDir.hpp"

static char emptySpace[1024] = {};
class ReferenceDirDummy : public dragenos::reference::ReferenceDir {
  dragenos::reference::HashtableConfig   hashtableConfig_;
  dragenos::reference::ReferenceSequence referenceSequence_;

public:
  ReferenceDirDummy() : hashtableConfig_(emptySpace, sizeof(emptySpace)) {}
  virtual const dragenos::reference::HashtableConfig& getHashtableConfig() const { return hashtableConfig_; };
  virtual const uint64_t*                             getHashtableData() const { return nullptr; }
  virtual const dragenos::reference::ReferenceSequence& getReferenceSequence() const
  {
    return referenceSequence_;
  }

  virtual const uint64_t* getExtendTableData() const { return 0; }
};

using namespace dragenos;
typedef sequences::Read             Read;
typedef sequences::ReadPair         ReadPair;
typedef Read::Name                  Name;
typedef Read::Bases                 Bases;
typedef Read::Qualities             Qualities;
typedef align::Alignment            Alignment;
typedef align::AlignmentPair        AlignmentPair;
typedef align::AlignmentPairs       AlignmentPairs;
typedef align::InsertSizeParameters InsertSizeParameters;
typedef align::Alignments           Alignments;
typedef map::SeedPosition           SeedPosition;
typedef sequences::Seed             Seed;

TEST(PairBuilder, updateMapq)
{
  const short                             match    = 1;
  const short                             mismatch = -1;
  const dragenos::align::SimilarityScores similarityScores(match, mismatch);
  const align::PairBuilder                pairBuilder(similarityScores, 19, 80, 25, 0, 0, 0, false, 50, 0);

  const unsigned char reference[] =
      "TCCATCGAGATGGACGCCGTTGGCGCTCTCCGTCTTTCTCCATTGCGTCGTGGCCTTGCTATTGACTCTACTGTAGACATTTTTACTTTTTATGTCCCTCTTATG";

  static const unsigned SEQUENCE_FLANKS = 20;
  const std::string     name            = "blah";
  const std::string     bases(
      (const char*)reference + SEQUENCE_FLANKS, (const char*)reference + sizeof(reference) - SEQUENCE_FLANKS);
  const std::string qualities =
      "#AAAAEEEEEEEEEEEEEEEEEEEEEEE6EEEEEEEEEEEEEEEEEAEAEEEEEEEEEEEEAEEE/EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE";

  uint64_t id       = 17;

  ReadPair readPair;
  readPair[0].init(
      Name(name.begin(), name.end()),
      Bases(bases.begin(), bases.end()),
      Qualities(qualities.begin() + SEQUENCE_FLANKS, qualities.begin() + SEQUENCE_FLANKS + bases.length()),
      id,
      0);

  readPair[1].init(
      Name(name.begin(), name.end()),
      Bases(bases.begin(), bases.end()),
      Qualities(qualities.begin() + SEQUENCE_FLANKS, qualities.begin() + SEQUENCE_FLANKS + bases.length()),
      id,
      1);

  AlignmentPairs            alignments;
  align::UnpairedAlignments unpairedAlignments;
  unpairedAlignments[0].addAlignment();
  unpairedAlignments[1].addAlignment();
  unpairedAlignments[0].at(0).setPosition(123);
  unpairedAlignments[1].at(0).setPosition(321);
  unpairedAlignments[0].at(0).setCigarOperations(std::string(bases.length()/3*2, 'M') + std::string(bases.length() - bases.length()/3*2, 'S'));
  unpairedAlignments[1].at(0).setCigarOperations(std::string(bases.length() - bases.length()/3*2, 'S') + std::string(bases.length()/3*2, 'M'));
  alignments.push_back(AlignmentPair(unpairedAlignments[0].at(0), unpairedAlignments[1].at(0)));
  alignments.push_back(AlignmentPair(unpairedAlignments[1].at(0), unpairedAlignments[0].at(0)));
  alignments[0].setScore(113);
  alignments[1].setScore(0);
  pairBuilder.pickBest(readPair, alignments, unpairedAlignments);
  // unique alignment pair. Both ends get mapq 60
  ASSERT_EQ(804, alignments.back()[0].getMapq());
  ASSERT_EQ(804, alignments.back()[1].getMapq());

  unpairedAlignments[0].addAlignment();
  unpairedAlignments[1].addAlignment();
  unpairedAlignments[0].at(1).setFlags(align::AlignmentHeader::REVERSE_COMPLEMENT);
  unpairedAlignments[0].at(1).setCigarOperations(std::string(bases.length()/3*2, 'M') + std::string(bases.length() - bases.length()/3*2, 'S'));
  unpairedAlignments[1].at(1).setFlags(align::AlignmentHeader::REVERSE_COMPLEMENT);
  unpairedAlignments[1].at(1).setCigarOperations(std::string(bases.length() - bases.length()/3*2, 'S') + std::string(bases.length()/3*2, 'M'));
  alignments.push_back(AlignmentPair(unpairedAlignments[0].at(1), unpairedAlignments[1].at(1)));
  alignments.back().setScore(108);
  pairBuilder.pickBest(readPair, alignments, unpairedAlignments);
  // one end gets suboptimal alignment
  ASSERT_EQ(48, alignments.front()[0].getMapq());
  ASSERT_EQ(48, alignments.front()[1].getMapq());
}

unsigned char encodeBase(const char base)
{
  switch (base) {
  case 'A':
    return 1;
    break;
  case 'C':
    return 2;
    break;
  case 'G':
    return 4;
    break;
  case 'T':
    return 8;
    break;
  default:
    throw std::invalid_argument(std::string("unknown base: ") + base);
  }
}

void text2Bases(char* begin, char* end)
{
  std::transform(begin, end, begin, [](const char c){return encodeBase(c);});
}

TEST(PairBuilder, pickBest1)
{
  const short                             match    = 1;
  const short                             mismatch = -1;
  const dragenos::align::SimilarityScores similarityScores(match, mismatch);
  const ReferenceDirDummy                 referenceDir;
  const align::PairBuilder                pairBuilder(similarityScores, 19, 80, 25, 0, 0, 0, false, 50, 0);

  char reference[] =
      "TCCATCGAGATGGACGCCGTTGGCGCTCTCCGTCTTTCTCCATTGCGTCGTGGCCTTGCTATTGACTCTACTGTAGACATTTTTACTTTTTATGTCCCTCTTATG";

  text2Bases(reference, reference + sizeof(reference) - 1);

  static const unsigned SEQUENCE_FLANKS = 10;
  const std::string     name            = "blah";
  const std::string     bases(
      (const char*)reference + SEQUENCE_FLANKS, (const char*)reference + sizeof(reference) - 1 - SEQUENCE_FLANKS * 2);
  const std::string qualities =
      "#AAAAEEEEEEEEEEEEEEEEEEEEEEE6EEEEEEEEEEEEEEEEEAEAEEEEEEEEEEEEAEEE/EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE";

  uint64_t id       = 17;
  unsigned position = 0;

  ReadPair readPair;
  readPair[0].init(
      Name(name.begin(), name.end()),
      Bases(bases.begin(), bases.end()),
      Qualities(qualities.begin() + SEQUENCE_FLANKS, qualities.begin() + SEQUENCE_FLANKS + bases.length()),
      id,
      position);

  readPair[1].init(
      Name(name.begin(), name.end()),
      Bases(bases.begin(), bases.end()),
      Qualities(qualities.begin() + SEQUENCE_FLANKS, qualities.begin() + SEQUENCE_FLANKS + bases.length()),
      id,
      position);

//  std::cerr << readPair[1] << std::endl;

  AlignmentPairs alignments;

  const InsertSizeParameters insertSizeParameters = {
      332, 486, 410, 332, 486, 14.0708, InsertSizeParameters::Orientation::pe_orient_fr_c};
  std::array<map::ChainBuilder, 2> chainBuilders{map::ChainBuilder(4.0), map::ChainBuilder(4.0)};

  chainBuilders[0].addSeedPosition(SeedPosition(Seed(&readPair[0], 0, 17), 0x2808c, 0), false, false);
  chainBuilders[0].addSeedPosition(SeedPosition(Seed(&readPair[0], 49, 17), 0x280bd, 0), false, false);

  chainBuilders[1].addSeedPosition(SeedPosition(Seed(&readPair[1], 0, 17), 0x2816b, 0), true, false);
  chainBuilders[1].addSeedPosition(SeedPosition(Seed(&readPair[1], 49, 17), 0x2813a, 0), true, false);

  chainBuilders[1].addSeedPosition(SeedPosition(Seed(&readPair[1], 30, 17), 0x281d9, 0), true, false);
  chainBuilders[1].addSeedPosition(SeedPosition(Seed(&readPair[1], 49, 17), 164294, 0), true, false);

  chainBuilders[1].addSeedPosition(SeedPosition(Seed(&readPair[1], 0, 17), 164063, 0), true, false);

  std::array<Alignments, 2> unpaired;
  Alignment                 a11, a21, a22, a23;
  a11.setScore(66);
  a11.setCigarOperations(std::to_string(bases.length()) + "M");
  a21.setScore(66);
  a21.setCigarOperations(std::to_string(bases.length()) + "M");
  a22.setScore(56);
  a22.setCigarOperations(std::to_string(bases.length()) + "M");
  a22.setFlags(align::AlignmentHeader::REVERSE_COMPLEMENT);
  a23.setScore(17);
  a23.setCigarOperations(std::to_string(bases.length()) + "M");
  unpaired[0].append(a11);
  unpaired[1].append(a21);
  unpaired[1].append(a22);
  unpaired[1].append(a23);

  for (unsigned i0 = 0; unpaired[0].size() > i0; ++i0) {
    for (unsigned i1 = 0; unpaired[1].size() > i1; ++i1) {
      alignments.insert(alignments.end(), AlignmentPair(unpaired[0].at(i0), unpaired[1].at(i1)));
      const bool isPair =
          pairMatch(insertSizeParameters, readPair, chainBuilders[0].at(i0), chainBuilders[1].at(i1));
      int       insert_len  = 0;
      int       insert_diff = 0;
      const int pair_pen = pairBuilder.computePairPenalty(
          insertSizeParameters, readPair, &unpaired[0].at(i0), &unpaired[1].at(i1), isPair, insert_len, insert_diff);
      alignments.back().setScore(
          alignments.back()[0].getScore() + alignments.back()[1].getScore() - pair_pen);
      alignments.back().setPotentialScore(
          alignments.back()[0].getPotentialScore() + alignments.back()[1].getPotentialScore() - pair_pen);
      alignments.back().setSeedChains(&chainBuilders[0].at(i0), &chainBuilders[1].at(i1));
    }
  }

  const auto best = pairBuilder.pickBest(readPair, alignments, unpaired);
  // unique alignment pair. Both ends get mapq 60
  ASSERT_NE(alignments.end(), best);
  ASSERT_EQ(460, best->at(0).getMapq());
  ASSERT_EQ(137, best->at(1).getMapq());
}

TEST(PairBuilder, pickBest2)
{
  const short                             match    = 1;
  const short                             mismatch = -1;
  const dragenos::align::SimilarityScores similarityScores(match, mismatch);
  const ReferenceDirDummy                 referenceDir;
  const align::PairBuilder                pairBuilder(similarityScores, 19, 80, 25, 0, 0, 0, false, 50, 0);

  char reference[] =
      "TCCATCGAGATGGACGCCGTTGGCGCTCTCCGTCTTTCTCCATTGCGTCGTGGCCTTGCTATTGACTCTACTGTAGACATTTTTACTTTTTATGTCCCTCTTATG";
  text2Bases(reference, reference + sizeof(reference) - 1);

  static const unsigned SEQUENCE_FLANKS = 20;
  const std::string     name            = "blah";
  const std::string     bases(
      (const char*)reference + SEQUENCE_FLANKS, (const char*)reference + sizeof(reference) - SEQUENCE_FLANKS);
  const std::string qualities =
      "#AAAAEEEEEEEEEEEEEEEEEEEEEEE6EEEEEEEEEEEEEEEEEAEAEEEEEEEEEEEEAEEE/EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE";

  uint64_t id       = 17;
  unsigned position = 0;

  ReadPair readPair;
  readPair[0].init(
      Name(name.begin(), name.end()),
      Bases(bases.begin(), bases.end()),
      Qualities(qualities.begin() + SEQUENCE_FLANKS, qualities.begin() + SEQUENCE_FLANKS + bases.length()),
      id,
      position);

  readPair[1].init(
      Name(name.begin(), name.end()),
      Bases(bases.begin(), bases.end()),
      Qualities(qualities.begin() + SEQUENCE_FLANKS, qualities.begin() + SEQUENCE_FLANKS + bases.length()),
      id,
      position);

  AlignmentPairs alignments;

  const InsertSizeParameters insertSizeParameters = {
      332, 486, 410, 332, 486, 14.0708, InsertSizeParameters::Orientation::pe_orient_fr_c};
  std::array<map::ChainBuilder, 2> chainBuilders{map::ChainBuilder(4.0), map::ChainBuilder(4.0)};

  chainBuilders[0].addSeedPosition(SeedPosition(Seed(&readPair[0], 0, 17), 0x2808c, 0), false, false);
  chainBuilders[0].addSeedPosition(SeedPosition(Seed(&readPair[0], 49, 17), 0x280bd, 0), false, false);

  chainBuilders[1].addSeedPosition(SeedPosition(Seed(&readPair[1], 0, 17), 0x2816b, 0), true, false);
  chainBuilders[1].addSeedPosition(SeedPosition(Seed(&readPair[1], 49, 17), 0x2813a, 0), true, false);

  chainBuilders[1].addSeedPosition(SeedPosition(Seed(&readPair[1], 0, 17), 164343, 0), true, false);
  chainBuilders[1].addSeedPosition(SeedPosition(Seed(&readPair[1], 49, 17), 164294, 0), true, false);

  chainBuilders[1].addSeedPosition(SeedPosition(Seed(&readPair[1], 30, 17), 0x28265, 0), true, false);
  chainBuilders[1].addSeedPosition(SeedPosition(Seed(&readPair[1], 49, 17), 0x28252, 0), true, false);

  chainBuilders[1].addSeedPosition(SeedPosition(Seed(&readPair[1], 0, 17), 0x280df, 0), true, false);

  //  chainBuilders[1].addSeedPosition(SeedPosition(Seed(&readPair[1], 30, 17), 164470, 0), false, false);
  //  chainBuilders[1].addSeedPosition(SeedPosition(Seed(&readPair[1], 49, 17), 164489, 0), false, false);
  //
  //
  //  chainBuilders[1].addSeedPosition(SeedPosition(Seed(&readPair[1], 0, 17), 164580, 0), true, false);

  std::array<Alignments, 2> unpaired;
  Alignment                 a11, a21, a22, a23, a24;
  a11.setScore(66);
  a11.setCigarOperations(std::to_string(bases.length()) + "M");
  a21.setScore(66);
  a21.setCigarOperations(std::to_string(bases.length()) + "M");
  a22.setScore(61);
  a22.setCigarOperations(std::to_string(bases.length()) + "M");
  a22.setFlags(align::AlignmentHeader::REVERSE_COMPLEMENT);
  a23.setScore(36);
  a23.setCigarOperations(std::to_string(bases.length()) + "M");
  a24.setScore(25);
  a24.setCigarOperations(std::to_string(bases.length()) + "M");
  unpaired[0].append(a11);
  unpaired[1].append(a21);
  unpaired[1].append(a22);
  unpaired[1].append(a23);
  unpaired[1].append(a24);

  for (unsigned i0 = 0; unpaired[0].size() > i0; ++i0) {
    for (unsigned i1 = 0; unpaired[1].size() > i1; ++i1) {
      alignments.insert(alignments.end(), AlignmentPair(unpaired[0].at(i0), unpaired[1].at(i1)));
      const bool isPair =
          pairMatch(insertSizeParameters, readPair, chainBuilders[0].at(i0), chainBuilders[1].at(i1));
      int       insert_len  = 0;
      int       insert_diff = 0;
      const int pair_pen = pairBuilder.computePairPenalty(
          insertSizeParameters, readPair, &unpaired[0].at(i0), &unpaired[1].at(i1), isPair, insert_len, insert_diff);
      alignments.back().setScore(
          alignments.back()[0].getScore() + alignments.back()[1].getScore() - pair_pen);
      alignments.back().setPotentialScore(
          alignments.back()[0].getPotentialScore() + alignments.back()[1].getPotentialScore() - pair_pen);
      alignments.back().setSeedChains(&chainBuilders[0].at(i0), &chainBuilders[1].at(i1));
    }
  }

  const auto best = pairBuilder.pickBest(readPair, alignments, unpaired);
  // unique alignment pair. Both ends get mapq 60
  ASSERT_NE(alignments.end(), best);
  ASSERT_EQ(488, best->at(0).getMapq());
  ASSERT_EQ(103, best->at(1).getMapq());
}

TEST(PairBuilder, pickBest3)
{
  const short                             match    = 1;
  const short                             mismatch = -1;
  const dragenos::align::SimilarityScores similarityScores(match, mismatch);
  const ReferenceDirDummy                 referenceDir;
  const align::PairBuilder                pairBuilder(similarityScores, 19, 80, 25, 0, 0, 0, false, 50, 0);

  char reference[] =
      "TCCATCGAGATGGACGCCGTTGGCGCTCTCCGTCTTTCTCCATTGCGTCGTGGCCTTGCTATTGACTCTACTGTAGACATTTTTACTTTTTATGTCCCTCTTATG";
  text2Bases(reference, reference + sizeof(reference) - 1);

  static const unsigned SEQUENCE_FLANKS = 20;
  const std::string     name            = "blah";
  const std::string     bases(
      (const char*)reference + SEQUENCE_FLANKS, (const char*)reference + sizeof(reference) - SEQUENCE_FLANKS);
  const std::string qualities =
      "#AAAAEEEEEEEEEEEEEEEEEEEEEEE6EEEEEEEEEEEEEEEEEAEAEEEEEEEEEEEEAEEE/EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE";

  uint64_t id       = 17;
  unsigned position = 0;

  ReadPair readPair;
  readPair[0].init(
      Name(name.begin(), name.end()),
      Bases(bases.begin(), bases.end()),
      Qualities(qualities.begin() + SEQUENCE_FLANKS, qualities.begin() + SEQUENCE_FLANKS + bases.length()),
      id,
      position);

  readPair[1].init(
      Name(name.begin(), name.end()),
      Bases(bases.begin(), bases.end()),
      Qualities(qualities.begin() + SEQUENCE_FLANKS, qualities.begin() + SEQUENCE_FLANKS + bases.length()),
      id,
      position);

  AlignmentPairs alignments;

  const InsertSizeParameters insertSizeParameters = {
      32, 486, 410, 332, 486, 14.0708, InsertSizeParameters::Orientation::pe_orient_fr_c};
  std::array<map::ChainBuilder, 2> chainBuilders{map::ChainBuilder(4.0), map::ChainBuilder(4.0)};

  chainBuilders[0].addSeedPosition(SeedPosition(Seed(&readPair[0], 0, 17), 0x2808c, 0), false, false);
  chainBuilders[0].addSeedPosition(SeedPosition(Seed(&readPair[0], 49, 17), 0x280bd, 0), false, false);

  chainBuilders[1].addSeedPosition(SeedPosition(Seed(&readPair[1], 0, 17), 0x2816b, 0), true, false);
  chainBuilders[1].addSeedPosition(SeedPosition(Seed(&readPair[1], 49, 17), 0x2813a, 0), true, false);

  chainBuilders[1].addSeedPosition(SeedPosition(Seed(&readPair[1], 0, 17), 164343, 0), true, false);
  chainBuilders[1].addSeedPosition(SeedPosition(Seed(&readPair[1], 49, 17), 164294, 0), true, false);

  chainBuilders[1].addSeedPosition(SeedPosition(Seed(&readPair[1], 30, 17), 0x28265, 0), true, false);
  chainBuilders[1].addSeedPosition(SeedPosition(Seed(&readPair[1], 49, 17), 0x28252, 0), true, false);

  chainBuilders[1].addSeedPosition(SeedPosition(Seed(&readPair[1], 0, 17), 0x280df, 0), true, false);

  //  chainBuilders[1].addSeedPosition(SeedPosition(Seed(&readPair[1], 30, 17), 164470, 0), false, false);
  //  chainBuilders[1].addSeedPosition(SeedPosition(Seed(&readPair[1], 49, 17), 164489, 0), false, false);
  //
  //
  //  chainBuilders[1].addSeedPosition(SeedPosition(Seed(&readPair[1], 0, 17), 164580, 0), true, false);

  std::array<Alignments, 2> unpaired;
  Alignment                 a11, a21, a22, a23, a24;
  a11.setScore(66);
  a11.setCigarOperations(std::string(bases.length(), 'M'));
  a21.setScore(66);
  a21.setCigarOperations(std::to_string(bases.length()) + "M");
  a22.setScore(61);
  a22.setCigarOperations(std::to_string(bases.length()) + "M");
  a22.setFlags(align::AlignmentHeader::REVERSE_COMPLEMENT);
  a23.setScore(36);
  a23.setCigarOperations(std::to_string(bases.length()) + "M");
  a24.setScore(25);
  a24.setCigarOperations(std::to_string(bases.length()) + "M");
  unpaired[0].append(a11);
  unpaired[1].append(a21);
  unpaired[1].append(a22);
  unpaired[1].append(a23);
  unpaired[1].append(a24);

  for (unsigned i0 = 0; unpaired[0].size() > i0; ++i0) {
    for (unsigned i1 = 0; unpaired[1].size() > i1; ++i1) {
      alignments.insert(alignments.end(), AlignmentPair(unpaired[0].at(i0), unpaired[1].at(i1)));
      const bool isPair =
          pairMatch(insertSizeParameters, readPair, chainBuilders[0].at(i0), chainBuilders[1].at(i1));
      int       insert_len  = 0;
      int       insert_diff = 0;
      const int pair_pen = pairBuilder.computePairPenalty(
          insertSizeParameters, readPair, &unpaired[0].at(i0), &unpaired[1].at(i1), isPair, insert_len, insert_diff);
      alignments.back().setScore(
          alignments.back()[0].getScore() + alignments.back()[1].getScore() - pair_pen);
      alignments.back().setPotentialScore(
          alignments.back()[0].getPotentialScore() + alignments.back()[1].getPotentialScore() - pair_pen);
      alignments.back().setSeedChains(&chainBuilders[0].at(i0), &chainBuilders[1].at(i1));
    }
  }

  const auto best = pairBuilder.pickBest(readPair, alignments, unpaired);
  // unique alignment pair. Both ends get mapq 60
  ASSERT_NE(alignments.end(), best);
  ASSERT_EQ(488, best->at(0).getMapq());
  ASSERT_EQ(155, best->at(1).getMapq());
}

TEST(PairBuilder, pickBest4)
{
  const short                             match    = 1;
  const short                             mismatch = -1;
  const dragenos::align::SimilarityScores similarityScores(match, mismatch);
  const ReferenceDirDummy                 referenceDir;
  const align::PairBuilder                pairBuilder(similarityScores, 19, 80, 25, 0, 0, 0, false, 50, 0);

  char reference[] =
      "TCCATCGAGATGGACGCCGTTGGCGCTCTCCGTCTTTCTCCATTGCGTCGTGGCCTTGCTATTGACTCTACTGTAGACATTTTTACTTTTTATGTCCCTCTTATG";
  text2Bases(reference, reference + sizeof(reference) - 1);

  static const unsigned SEQUENCE_FLANKS = 20;
  const std::string     name            = "blah";
  const std::string     bases(
      (const char*)reference + SEQUENCE_FLANKS, (const char*)reference + sizeof(reference) - SEQUENCE_FLANKS);
  const std::string qualities =
      "#AAAAEEEEEEEEEEEEEEEEEEEEEEE6EEEEEEEEEEEEEEEEEAEAEEEEEEEEEEEEAEEE/EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE";

  uint64_t id       = 17;
  unsigned position = 0;

  ReadPair readPair;
  readPair[0].init(
      Name(name.begin(), name.end()),
      Bases(bases.begin(), bases.end()),
      Qualities(qualities.begin() + SEQUENCE_FLANKS, qualities.begin() + SEQUENCE_FLANKS + bases.length()),
      id,
      position);

  readPair[1].init(
      Name(name.begin(), name.end()),
      Bases(bases.begin(), bases.end()),
      Qualities(qualities.begin() + SEQUENCE_FLANKS, qualities.begin() + SEQUENCE_FLANKS + bases.length()),
      id,
      position);

  AlignmentPairs alignments;

  const InsertSizeParameters insertSizeParameters = {
      32, 486, 310, 332, 486, 14.0708, InsertSizeParameters::Orientation::pe_orient_fr_c};
  std::array<map::ChainBuilder, 2> chainBuilders{map::ChainBuilder(4.0), map::ChainBuilder(4.0)};

  chainBuilders[0].addSeedPosition(SeedPosition(Seed(&readPair[0], 0, 17), 0x2808c, 0), false, false);
  chainBuilders[0].addSeedPosition(SeedPosition(Seed(&readPair[0], 49, 17), 0x280bd, 0), false, false);

  chainBuilders[1].addSeedPosition(SeedPosition(Seed(&readPair[1], 0, 17), 0x2816b, 0), true, false);
  chainBuilders[1].addSeedPosition(SeedPosition(Seed(&readPair[1], 49, 17), 0x2813a, 0), true, false);

  chainBuilders[1].addSeedPosition(SeedPosition(Seed(&readPair[1], 0, 17), 164343, 0), true, false);
  chainBuilders[1].addSeedPosition(SeedPosition(Seed(&readPair[1], 49, 17), 164294, 0), true, false);

  chainBuilders[1].addSeedPosition(SeedPosition(Seed(&readPair[1], 30, 17), 0x28265, 0), true, false);
  chainBuilders[1].addSeedPosition(SeedPosition(Seed(&readPair[1], 49, 17), 0x28252, 0), true, false);

  chainBuilders[1].addSeedPosition(SeedPosition(Seed(&readPair[1], 0, 17), 0x280df, 0), true, false);

  //  chainBuilders[1].addSeedPosition(SeedPosition(Seed(&readPair[1], 30, 17), 164470, 0), false, false);
  //  chainBuilders[1].addSeedPosition(SeedPosition(Seed(&readPair[1], 49, 17), 164489, 0), false, false);
  //
  //
  //  chainBuilders[1].addSeedPosition(SeedPosition(Seed(&readPair[1], 0, 17), 164580, 0), true, false);

  std::array<Alignments, 2> unpaired;
  Alignment                 a11, a21, a22, a23, a24;
  a11.setScore(66);
  a11.setCigarOperations(std::to_string(bases.length()) + "M");
  a21.setScore(66);
  a21.setCigarOperations(std::to_string(bases.length()) + "M");
  a22.setScore(61);
  a22.setCigarOperations(std::to_string(bases.length()) + "M");
  a22.setFlags(align::AlignmentHeader::REVERSE_COMPLEMENT);
  a23.setScore(36);
  a23.setCigarOperations(std::to_string(bases.length()) + "M");
  a24.setScore(25);
  a24.setCigarOperations(std::to_string(bases.length()) + "M");
  unpaired[0].append(a11);
  unpaired[1].append(a21);
  unpaired[1].append(a22);
  unpaired[1].append(a23);
  unpaired[1].append(a24);

  for (unsigned i0 = 0; unpaired[0].size() > i0; ++i0) {
    for (unsigned i1 = 0; unpaired[1].size() > i1; ++i1) {
      alignments.insert(alignments.end(), AlignmentPair(unpaired[0].at(i0), unpaired[1].at(i1)));
      const bool isPair =
          pairMatch(insertSizeParameters, readPair, chainBuilders[0].at(i0), chainBuilders[1].at(i1));
      int       insert_len  = 0;
      int       insert_diff = 0;
      const int pair_pen = pairBuilder.computePairPenalty(
          insertSizeParameters, readPair, &unpaired[0].at(i0), &unpaired[1].at(i1), isPair, insert_len, insert_diff);
      alignments.back().setScore(
          alignments.back()[0].getScore() + alignments.back()[1].getScore() - pair_pen);
      alignments.back().setPotentialScore(
          alignments.back()[0].getPotentialScore() + alignments.back()[1].getPotentialScore() - pair_pen);
      alignments.back().setSeedChains(&chainBuilders[0].at(i0), &chainBuilders[1].at(i1));
    }
  }

  const auto best = pairBuilder.pickBest(readPair, alignments, unpaired);
  // unique alignment pair. Both ends get mapq 60
  ASSERT_NE(alignments.end(), best);
  ASSERT_EQ(488, best->at(0).getMapq());
  ASSERT_EQ(51, best->at(1).getMapq());
}
