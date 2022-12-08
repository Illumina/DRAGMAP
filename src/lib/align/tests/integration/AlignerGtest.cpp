
#include <boost/filesystem.hpp>
#include <sstream>

#include "gtest/gtest.h"

#include "align/Aligner.hpp"
#include "align/SinglePicker.hpp"
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

TEST(Aligner, getAlignments)
{
  using namespace dragenos;
  typedef align::Aligner::Read Read;
  //  typedef Read::Name Name;
  //  typedef Read::Bases Bases;
  //  typedef Read::Qualities Qualities;
  //  typedef align::Aligner::ReadPair ReadPair;
  typedef align::Aligner::Alignments Alignments;
  //  typedef align::Aligner::AlignmentPair AlignmentPair;
  const short                             match    = 1;
  const short                             mismatch = -1;
  const dragenos::align::SimilarityScores similarityScores(match, mismatch);
  //  const ReferenceDirDummy referenceDir;
  //  const align::Aligner aligner(referenceDir, similarityScores, gapInit, gapExtend);
  //  std::istringstream fastq("");
  Read              singleRead;
  const std::string name = "blah";
  //  const unsigned char reference[] =
  //    "ACGTACGTACGTACGTACGT"
  //    "ACGTACGTACGTACGTACGT"
  //    "ACGTACGTACGTACGTACGT";
  //  const std::string bases = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";

  //  const unsigned char reference[] = "AAAACCCCCGGGGGTTTTACGT""AACGTAAAAA";
  //  const std::string bases =         "AAAACCCCCGGGGGTTTTACGT"  "CGTAAAAA";
  //  const unsigned char reference[] = "GTTCCG""ACGTA";
  //  const std::string bases =         "GTTCCG"  "GTA";
  const std::string bases =
      "GTTCCG"
      "ACGTA";

  const std::string qualities;
  // init requires rvalue references for all parameters
  //  uint64_t id = 17;
  //  unsigned position = 0;
  //  singleRead.init(
  //    Name(name.begin(), name.end()), Bases(bases.begin(), bases.end()),
  //    Qualities(qualities.begin(), qualities.end()), id, position);
  Alignments   alignments;
  align::Query q(bases.begin(), bases.end());
  alignments.addAlignment();

  //  aligner.getAlignments(singleRead, alignment);
  //  ReadPair readPair;
  //  readPair[0].init(Name(), Bases(), Qualities(), id, 0);
  //  readPair[1].init(Name(), Bases(), Qualities(), id, 1);
  //  AlignmentPair alignmentPair;
  //  aligner.getAlignments(readPair, alignmentPair);
}

using namespace dragenos;
typedef align::Aligner::Read       Read;
typedef Read::Name                 Name;
typedef Read::Bases                Bases;
typedef Read::Qualities            Qualities;
typedef align::Aligner::Alignments Alignments;

TEST(Aligner, updateMapqSingleEnded)
{
  const short                             match    = 1;
  const short                             mismatch = -1;
  const dragenos::align::SimilarityScores similarityScores(match, mismatch);
  const ReferenceDirDummy                 referenceDir;
  const align::SinglePicker               picker(similarityScores, 19, 8, 0, 0, 0, false, 50, 0);

  const unsigned char reference[] =
      "TCCATCGAGATGGACGCCGTTGGCGCTCTCCGTCTTTCTCCATTGCGTCGTGGCCTTGCTATTGACTCTACTGTAGACATTTTTACTTTTTATGTCCCTC";

  static const unsigned SEQUENCE_FLANKS = 20;
  const std::string     name            = "blah";
  const std::string     bases(
      (const char*)reference + SEQUENCE_FLANKS, (const char*)reference + sizeof(reference) - SEQUENCE_FLANKS);
  //  const std::string qualities = "#AAAAEEEEEEEEEEEEEEEEEEEEEEE6EEEEEEEEEEEEEEEEEAEAEEEEEEEEEEEEAEEE/EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE";
  const std::string qualities =
      "#AAAAEEEEEEEEEEEEEEEEEEEEEEE6EEEEEEEEEEEEEEEEEAEAEEEEEEEEEEEEAEEE/EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE";

  uint64_t id       = 17;
  unsigned position = 0;

  Read singleRead;
  singleRead.init(
      Name(name.begin(), name.end()),
      Bases(bases.begin(), bases.end()),
      Qualities(qualities.begin() + SEQUENCE_FLANKS, qualities.begin() + SEQUENCE_FLANKS + bases.length()),
      id,
      position);

  Alignments alignments;

  alignments.addAlignment();
  alignments.back().setScore(100);
  Alignments::iterator best = picker.pickBest(singleRead.getLength(), alignments);
  ASSERT_NE(alignments.end(), best);
  ASSERT_EQ(875, alignments.back().getMapq());

  alignments.addAlignment();
  alignments.back().setScore(100);
  alignments.back().setFlags(align::AlignmentHeader::REVERSE_COMPLEMENT);
  alignments.back().setCigarOperations(
      std::string('M', bases.length() / 2) + std::string('S', bases.length() / 2));
  best = picker.pickBest(singleRead.getLength(), alignments);
  ASSERT_NE(alignments.end(), best);
  ASSERT_EQ(0, alignments.front().getMapq());
}

TEST(Aligner, updateMapqSingleEnded1XRepeat)
{
  const short                             match    = 1;
  const short                             mismatch = -1;
  const dragenos::align::SimilarityScores similarityScores(match, mismatch);
  const ReferenceDirDummy                 referenceDir;
  const align::SinglePicker               picker(similarityScores, 19, 8, 0, 0, 0, false, 50, 0);

  const unsigned char reference[] =
      "TCCATCGAGATGGACGCCGTTGGCGCTCTCCGTCTTTCTCCATTGCGTCGTGGCCTTGCTATTGACTCTACTGTAGACATTTTTACTTTTTATGTCCCTC";

  static const unsigned SEQUENCE_FLANKS = 20;
  const std::string     name            = "blah";
  const std::string     bases(
      (const char*)reference + SEQUENCE_FLANKS, (const char*)reference + sizeof(reference) - SEQUENCE_FLANKS);
  //  const std::string qualities = "#AAAAEEEEEEEEEEEEEEEEEEEEEEE6EEEEEEEEEEEEEEEEEAEAEEEEEEEEEEEEAEEE/EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE";
  const std::string qualities =
      "#AAAAEEEEEEEEEEEEEEEEEEEEEEE6EEEEEEEEEEEEEEEEEAEAEEEEEEEEEEEEAEEE/EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE";

  uint64_t id       = 17;
  unsigned position = 0;

  Read singleRead;
  singleRead.init(
      Name(name.begin(), name.end()),
      Bases(bases.begin(), bases.end()),
      Qualities(qualities.begin() + SEQUENCE_FLANKS, qualities.begin() + SEQUENCE_FLANKS + bases.length()),
      id,
      position);

  Alignments alignments;

  align::Alignment primary(0, 101);

  alignments.append(primary);
  Alignments::iterator best = picker.pickBest(singleRead.getLength(), alignments);
  ASSERT_NE(alignments.end(), best);
  ASSERT_EQ(886, alignments.back().getMapq());

  align::Alignment secondBest(0, 99);
  secondBest.setCigarOperations(std::to_string(bases.length()) + "M");
  alignments.append(secondBest);
  alignments.back().setFlags(align::AlignmentHeader::REVERSE_COMPLEMENT);
  best = picker.pickBest(singleRead.getLength(), alignments);
  ASSERT_NE(alignments.end(), best);
  ASSERT_EQ(21, alignments.front().getMapq());
}

TEST(Aligner, isPerfectAlignment)
{
  dragenos::align::Alignment alignment;
  using dragenos::align::Aligner;
  constexpr char A   = 1;
  constexpr char C   = 2;
  constexpr char G   = 4;
  constexpr char T   = 8;
  constexpr char N   = 0xF;
  char           q[] = {N, A, C, G, T, A, A, C, C, G, G, T, T};
  char           d[] = {N, A, C, G, T, A, A, C, C, G, G, T, T};
  // both empty
  ASSERT_FALSE(Aligner::isPerfectAlignment(q, q, d, d, alignment));
  // different lengths
  ASSERT_FALSE(Aligner::isPerfectAlignment(q + 1, q + 2, d + 1, d + 3, alignment));
  ASSERT_FALSE(Aligner::isPerfectAlignment(q + 1, q + 3, d + 1, d + 2, alignment));
  // ok
  ASSERT_TRUE(Aligner::isPerfectAlignment(q + 1, q + 3, d + 1, d + 3, alignment));
  ASSERT_EQ(alignment.getScore(), 2);
}
