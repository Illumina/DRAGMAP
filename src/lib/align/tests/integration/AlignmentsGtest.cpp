
#include <debug/vector>

#include "gtest/gtest.h"

#include "align/Alignments.hpp"

TEST(Alignments, OneBaseDeletion)
{
  using namespace dragenos;

  const short                             match    = 1;
  const short                             mismatch = -1;
  const dragenos::align::SimilarityScores similarityScores(match, mismatch);
  const int                               gapInit   = 2;
  const int                               gapExtend = 1;

  const char reference[] =
      "GTTCCG"
      "ACGTAAA";
  const std::string q =
      "GTTCCG"
      "CGTA";
  align::SmithWatermanT<char, short, 8, 16, 4> sw(similarityScores, gapInit, gapExtend);
  std::string                                  result;
  sw.align(q.data(), q.data() + q.size(), reference, reference + sizeof(reference) - 1, 10, 2, false, result);
  ASSERT_EQ(std::string("MMMMMMDMMMM"), result);
}

TEST(Alignments, TwoBaseDeletion)
{
  using namespace dragenos;

  const short                             match    = 1;
  const short                             mismatch = -1;
  const dragenos::align::SimilarityScores similarityScores(match, mismatch);
  const int                               gapInit   = 2;
  const int                               gapExtend = 1;

  const char reference[] =
      "GTTCCG"
      "ACGTAAA";
  const std::string expected =
      "MMMMMM"
      "DDMMMM";
  const std::string q =
      "GTTCCG"
      "GTAA";

  align::SmithWatermanT<char, short, 8, 16, 4> sw(similarityScores, gapInit, gapExtend);
  std::string                                  result;
  sw.align(q.data(), q.data() + q.size(), reference, reference + sizeof(reference) - 1, 10, 2, false, result);
  ASSERT_EQ(expected, result);
}

TEST(Alignments, OneBaseInsertion)
{
  using namespace dragenos;

  const short                             match    = 1;
  const short                             mismatch = -1;
  const dragenos::align::SimilarityScores similarityScores(match, mismatch);
  const int                               gapInit   = 2;
  const int                               gapExtend = 1;

  const char reference[] =
      "GTTCCG"
      "CGTTTT";
  const std::string expected =
      "MMMMMM"
      "IMMM";
  const std::string q =
      "GTTCCG"
      "ACGT";

  align::SmithWatermanT<char, short, 8, 16, 4> sw(similarityScores, gapInit, gapExtend);
  std::string                                  result;
  sw.align(q.data(), q.data() + q.size(), reference, reference + sizeof(reference) - 1, 10, 2, false, result);
  ASSERT_EQ(expected, result);
}

TEST(Alignments, TwoBaseInsertion)
{
  using namespace dragenos;

  const short                             match    = 1;
  const short                             mismatch = -1;
  const dragenos::align::SimilarityScores similarityScores(match, mismatch);
  const int                               gapInit   = 2;
  const int                               gapExtend = 1;

  const char reference[] =
      "GTTCCG"
      "GTAAAT";
  const std::string expected = "MMMMMMIIMMMM";
  const std::string q        = "GTTCCGACGTAA";

  align::SmithWatermanT<char, short, 8, 16, 4> sw(similarityScores, gapInit, gapExtend);
  std::string                                  result;
  sw.align(q.data(), q.data() + q.size(), reference, reference + sizeof(reference) - 1, 10, 2, false, result);
  ASSERT_EQ(expected, result);
}

TEST(Alignments, SoftClipAtEnd)
{
  using namespace dragenos;

  const short                             match    = 1;
  const short                             mismatch = -1;
  const dragenos::align::SimilarityScores similarityScores(match, mismatch);
  const int                               gapInit   = 2;
  const int                               gapExtend = 1;

  const char reference[] =
      "GTTCCG"
      "GTAAATTTTTTTTTTT";
  const std::string expected = "MMMMMMIIMMMMSSSSSS";
  const std::string q        = "GTTCCGACGTAAGGGGGG";

  align::SmithWatermanT<char, short, 8, 16, 4> sw(similarityScores, gapInit, gapExtend);
  std::string                                  result;
  sw.align(q.data(), q.data() + q.size(), reference, reference + sizeof(reference) - 1, 10, 3, false, result);
  ASSERT_EQ(expected, result);
}

TEST(Alignments, MismatchSoftClipAtStart)
{
  using namespace dragenos;

  const short                             match    = 1;
  const short                             mismatch = -1;
  const dragenos::align::SimilarityScores similarityScores(match, mismatch);
  const int                               gapInit   = 2;
  const int                               gapExtend = 1;

  const char reference[] =
      "AAAAA"
      "GTTCCG"
      "GTAAATTTTTTTTTTT";
  const std::string expected = "NNNNNSSSSSMMMMMMIIMMMM";
  const std::string q        = "GGGGGGTTCCGACGTAA";

  align::SmithWatermanT<char, short, 8, 16, 4> sw(similarityScores, gapInit, gapExtend);
  std::string                                  result;
  sw.align(q.data(), q.data() + q.size(), reference, reference + sizeof(reference) - 1, 10, 1, false, result);
  ASSERT_EQ(expected, result);
}

TEST(Alignments, ShortInsertionSoftClipAtStart)
{
  using namespace dragenos;

  const short                             match    = 1;
  const short                             mismatch = -1;
  const dragenos::align::SimilarityScores similarityScores(match, mismatch);
  const int                               gapInit   = 2;
  const int                               gapExtend = 1;

  const char reference[] =     "G"    "TTCCG""GTAAATTTTTTTTTTT";
  const std::string expected = "NSSSSSSMMMMMIIMMMM";
  const std::string q        =  "GAAAGGTTCCGACGTAA";

  std::string                                  result;
  align::SmithWatermanT<char, short, 8, 16, 4> sw(similarityScores, gapInit, gapExtend);
  sw.align(q.data(), q.data() + q.size(), reference, reference + sizeof(reference) - 1, 10, 1, false, result);
  ASSERT_EQ(expected, result);
}

TEST(Alignments, InsertionSoftClipAtStart)
{
  using namespace dragenos;

  const short                             match    = 1;
  const short                             mismatch = -1;
  const dragenos::align::SimilarityScores similarityScores(match, mismatch);
  const int                               gapInit   = 2;
  const int                               gapExtend = 1;

  const char reference[] =        "G"    "TTCCG""GTAAATTTTTTTTTTT";
  const std::string expected =   "NSSSSSSSMMMMMIIMMMM";
  const std::string q        =    "GAAAGGGTTCCGACGTAA";

  align::SmithWatermanT<char, short, 8, 16, 4> sw(similarityScores, gapInit, gapExtend);
  std::string                                  result;
  sw.align(q.data(), q.data() + q.size(), reference, reference + sizeof(reference) - 1, 0, 0, false, result);
  ASSERT_EQ(expected, result);
}

TEST(Alignments, SkipAtStart)
{
  using namespace dragenos;

  const short                             match    = 1;
  const short                             mismatch = -1;
  const dragenos::align::SimilarityScores similarityScores(match, mismatch);
  const int                               gapInit   = 2;
  const int                               gapExtend = 1;

  const char reference[] =
      "GAAAGGGTTCCG"
      "GTAAATTTTTTTTTTT";
  const std::string expected = "NNNNMMMMMMMMIIMMMM";
  const std::string q        = "GGGTTCCGACGTAA";

  align::SmithWatermanT<char, short, 8, 16, 4> sw(similarityScores, gapInit, gapExtend);
  std::string                                  result;
  sw.align(q.data(), q.data() + q.size(), reference, reference + sizeof(reference) - 1, 10, 8, false, result);
  ASSERT_EQ(expected, result);
}

TEST(Alignments, SkipAtStart2)
{
  using namespace dragenos;

  const short                             match    = 1;
  const short                             mismatch = -1;
  const dragenos::align::SimilarityScores similarityScores(match, mismatch);
  const int                               gapInit   = 2;
  const int                               gapExtend = 1;

  const char reference[] =     "AAAGGGTTCCG""GTAAATTTTTTTTTTT";
  const std::string expected = "NNNMMMMMMMMIIMMMM";
  const std::string q        =    "GGGTTCCGACGTAA";

  align::SmithWatermanT<char, short, 8, 16, 4> sw(similarityScores, gapInit, gapExtend);
  std::string                                  result;
  sw.align(q.data(), q.data() + q.size(), reference, reference + sizeof(reference) - 1, 10, 8, false, result);
  ASSERT_EQ(expected, result);
}

TEST(Alignments, MMMMMMMMMMMMM)
{
  using namespace dragenos;

  const short                             match    = 1;
  const short                             mismatch = -1;
  const dragenos::align::SimilarityScores similarityScores(match, mismatch);
  const int                               gapInit   = 2;
  const int                               gapExtend = 1;

  const char        reference[] = "ATTCGACCTATCC";
  const std::string expected    = "MMMMMMMMMMMMM";
  const std::string q           = "ATTCGACCTATCC";

  align::SmithWatermanT<char, short, 8, 16, 4> sw(similarityScores, gapInit, gapExtend);
  std::string                                  result;
  sw.align(q.data(), q.data() + q.size(), reference, reference + sizeof(reference) - 1, 0, 0, false, result);
  ASSERT_EQ(expected, result);
}

TEST(Alignments, NNMMMMMMMMMMMMM)
{
  using namespace dragenos;

  const short                             match    = 1;
  const short                             mismatch = -1;
  const dragenos::align::SimilarityScores similarityScores(match, mismatch);
  const int                               gapInit   = 2;
  const int                               gapExtend = 1;

  const char        reference[] = "CCATTCGACCTATCC";
  const std::string expected    = "NNMMMMMMMMMMMMM";
  const std::string q           = "ATTCGACCTATCC";

  align::SmithWatermanT<char, short, 8, 16, 4> sw(similarityScores, gapInit, gapExtend);
  std::string                                  result;
  sw.align(q.data(), q.data() + q.size(), reference, reference + sizeof(reference) - 1, 2, 2, false, result);
  ASSERT_EQ(expected, result);
}

TEST(Alignments, NNNNNNNNMMMMMMMMMM)
{
  using namespace dragenos;

  const short                             match    = 1;
  const short                             mismatch = -1;
  const dragenos::align::SimilarityScores similarityScores(match, mismatch);
  const int                               gapInit   = 2;
  const int                               gapExtend = 1;

  const char        reference[] = "CCCCCCCCATTCGACCTA";
  const std::string expected    = "NNNNNNNNMMMMMMMMMM";
  const std::string q           = "ATTCGACCTA";

  align::SmithWatermanT<char, short, 8, 16, 4> sw(similarityScores, gapInit, gapExtend);
  std::string                                  result;
  sw.align(q.data(), q.data() + q.size(), reference, reference + sizeof(reference) - 1, 4, 4, false, result);
  ASSERT_EQ(expected, result);
}

TEST(Alignments, NNNNNNNNMMMMMMMMMMMMMMMMMMMMMMM)
{
  using namespace dragenos;

  const short                             match    = 1;
  const short                             mismatch = -1;
  const dragenos::align::SimilarityScores similarityScores(match, mismatch);
  const int                               gapInit   = 2;
  const int                               gapExtend = 1;

  const char        reference[] = "CCGCATGCACCGTTGTCTGTGCTGTGACTTC";
  const std::string expected    = "NNNNNNNNMMMMMMMMMMMMMMMMMMMMMMM";
  const std::string q           = "ACCGTTGTCTGTGCTGTGACTTC";

  align::SmithWatermanT<char, short, 8, 16, 10> sw(similarityScores, gapInit, gapExtend);
  std::string                                   result;
  sw.align(q.data(), q.data() + q.size(), reference, reference + sizeof(reference) - 1, 10, 2, false, result);
  ASSERT_EQ(expected, result);
}

TEST(Alignments, OffsetOutOfBoundsWhenBacktracking)
{
  using namespace dragenos;

  constexpr int                           MATCH    = 1;
  constexpr int                           MISMATCH = -4;
  const dragenos::align::SimilarityScores similarityScores(MATCH, MISMATCH);
  constexpr int                           gapExtend = 1;
  constexpr int                           gapInit   = 7;

  const char reference[] = "TCACACAGACTGGCCCCCCAAACACCCACACAGACCATCCCCCCAACACCCACATAGACCATGCCCCCAAACACCCACACAGACTGTCCCCCCAGACACCCGCACAGACCA";
  //  const std::string expected = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNMMMMMXMMMMMMMXMMMMMMXMMMMMMMMMMMMMMXXMMMMMMMMMMMMMMMMMMMMMMMMMSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS";
  const std::string expected =
      "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS";
  const std::string q =
      "CCACACAGACCATCCCCCCAGACACCCACACAGACCATCCCCCCAGACACCCGCACAGACCATCCCCCCAGACAACCTACACAGAATAAACTGTCCCC";

  align::SmithWatermanT<char, short, 48, 16, 10> sw(similarityScores, gapInit, gapExtend);
  std::string                                    result;
  sw.align(q.data(), q.data() + q.size(), reference, reference + sizeof(reference) - 1, 1, 48, false, result);
  ASSERT_EQ(expected, result);
}

TEST(Alignments, Short)
{
  using namespace dragenos;

  constexpr int                           MATCH    = 1;
  constexpr int                           MISMATCH = -4;
  const dragenos::align::SimilarityScores similarityScores(MATCH, MISMATCH);
  constexpr int                           gapExtend = 1;
  constexpr int                           gapInit   = 7;

  const int  width       = 48;
  const char reference[] = "GCAAAAGTCCCTTCACTTCCACCACCTCCCAGAGGCATGCCTTCTTTCACCAACTCGC";
  //  const std::string expected = "MMMMMMMMMMMMMMMMMMMMMMMMXMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM";
  const std::string expected = "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM";
  const std::string q        = "GCAAAAGTCCCTTCACTTCCACCATCTCCCAGAGGCATGCCTTCTTTCACCAACTCGC";

  align::SmithWatermanT<char, short, width, 16, 10> sw(similarityScores, gapInit, gapExtend);
  std::string                                       result;
  sw.align(q.data(), q.data() + q.size(), reference, reference + sizeof(reference) - 1, 11, 2, false, result);
  ASSERT_EQ(expected, result);
}

TEST(Alignments, VeryShortQuery)
{
  using namespace dragenos;

  constexpr int                           MATCH    = 1;
  constexpr int                           MISMATCH = -4;
  const dragenos::align::SimilarityScores similarityScores(MATCH, MISMATCH);
  constexpr int                           gapExtend = 1;
  constexpr int                           gapInit   = 7;

  const int         width       = 48;
  const char        reference[] = "AAAAAAAGCAAAAGTCCCTTCACTTCCACCACCTCCCAGAGGCATGCCTTCTTTCACCAACTCGC";
  const std::string expected    = "NNNNNNNMMMMMMMMMMMMMMM";
  const std::string q           = "GCAAAAGTCCCTTCA";

  align::SmithWatermanT<char, short, width, 16, 10> sw(similarityScores, gapInit, gapExtend);
  std::string                                       result;
  sw.align(q.data(), q.data() + q.size(), reference, reference + sizeof(reference) - 1, 11, 4, false, result);
  ASSERT_EQ(expected, result);
}

TEST(Alignments, TinyQuery)
{
  using namespace dragenos;

  constexpr int                           MATCH    = 1;
  constexpr int                           MISMATCH = -4;
  const dragenos::align::SimilarityScores similarityScores(MATCH, MISMATCH);
  constexpr int                           gapExtend = 1;
  constexpr int                           gapInit   = 7;

  const int         width       = 48;
  const char        reference[] = "AAAAAAAGCAAAAGTCCCTTCACTTCCACCACCTCCCAGAGGCATGCCTTCTTTCACCAACTCGC";
  const std::string expected    = "NNNNNNNM";
  const std::string q           = "G";

  align::SmithWatermanT<char, short, width, 16, 10> sw(similarityScores, gapInit, gapExtend);
  std::string                                       result;
  sw.align(q.data(), q.data() + q.size(), reference, reference + sizeof(reference) - 1, 11, 2, false, result);
  ASSERT_EQ(expected, result);
}

TEST(Alignments, TinyQuery2)
{
  using namespace dragenos;

  constexpr int                           MATCH    = 1;
  constexpr int                           MISMATCH = -4;
  const dragenos::align::SimilarityScores similarityScores(MATCH, MISMATCH);
  constexpr int                           gapExtend = 1;
  constexpr int                           gapInit   = 7;

  const int         width       = 48;
  const char        reference[] = "CAATACCTCACTCAGATTCCATTATGCCAAATAATTAGCAAGGTGACAAAAGCTCTGCATGAGCTG";
  const std::string expected    = "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM";
  const std::string q           = "CAATACCTCACTCAGATTCCATTATGCCAAATAATTAGCAAGGTGACAAAAGCTCTGCATGAGCTG";

  align::SmithWatermanT<char, short, width, 16, 10> sw(similarityScores, gapInit, gapExtend);
  std::string                                       result;
  sw.align(q.data(), q.data() + q.size(), reference, reference + sizeof(reference) - 1, 11, 2, false, result);
  ASSERT_EQ(expected, result);
}

TEST(Alignments, wrongNSatStart)
{
  using namespace dragenos;

  constexpr int                           MATCH    = 1;
  constexpr int                           MISMATCH = -4;
  const dragenos::align::SimilarityScores similarityScores(MATCH, MISMATCH);
  constexpr int                           gapExtend = 1;
  constexpr int                           gapInit   = 7;

  const int  width = 48;
  const char reference[] =
      "ACGGTATGTGATG"
      "TTTTTTTTTTTTTTCCTCACAAGATACTTTTC";
  const std::string expected =
      "NMMMMMMMMMMMM"
      "IMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM";
  const std::string q =
      "CGGTATGTGATG"
      "TTTTTTTTTTTTTTTCCTCACAAGATACTTTTC";

  align::SmithWatermanT<char, short, width, 16, 10> sw(similarityScores, gapInit, gapExtend);
  std::string                                       result;
  sw.align(q.data(), q.data() + q.size(), reference, reference + sizeof(reference) - 1, 11, 2, false, result);
  ASSERT_EQ(expected, result);
}

// when match is a reverse complement of the reference there is no point
// in reversing the query as s-w needs a query reversed in relation to the reference anyway.
// just inform s-w that the query is already reversed
TEST(Alignments, reverseQuery)
{
  using namespace dragenos;

  constexpr int                           MATCH    = 1;
  constexpr int                           MISMATCH = -4;
  const dragenos::align::SimilarityScores similarityScores(MATCH, MISMATCH);
  constexpr int                           gapExtend = 1;
  constexpr int                           gapInit   = 7;

  const int  width = 48;
  const char reference[] =
      "ACGGTATGTGATG"
      "TTTTTTTTTTTTTTCCTCACAAGATACTTTTC";
  const std::string expected =
      "NMMMMMMMMMMMMI"
      "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM";
  // original query              CGGTATGTGATGT""TTTTTTTTTTTTTTCCTCACAAGATACTTTTC
  const std::string q = "CTTTTCATAGAACACTCCTTTTTTTTTTTTTTTGTAGTGTATGGC";

  align::SmithWatermanT<char, short, width, 16, 10> sw(similarityScores, gapInit, gapExtend);
  std::string                                       result;
  sw.align(q.data(), q.data() + q.size(), reference, reference + sizeof(reference) - 1, 11, 2, true, result);
  ASSERT_EQ(expected, result);
}

TEST(Alignments, indelsTogether)
{
  using namespace dragenos;

  constexpr int                           MATCH    = 1;
  constexpr int                           MISMATCH = -4;
  const dragenos::align::SimilarityScores similarityScores(MATCH, MISMATCH);
  constexpr int                           gapExtend = 1;
  constexpr int                           gapInit   = 7;

  const char reference[] = "GAGTGAAAATTTGATCTACTCACAG";
  // prioritise gap when backtracking to avoid gaps interspersed with matches
  const std::string expected = "MMMMMMMMMMMDDDMMMMMMMMMMM";
  const std::string q =
      "GAGTGAAAATT"
      "TCTACTCACAG";

  align::SmithWatermanT<char, short, 8, 16, 4> sw(similarityScores, gapInit, gapExtend);
  std::string                                  result;
  sw.align(q.data(), q.data() + q.size(), reference, reference + sizeof(reference) - 1, 20, 0, false, result);
  ASSERT_EQ(expected, result);
}

TEST(Alignments, indelsTogetherSmall)
{
  using namespace dragenos;

  constexpr int                           MATCH    = 1;
  constexpr int                           MISMATCH = -4;
  const dragenos::align::SimilarityScores similarityScores(MATCH, MISMATCH);
  constexpr int                           gapExtend = 1;
  constexpr int                           gapInit   = 7;

  const char        reference[] = "ACGTACGTGAAAATTTGATCTACTT";
  const std::string expected    = "MMMMMMMMMMMMMMDDDDMMMMMMM";
  const std::string q =           "ACGTACGTGAAAAT"  "TCTACTT";

  align::SmithWatermanT<char, short, 8, 16, 10> sw(similarityScores, gapInit, gapExtend, 5);
  std::string                                   result;
  sw.align(q.data(), q.data() + q.size(), reference, reference + sizeof(reference) - 1, 11, 2, false, result);
  ASSERT_EQ(expected, result);
}

TEST(Alignments, unclipScore)
{
  using namespace dragenos;

  constexpr int                           MATCH    = 1;
  constexpr int                           MISMATCH = -4;
  const dragenos::align::SimilarityScores similarityScores(MATCH, MISMATCH);
  constexpr int                           gapExtend   = 1;
  constexpr int                           gapInit     = 7;
  constexpr int                           unclipScore = 5;

  const int  width       = 8;
  const char reference[] =     "TTTAAACTCACGAG";
  //  const std::string expected = "MMMXMMMMMMM";
  const std::string expected = "MMMMMMMMMMMMMM";
  const std::string q        = "TTTAAAATCACGAG";

  align::SmithWatermanT<char, short, width, 16, 10> sw(similarityScores, gapInit, gapExtend, unclipScore);
  std::string                                       result;
  sw.align(q.data(), q.data() + q.size(), reference, reference + sizeof(reference) - 1, 5, 2, false, result);
  ASSERT_EQ(expected, result);
}

TEST(Alignments, unclipScore0)
{
  using namespace dragenos;

  constexpr int                           MATCH    = 1;
  constexpr int                           MISMATCH = -4;
  const dragenos::align::SimilarityScores similarityScores(MATCH, MISMATCH);
  constexpr int                           gapExtend   = 1;
  constexpr int                           gapInit     = 7;
  constexpr int                           unclipScore = 0;

  const int  width = 8;
  const char reference[] =
      "TCTTTC"
      "TTTCTTT";
  const std::string expected = "NNNNNNSSSSSSSSSMMMS";
  const std::string q        = "GGGCAAAGGTTTG";

  align::SmithWatermanT<char, short, width, 16, 10> sw(similarityScores, gapInit, gapExtend, unclipScore);
  std::string                                       result;
  sw.align(q.data(), q.data() + q.size(), reference, reference + sizeof(reference) - 1, 11, 2, false, result);
  ASSERT_EQ(expected, result);
}

// TEST(Alignments, unclipScore5)
//{
//  using namespace dragenos;
//
//  constexpr int MATCH = 1;
//  constexpr int MISMATCH = -4;
//  const dragenos::align::SimilarityScores similarityScores(MATCH, MISMATCH);
//  constexpr int gapExtend = 1;
//  constexpr int gapInit = 7;
//  constexpr int unclipScore = 5;
//
//  const int width = 8;
//  const char reference[] =     "TCTTTCTTTCTTT";
//  const std::string expected = "MMMMMMMMMMMMM";
//  const std::string q =        "GGGCAAAGGTTTG";
//
//  align::SmithWatermanT<char, short, width, 16, 10> sw(similarityScores, gapInit, gapExtend, unclipScore);
//  const auto result = sw.align(q.data(), q.data() + q.size(), reference, reference + sizeof(reference) - 1, 11, 2);
//  ASSERT_EQ(expected, result);
//}

TEST(Alignments, unclipScore0Start)
{
  using namespace dragenos;

  const short                             match    = 1;
  const short                             mismatch = -1;
  const dragenos::align::SimilarityScores similarityScores(match, mismatch);
  const int                               gapInit     = 2;
  const int                               gapExtend   = 1;
  const int                               unclipScore = 0;

  const char reference[] =
      "GTTCCG"
      "ACGTAAA";
  //  const std::string expected = "NS""MMMMM""DMMMM";
  const std::string expected =
      "MMMMMM"
      "DMMMM";
  const std::string q =
      "ATTCCG"
      "CGTA";
  align::SmithWatermanT<char, short, 8, 16, 4> sw(similarityScores, gapInit, gapExtend, unclipScore);
  std::string                                  result;
  sw.align(q.data(), q.data() + q.size(), reference, reference + sizeof(reference) - 1, 2, 2, false, result);
  ASSERT_EQ(expected, result);
}

TEST(Alignments, unclipScore5Start)
{
  using namespace dragenos;

  const short                             match    = 1;
  const short                             mismatch = -1;
  const dragenos::align::SimilarityScores similarityScores(match, mismatch);
  const int                               gapInit     = 2;
  const int                               gapExtend   = 1;
  const int                               unclipScore = 5;

  const char reference[] =     "TGTTCCGA""CGTAAA";
  const std::string expected = "NMMMMMMD""MMMM";
  const std::string q =         "ATTCCG" "CGTA";
  align::SmithWatermanT<char, short, 8, 16, 4> sw(similarityScores, gapInit, gapExtend, unclipScore);
  std::string                                  result;
  sw.align(q.data(), q.data() + q.size(), reference, reference + sizeof(reference) - 1, 2, 2, false, result);
  ASSERT_EQ(expected, result);
}

TEST(Alignments, unclipScore0End)
{
  using namespace dragenos;

  const short                             match    = 1;
  const short                             mismatch = -1;
  const dragenos::align::SimilarityScores similarityScores(match, mismatch);
  const int                               gapInit     = 2;
  const int                               gapExtend   = 1;
  const int                               unclipScore = 0;

  const char reference[] =
      "GTTCCG"
      "ACGTAAA";
  const std::string expected =
      "MMMMMM"
      "DMMMS";
  const std::string q =
      "GTTCCG"
      "CGTT";
  align::SmithWatermanT<char, short, 8, 16, 4> sw(similarityScores, gapInit, gapExtend, unclipScore);
  std::string                                  result;
  sw.align(q.data(), q.data() + q.size(), reference, reference + sizeof(reference) - 1, 2, 2, false, result);
  ASSERT_EQ(expected, result);
}

TEST(Alignments, unclipScore5End)
{
  using namespace dragenos;

  const short                             match    = 1;
  const short                             mismatch = -1;
  const dragenos::align::SimilarityScores similarityScores(match, mismatch);
  const int                               gapInit     = 2;
  const int                               gapExtend   = 1;
  const int                               unclipScore = 5;

  const char reference[] =     "GGTTCCG""ACGTAAA";
  const std::string expected = "NMMMMMM""DMMMM";
  const std::string q =         "GTTCCG" "CGTT";
  align::SmithWatermanT<char, short, 8, 16, 4> sw(similarityScores, gapInit, gapExtend, unclipScore);
  std::string                                  result;
  sw.align(q.data(), q.data() + q.size(), reference, reference + sizeof(reference) - 1, 2, 2, false, result);
  ASSERT_EQ(expected, result);
}

TEST(Alignments, SSMMMMMMMMMMDMMMMMMMMMMMMMMMMMMMMMM)
{
  using namespace dragenos;

  constexpr int                           MATCH    = 1;
  constexpr int                           MISMATCH = -4;
  const dragenos::align::SimilarityScores similarityScores(MATCH, MISMATCH);
  constexpr int                           gapExtend = 1;
  constexpr int                           gapInit   = 7;

  const int  width = 8;
  const char reference[] =     "T""TGGGCCGATTA""AAAAAAAAAAAAAAACATAAAA";
  const std::string expected = "NSSMMMMMMMMMMD""MMMMMMMMMMMMMMMMMMMMMM";
  const std::string q =         "TATGGGCCGATT" "AAAAAAAAAAAAAAACATAAAA";

  align::SmithWatermanT<char, short, width, 16, 10> sw(similarityScores, gapInit, gapExtend, 5);
  std::string                                       result;
  sw.align(q.data(), q.data() + q.size(), reference, reference + sizeof(reference) - 1, 11, 3, false, result);
  ASSERT_EQ(expected, result);
}

TEST(Alignments, insAtStart)
{
  using namespace dragenos;

  constexpr int                           MATCH    = 1;
  constexpr int                           MISMATCH = -4;
  const dragenos::align::SimilarityScores similarityScores(MATCH, MISMATCH);
  constexpr int                           gapExtend = 1;
  constexpr int                           gapInit   = 7;

  const char reference[] =
      "G"
      "TGGGCCGATTAAAAAAAAAAAAAAAACATAAAAAAAATCATCTCTACCCCAAAGTGGTACAATCGGTCCTACC";
  const std::string expected = "NSSMMMMMMMMMMDMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM";
  // reverse                     TATGGGCCGATT AAAAAAAAAAAAAAACATAAAAAAAATCATCTCTACCCCAAAGTGGTACAATCGGT
  const std::string q = "TGGCTAACATGGTGAAACCCCATCTCTACTAAAAAAAATACAAAAAAAAAAAAAAATTAGCCGGGTAT";

  align::SmithWatermanT<char, short, 48, 16, 10> sw(similarityScores, gapInit, gapExtend);
  std::string                                    result;
  sw.align(q.data(), q.data() + q.size(), reference, reference + sizeof(reference) - 1, 11, 2, true, result);
  ASSERT_EQ(expected, result);
}

TEST(Alignments, backtrackingEdge)
{
  using namespace dragenos;

  constexpr int                           MATCH    = 1;
  constexpr int                           MISMATCH = -4;
  const dragenos::align::SimilarityScores similarityScores(MATCH, MISMATCH);
  constexpr int                           gapExtend = 1;
  constexpr int                           gapInit   = 7;

  const char        reference[] = "GCCGGGCGTAGTGGCACA"     "CATCTGTAGTCCCAGATACTCAAGAGGCTGAGGCAGGAGA";
  const std::string expected    = "NNNNNNNNNNNNNNNNNNSSSSSSSMMMMMSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS";
//  const std::string expected    = "SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSMMMMMMSSS";
  const std::string q           =                   "GAAACCCCATCTCTACTAAAAAAAATACAAAAAAAAAAAAAAATTAGCCGGGTAT";

  align::SmithWatermanT<char, short, 48, 16, 10> sw(similarityScores, gapInit, gapExtend);
  std::string                                    result;
  sw.align(q.data(), q.data() + q.size(), reference, reference + sizeof(reference) - 1, 11, 0, false, result);
  ASSERT_EQ(expected, result);
}

TEST(Alignments, backtrackingEdge2)
{
  using namespace dragenos;

  constexpr int                           MATCH    = 1;
  constexpr int                           MISMATCH = -4;
  const dragenos::align::SimilarityScores similarityScores(MATCH, MISMATCH);
  constexpr int                           gapExtend = 1;
  constexpr int                           gapInit   = 7;

  const char reference[] =     "A"                                    "AAAAATACAAAAAAAAAAAAAAAAATTAGCCGGGCGTAGTGGCACACATCTGTAGTCCCAGATACTCAAGAGGCTGAGGCA";
//  const std::string expected = "SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSMMMMMMMMMMMMMMMMMMMMMMMMSSSSS";
  const std::string expected = "NSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSMMMMMMMMMMMMMMMMMMMMMMMSSSSS";
  const std::string q        =  "CATCCTGGCTAACATGGTGAAACCCCATCTCTACTAAAAAAAATACAAAAAAAAAAAAAAATTAGC";

  align::SmithWatermanT<char, short, 48, 16, 10> sw(similarityScores, gapInit, gapExtend);
  std::string                                    result;
  sw.align(q.data(), q.data() + q.size(), reference, reference + sizeof(reference) - 1, 11, 2, false, result);
  ASSERT_EQ(expected, result);
}

TEST(Alignments, horizontal3Forward)
{
  using namespace dragenos;

  constexpr int                           MATCH    = 1;
  constexpr int                           MISMATCH = -4;
  const dragenos::align::SimilarityScores similarityScores(MATCH, MISMATCH);
  constexpr int                           gapExtend = 1;
  constexpr int                           gapInit   = 6;

  const int  width = 48;
  const char reference[] =
      "TAGTGGTGTGGGCCGATTA"
      "AAAAAAAAAAAAAAACATAAAAAAAATCATCTCTACCCCAAAGTGG";
  const std::string expected =
      "NNNNNNMMMMMMMMMMMMD"
      "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM";
  const std::string q =
      "TATGGGCCGATT"
      "AAAAAAAAAAAAAAACATAAAAAAAATCATCTCTACCCCAAA";

  align::SmithWatermanT<char, short, width, 16, 10> sw(similarityScores, gapInit, gapExtend, 5);
  std::string                                       result;
  sw.align(q.data(), q.data() + q.size(), reference, reference + sizeof(reference) - 1, 11, 3, false, result);
  ASSERT_EQ(expected, result);
}

// TEST(Alignments, unnecessaryIndels)
//{
//  using namespace dragenos;
//
//  constexpr int MATCH = 1;
//  constexpr int MISMATCH = -4;
//  const dragenos::align::SimilarityScores similarityScores(MATCH, MISMATCH);
//  constexpr int gapExtend = 1;
//  constexpr int gapInit = 7;
//
//  const int width = 8;
////  const char reference[] =     "CTGGGATTTTATAGGACCAGTGTTGAGAAAATACAACCTGGAAACCCTTGGTATATAAGGTAAATTATGGCAGGGCAAGCTGCATAAGGAAAGAAACCAAAAGTGAAAATTCAATATTGGAGTCATTTTTTTATCTCTTACCTCTTCTTCTTCCTGTGCTTCTACTGCTGGTCCATTATGTGCCTCCTGGAAAAGGAGTTTCTTCATCTTTCGATACTGCAGATTGTCCAGCTCTCTTACTGCATC";
//  const char reference[] =     "GCATAAGGAAAGAAACCAAAAGTGAAAATTCAATATTGGAGTCATTTTT";
////  const char reference[] =     "CATAAGGAAAGAAACCAAAAGTGAAAATTCAATATTGGAGTCATTT";
//  const std::string expected = "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM";
////  const std::string q =        "CATAAGGGAAAAAAACAAAAAAGGAAAATTAATATTTGAGGCATTT";
////  const std::string q =        "CATAAGGGAAAAAAACAAAAAAGGAAAATT";
//  const std::string q =        "GCATAAGGGAAAAAAACAAAAAAGGAAAATTAATATTTGAGGCATTTTT";
//
////  const std::string q =              "TTGGTATATAAGGTAAATTATGGCAGGGGAAGCTGCATAAGGGAAAAAAACAAAAAAGGAAAATTAATATTTGAGGCATTTTTTTTTTTCTTTCCTCTTCTTCTTCCTGTGCTTCTACTGCTGGGCCCTTATTTGCCTCCTGGGAAAAGA";
//
//  align::SmithWatermanT<char, short, width, 16, 10> sw(similarityScores, gapInit, gapExtend, 5);
//  const auto result = sw.align(q.data(), q.data() + q.size(), reference, reference + sizeof(reference) - 1, 11, 6, false);
//  ASSERT_EQ(expected, result);
//}

// TEST(Alignments, unnecessaryIndels)
//{
//  using namespace dragenos;
//
//  constexpr int MATCH = 1;
//  constexpr int MISMATCH = -4;
//  const dragenos::align::SimilarityScores similarityScores(MATCH, MISMATCH);
//  constexpr int gapExtend = 1;
//  constexpr int gapInit = 6;
//
//  const int width = 48;
////  const char reference[] =     "CTGGGATTTTATAGGACCAGTGTTGAGAAAATACAACCTGGAAACCCTTGGTATATAAGGTAAATTATGGCAGGGCAAGCTGCATAAGGAAAGAAACCAAAAGTGAAAATTCAATATTGGAGTCATTTTTTTATCTCTTACCTCTTCTTCTTCCTGTGCTTCTACTGCTGGTCCATTATGTGCCTCCTGGAAAAGGAGTTTCTTCATCTTTCGATACTGCAGATTGTCCAGCTCTCTTACTGCATC";
//  const char reference[] =     "ATGGCAGGGCAAGCTGCATAAGGAAAGAAACCAAAAGTGAAAATTCAATATTGG";
//  const std::string expected = "NNNNNNMMMMMMMMMMMMD""MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM";
//  const std::string q =        "ATGGCAGGGGAAGCTGCATAAGGGAAAAAAACAAAAAAGGAAAATTAATAT";
////  const std::string q =              "TTGGTATATAAGGTAAATTATGGCAGGGGAAGCTGCATAAGGGAAAAAAACAAAAAAGGAAAATTAATATTTGAGGCATTTTTTTTTTTCTTTCCTCTTCTTCTTCCTGTGCTTCTACTGCTGGGCCCTTATTTGCCTCCTGGGAAAAGA";
//
//  align::SmithWatermanT<char, short, width, 16, 10> sw(similarityScores, gapInit, gapExtend, 5);
//  const auto result = sw.align(q.data(), q.data() + q.size(), reference, reference + sizeof(reference) - 1, 11, 3, false);
//  ASSERT_EQ(expected, result);
//}

//                                               MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMIMMMMMMMMMMMMMMMMMMMMMMDMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN

//
// TEST(Alignments, horizontal3Reverse)
//{
//  using namespace dragenos;
//
//  constexpr int MATCH = 1;
//  constexpr int MISMATCH = -4;
//  const dragenos::align::SimilarityScores similarityScores(MATCH, MISMATCH);
//  constexpr int gapExtend = 1;
//  constexpr int gapInit = 7;
//
//  const int width = 48;
//  const char reference[] =     "TCGCCTCCCTGGAACACCATTTGGTAACTTATGAGGCATAACCCTGTTCAGGCTCCCAGGGCTATTATGCACATTTTCTAAAATTTCAGGCATGTTGATCTTTGCACTGTGATTACTTTTTCATCAAAAGCCACACAGAGGGATGTGGAGTGACCGTAATGTGAGTGCTGCTGGGGCAGGGGGTACCGGCCATCCCGGAGGTGTGAGGGGCAGGTACCTGGAGCCTGGCTTCTGGCTACACCGGGC";
//  const std::string expected =   "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNMMMMMMMMMMMMDMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM";
//  const std::string q =                                             "CAGCTCCCAGGGCTATTATGCACATTTTCTAAAATTTCAGGCATGTTGATCTTTGCACTGTGATTACTTTTTCATCAAAAGCCACACAGAGGGATGTGGAGTGACCGTAATGTGAGTGCTGCTGGGGCAGGGGGTACCGGCCATCCCGGA";
//
//  align::SmithWatermanT<char, short, width, 16, 10> sw(similarityScores, gapInit, gapExtend, 5);
//  const auto result = sw.align(q.data(), q.data() + q.size(), reference, reference + sizeof(reference) - 1, 11, width, false);
//  ASSERT_EQ(expected, result);
//}

TEST(Alignments, ambigDeletion)
{
  using namespace dragenos;

  constexpr int                           MATCH    = 1;
  constexpr int                           MISMATCH = -4;
  const dragenos::align::SimilarityScores similarityScores(MATCH, MISMATCH);
  constexpr int                           gapExtend = 1;
  constexpr int                           gapInit   = 7;

  const int  width = 48;
  const char reference[] =
      "TTGTCGAGGCTGTACAAAGCATTCACTCTGTTCGGTCACGTTCAAAAAAAAAAAGGAAAAAAAAAAAAACAGAAAACGAATGGAAGAACGAATTACCTTAACAATACCGATTCGTGTATCTTCCGGTTTTTTCCTCAAAAGGTTTGGGTCGTTTAGTTCACGAACCTAAGACTTGACGGTTTTCTTTTGACGTGAAGGGGAGAATTCATTTTGCTTTACTCAAAGAATCCATTTACATAAGTAGTC";
  const std::string expected =
      "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNMMMMMMMMDDMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM";
  const std::string q =
      "AAAAAAAG"
      "AAAAAAAAAAAACAGAAAACGAATGGAAGAACGAATTACCTTAACAATACCGATTCGTGTATCTTCCGGTTTTTTCCTCAAAAGGTTTGGGTCGTTTAGTTCACGAACCTAAGACTTGACGGTTTTCTTTTGACGTGAAGGG";

  align::SmithWatermanT<char, short, width, 16, 10> sw(similarityScores, gapInit, gapExtend, 5);
  std::string                                       result;
  sw.align(
      q.data(), q.data() + q.size(), reference, reference + sizeof(reference) - 1, 11, width, false, result);
  ASSERT_EQ(expected, result);
}

TEST(Alignments, ambigDeletionRev)
{
  using namespace dragenos;

  constexpr int                           MATCH    = 1;
  constexpr int                           MISMATCH = -4;
  const dragenos::align::SimilarityScores similarityScores(MATCH, MISMATCH);
  constexpr int                           gapExtend = 1;
  constexpr int                           gapInit   = 7;

  const int  width = 48;
  const char reference[] =
      "CTGATGAATACATTTACCTAAGAAACTCATTTCGTTTTACTTAAGAGGGGAAGTGCAGTTTTCTTTTGGCAGTTCAGAATCCAAGCACTTGATTTGCTGGGTTTGGAAAACTCCTTTTTTGGCCTTCTATGTGCTTAGCCATAACAATTCCATTAAGCAAGAAGGTAAGCAAAAGACAAAAAAAAAAAAAGGAAAAAAAAAAACTTGCACTGGCTTGTCTCACTTACGAAACATGTCGGAGCTGTT";
  const std::string expected =
      "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMDDMMMMMMMM";
  const std::string q =
      "GGGAAGTGCAGTTTTCTTTTGGCAGTTCAGAATCCAAGCACTTGATTTGCTGGGTTTGGAAAACTCCTTTTTTGGCCTTCTATGTGCTTAGCCATAACAATTCCATTAAGCAAGAAGGTAAGCAAAAGACAAAAAAAAAAAAGAAAAAAA";

  align::SmithWatermanT<char, short, width, 16, 10> sw(similarityScores, gapInit, gapExtend, 5);
  std::string                                       result;
  sw.align(
      q.data(), q.data() + q.size(), reference, reference + sizeof(reference) - 1, 11, width, false, result);
  ASSERT_EQ(expected, result);
}

TEST(Alignments, ambigDeletionShort)
{
  using namespace dragenos;

  constexpr int                           MATCH    = 1;
  constexpr int                           MISMATCH = -4;
  const dragenos::align::SimilarityScores similarityScores(MATCH, MISMATCH);
  constexpr int                           gapExtend = 1;
  constexpr int                           gapInit   = 7;

  const int         width       = 8;
  const char        reference[] = "AAAAAAAGGA""AAAAAAAAAAAACAG";
  const std::string expected    = "NMMMMMMMMDMMMMMMMMMMMMMMM";
//  const std::string expected    = "MMMMMMMMDDMMMMMMMMMMMMMMM";
  const std::string q =            "AAAAAAAG" "AAAAAAAAAAAACAG";

  align::SmithWatermanT<char, short, width, 16, 10> sw(similarityScores, gapInit, gapExtend, 5);
  std::string                                       result;
  sw.align(q.data(), q.data() + q.size(), reference, reference + sizeof(reference) - 1, 11, 5, false, result);
  ASSERT_EQ(expected, result);
}

TEST(Alignments, ambigDeletionRevShort)
{
  using namespace dragenos;

  constexpr int                           MATCH    = 1;
  constexpr int                           MISMATCH = -4;
  const dragenos::align::SimilarityScores similarityScores(MATCH, MISMATCH);
  constexpr int                           gapExtend = 1;
  constexpr int                           gapInit   = 7;

  const int         width       = 8;
  const char        reference[] = "GGGACAAAAAAAAAAAAAGGAAAAAAA";
  const std::string expected    = "MMMMMMMMMMMMMMMMMDDMMMMMMMM";
  const std::string q =           "GGGACAAAAAAAAAAAA""GAAAAAAA";

  align::SmithWatermanT<char, short, width, 16, 10> sw(similarityScores, gapInit, gapExtend, 5);
  std::string                                       result;
  sw.align(q.data(), q.data() + q.size(), reference, reference + sizeof(reference) - 1, 11, 5, false, result);
  ASSERT_EQ(expected, result);
}

TEST(Alignments, dbSequenceOverrunRandomCigarFix)
{
  using namespace dragenos;

  constexpr int                           MATCH    = 1;
  constexpr int                           MISMATCH = -4;
  const dragenos::align::SimilarityScores similarityScores(MATCH, MISMATCH);
  constexpr int                           gapExtend = 1;
  constexpr int                           gapInit   = 7;

  const int width = 48;

//  const char reference[] =
//      "GTCCAGCCCCAGGGCAGAGGTCAGAAGTAAGCGTCCTGCCGGGCCTCGGCTGAGGGCCTGTCCCCACCATCATGGGGGCCTGGGGGCCCATCTGCCTTTCGACGCAGCGAAAGCTTCCGGCCTTGAGCCTTCTTCTCCTGCTTCTGCCGCTCCTTCTCCTGCTTCTCTCGCTCCTTTTCCTGCTTCTCTC"
//      "GCTCTTTCTCCTGCTTCTGCCGCTCCTTCTCCCGCTCCTTCTCCTGTTTCTGCC";
  const char reference[] =
      "GTCCAGCCCCAGGGCAGAGGTCAGAAGTAAGCGTCCTGCCGGGCCTCGGCTGAGGGCCTGTCCCCACCATCATGGGGGCCTGGGGGCCCATCTGCCTTTCGACGCAGCGAAAGCTTCCGGCCTTGAGCCTTCTTCTCCTGCTTCTGCCGCTCCTTCTCCTGCTTCTCTCGCTCCTTTTCCTGCTTCTCTCGC"
      "TCTTTCTCCTGCTTCTGCCGCTCCTTCTCCCGCTCCTTCTCCTGTTTCTGCC";
//  const std::string expected =
//      "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
//      "SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS"
//      "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMSSSSSSSSSSSSSSSS";
//  const std::string expected =
//      "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
//      "SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS"
//      "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMSSSSSSSSSSSSSSSS";

  const std::string expected =
      "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
      "SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS"
      "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMSSSSSSSSSSSSSSSS";
//  const std::string expected =
//      "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
//      "SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS"
//      "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMSSSSSSSSSSSSSSSS";

  const std::string q =
      "GCTTCCGGCCTTGAGCCTTCTTCTCCTGCTTCTGCCGCTACTTCTCCTGCTTCTCTCGCTCCTTTTCCTGCTTCTCTCGC"
      "TCTTTCTCCTGCTTCTGCCGCTCCTTCTCCCGCTCCTTCTCCTGTTTCTGCCGCTCCTTCTCCTGTTT";

  align::SmithWatermanT<char, short, width, 16, 11> sw(similarityScores, gapInit, gapExtend, 5);
  static constexpr size_t                           forcedDiagonalMotion   = 2;
  static constexpr size_t                           forcedHorizontalMotion = sw.width + 1;
  std::string                                       result;
  sw.align(
      q.data(),
      q.data() + q.size(),
      reference,
      reference + sizeof(reference) - 1,
      forcedDiagonalMotion,
      forcedHorizontalMotion,
      false,
      result);
  ASSERT_EQ(expected, result);
}

TEST(Alignments, dbSequenceOverrunRandomCigarFix2)
{
  using namespace dragenos;

  constexpr int                           MATCH    = 1;
  constexpr int                           MISMATCH = -4;
  const dragenos::align::SimilarityScores similarityScores(MATCH, MISMATCH);
  constexpr int                           gapExtend = 1;
  constexpr int                           gapInit   = 7;

  const int width = 48;

  const char reference[] =
      "TGACCTGTACGACGACCTGAACGGTGACCTGAACTGTGACCTGGACGGGGG"
      "CCTGGACGGCGACCTGGACGGTGACCTGGACGGTGACCTGGACTGCGACCTGGAC"
      "TACGACCTGGAGAGGGACCTGGACTACGACCTGGACGGCGACCTGAACTGG"
      "GACCTGGACGGGGACCTGGACGGCGACCTGGACTACGACCTGGACAGGGACCTGGACGGGGACCTGGACGGTGACCTGGACAGGGACCT";
  const std::string expected =
      "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
      "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
      "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
      "SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS"
      "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMDDDDDDDDDDDDMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMS";
  const std::string q =
      "AACGGGAGCTGGGCCCCGGACGGGGCCCAGGACGGAAACCGGCACAGAGACCTGGACGGGTACCGGGACGGC"
      "GACCTGGACGGGGACCTGGACAGCGCCCTGG"
      "ACAGGGACCTGGACGGGGACCTGGACGGCGACCTGGACAGGGACCTG";

  align::SmithWatermanT<char, short, width, 16, 11> sw(similarityScores, gapInit, gapExtend, 5);
  static constexpr size_t                           forcedDiagonalMotion   = 2;
  static constexpr size_t                           forcedHorizontalMotion = sw.width + 1;
  std::string                                       result;
  sw.align(
      q.data(),
      q.data() + q.size(),
      reference,
      reference + sizeof(reference) - 1,
      forcedDiagonalMotion,
      forcedHorizontalMotion,
      false,
      result);
  ASSERT_EQ(expected, result);
}

TEST(Alignments, dbSequenceOverrunRandomCigarFix3)
{
  using namespace dragenos;

  constexpr int                           MATCH    = 1;
  constexpr int                           MISMATCH = -4;
  const dragenos::align::SimilarityScores similarityScores(MATCH, MISMATCH);
  constexpr int                           gapExtend = 1;
  constexpr int                           gapInit   = 7;

  const int width = 48;

  const char        reference[] = "GGGGGGGGGGGAGGGAGGGGGGGGTGAGAGGGGGGGGGGGGGGGGGGGG";
  const std::string expected    = "NSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS";
  const std::string q           = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";

  align::SmithWatermanT<char, short, width, 16, 11> sw(similarityScores, gapInit, gapExtend, 5);
  static constexpr size_t                           forcedDiagonalMotion   = 2;
  static constexpr size_t                           forcedHorizontalMotion = 1;
  std::string                                       result;
  sw.align(
      q.data(),
      q.data() + q.size(),
      reference,
      reference + sizeof(reference) - 1,
      forcedDiagonalMotion,
      forcedHorizontalMotion,
      false,
      result);
  ASSERT_EQ(expected, result);
}

// Another couple of tests to highlight unclipScore behavior
TEST(Alignments, 94M2D56MInsteadOf100M50S)
{
  using namespace dragenos;

  constexpr int                           MATCH    = 1;
  constexpr int                           MISMATCH = -4;
  const dragenos::align::SimilarityScores similarityScores(MATCH, MISMATCH);
  constexpr int                           gapExtend = 1;
  constexpr int                           gapInit   = 7;

  const int width = 48;

  const char reference[] =
      "TGACCCATTAAATACTTCGTTTCTCCAAATTAACTGAGTGTCAAGACGTCCGAATTGTCCTTCGTACTGACCCTCCGGAGTTCTTTGAATGTTAGTACCGTCTTCCACTTCCCCTTCGTTCGTGGAAGAAGTGTACCACCGTCTTCTCTGTCACTTCCCTCTCCACGGTGTGTGAAAATTTCCTAGTCTAGAGGACTCTTGAGTGAGTGAAAATGTTCTTGTCGTTCCCTTTATAGGCGGGGTACT";
  const std::string expected =
      "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMDDMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM";
  //  const std::string q = "GTTCTCAGGAGATCTGGTGCTTTAAAAGTGTGTGTGACTTTGACCTTCTCTCTCTCCTGCCACCCTGTGAAGAAGGTGCTTGCTTCCCCTTCACCTTGTGCTATGATTGTAAGTTTCCGGAGGACTCCCAATCATGCTTCCTGTTAAGCA";
  const std::string q =
      "ACGAATTGTCCTTCGTACTAACCCTCAGGAGGCCTTTGAATGTTAGTATCGTGTTCCACTTCCCCTTCGTTCGTGGAAGAAGTGTCCCACCGTCCTCTCTCTCTTCCAGTTTCAGTGTGTGTGAAAATTTCGTGGTCTAGAGGACTCTTG";

  align::SmithWatermanT<char, short, width, 16, 11> sw(similarityScores, gapInit, gapExtend, 5);
  static constexpr size_t                           forcedDiagonalMotion   = 2;
  static constexpr size_t                           forcedHorizontalMotion = sw.width + 1;
  std::string                                       result;
  sw.align(
      q.data(),
      q.data() + q.size(),
      reference,
      reference + sizeof(reference) - 1,
      forcedDiagonalMotion,
      forcedHorizontalMotion,
      false,
      result);
  ASSERT_EQ(expected, result);
}

TEST(Alignments, 100M50SInsteadOf94M2D56M)
{
  using namespace dragenos;

  constexpr int                           MATCH    = 1;
  constexpr int                           MISMATCH = -4;
  const dragenos::align::SimilarityScores similarityScores(MATCH, MISMATCH);
  constexpr int                           gapExtend = 1;
  constexpr int                           gapInit   = 7;

  const int width = 48;

  const char reference[] =
      "TGACCCATTAAATACTTCGTTTCTCCAAATTAACTGAGTGTCAAGACGTCCGAATTGTCCTTCGTACTGACCCTCCGGAGTTCTTTGAATGTTAGTACCGTCTTCCACTTCCCCTTCGTTCGTGGAAGAAGTGTACCACCGTCTTCTCTGTCACTTCCCTCTCCACGGTGTGTGAAAATTTCCTAGTCTAGAGGACTCTTGAGTGAGTGAAAATGTTCTTGTCGTTCCCTTTATAGGCGGGGTACT";
  const std::string expected =
      "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNSMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS";
  //  const std::string q = "GTTCTCAGGAGATCTGGTGCTTTAAAAGTGTGTGTGACTTTGACCTTCTCTCTCTCCTGCCACCCTGTGAAGAAGGTGCTTGCTTCCCCTTCACCTTGTGCTATGATTGTAAGTTTCCGGAGGACTCCCAATCATGCTTCCTGTTAAGCA";
  const std::string q =
      "ACGAATTGTCCTTCGTACTAACCCTCAGGAGGCCTTTGAATGTTAGTATCGTGTTCCACTTCCCCTTCGTTCGTGGAAGAAGTGTCCCACCGTCCTCTCTCTCTTCCAGTTTCAGTGTGTGTGAAAATTTCGTGGTCTAGAGGACTCTTG";

  align::SmithWatermanT<char, short, width, 16, 11> sw(similarityScores, gapInit, gapExtend, 0);
  static constexpr size_t                           forcedDiagonalMotion   = 2;
  static constexpr size_t                           forcedHorizontalMotion = sw.width + 1;
  std::string                                       result;
  sw.align(
      q.data(),
      q.data() + q.size(),
      reference,
      reference + sizeof(reference) - 1,
      forcedDiagonalMotion,
      forcedHorizontalMotion,
      false,
      result);
  ASSERT_EQ(expected, result);
}


TEST(Alignments, fixForMInsteadOfSAtEnd)
{
  using namespace dragenos;

  constexpr int                           MATCH    = 1;
  constexpr int                           MISMATCH = -4;
  const dragenos::align::SimilarityScores similarityScores(MATCH, MISMATCH);
  constexpr int                           gapExtend = 1;
  constexpr int                           gapInit   = 7;

//  const char reference[] =
//      "TCACTACTTGTATATGGTGTGTCTCCTCTTCGTGAGGAGTTAAGTGTGACACTACTTGTATATGGTGTGTCTCCTCTTCGTGAGGAGTTAAGTGTGACACTACTTGTATATGGTGTGTCTCCTCTTCGTGAGGAAGTTAA";
  const char reference[] =     "GTGACACTACTTGTATATGGTGTGTCTCCTCTTCGTGAGGAAGTTAA";
  const std::string expected = "NNNNMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMS";
//  const std::string expected =
//      "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMS";
  const std::string q =            "CACTACTTGTATATGGTGTGTCTCCTCTTCGTGAGGAAGTTAAG";
//  const std::string q = "GAATTGAAGGAGTGCTTCTCCTCTGTGTGGTATATGTTCATCAC";

  align::SmithWatermanT<char, short, 48, 16, 11> sw(similarityScores, gapInit, gapExtend, 5);
  static constexpr size_t                           forcedDiagonalMotion   = 2;
  static constexpr size_t                           forcedHorizontalMotion = sw.width + 1;
  std::string                                       result;
  sw.align(
      q.data(),
      q.data() + q.size(),
      reference,
      reference + sizeof(reference) - 1,
      forcedDiagonalMotion,
      forcedHorizontalMotion,
      false,
      result);
  ASSERT_EQ(expected, result);
}

TEST(Alignments, fixForWrongSoftClipAtEnd)
{
  using namespace dragenos;

  constexpr int                           MATCH    = 1;
  constexpr int                           MISMATCH = -4;
  const dragenos::align::SimilarityScores similarityScores(MATCH, MISMATCH);
  constexpr int                           gapExtend = 1;
  constexpr int                           gapInit   = 7;

//  const char reference[] =     "CCCCGACAAGGGGAGAAACCCCTTGGACATCCCTCACGACTCCGCCGTACCAAGACTCAGTGTCCCCTGGACTCCTGTGTCCCTACCCCGTACCACTCGTCGGAGACACACCCTGTGCCCCTGCCCCGTCCCACCCGTCGGAGATCCACCCTGTACCCCTGCCCCGTCCCACCCGTCGGAGATGCACCC"
//                               "TGTGCCCCTACCCCGTCCCACCCGTCGGAGACGACCCTGTGCCCCTGCCCCGTCCCACC";
//  const std::string expected = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
//                               "SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS"
//                               "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM";
//  //  const std::string q =        "GTCCCCGTGTCCCAGCAGAGGCTGCCCACCCTGCCCCATCCCCGTGTCCCCCGGAGAGGGTGCCCACCCCCCCCCCTCCCCCAGTCCCCCCCAAAGGCTTGCCACCCCCCCCCGCCCCCCGGTCCCCCCCAACAGGTCGCCCCAAAACC";
//  const std::string q =        "CCAAAACCCCGCTGGACAACCCCCCCTGGCCCCCCGCCCCCCCCCACCGTTCGGAAACCCCCCCTGACCCCCTCCCCCCCCCCACCCGTGGGAGAGGCCCCC"
//                               "TGTGCCCCTACCCCGTCCCACCCGTCGGAGACGACCCTGTGCCCCTG";


  const char reference[] =     "CCCCGACAAGGGGAGAAACCCCTTGGACATCCCTCACGACTCCGCCGTACCAAGACTCAGTGTCCCCTGGACTCCTGTGTCCCTACCCCGTACCACTCGTCGGAGACACACCCTGTGCCCCTGCCCCGTCCCACCCGTCGGAGATCCACCCTGTACCCCTGCCCCGT"
                               "CCCACCCGTCGGAGATGCACCCTGTGCCCCTACCCCGTCCCACCCGTCGGAGACGACCCTGTGCCCCTGCCCCGTCCCACC";

  const std::string expected = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
                               "SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS"
                               "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM";
  const std::string q =        "CCAAAACCCCGCTGGACAACCCCCCCTGGCCCCCCGCCCCCCCCCACCGTTCGGAAACCCCCCCTGACCCCCTCCCCCCC"
                               "CCCACCCGTGGGAGAGGCCCCCTGTGCCCCTACCCCGTCCCACCCGTCGGAGACGACCCTGTGCCCCTG";

  align::SmithWatermanT<char, short, 48, 16, 9> sw(similarityScores, gapInit, gapExtend, 5);
  static constexpr size_t                           forcedDiagonalMotion   = 99;
  static constexpr size_t                           forcedHorizontalMotion = sw.width;
  std::string                                       result;
  sw.align(
      q.data(),
      q.data() + q.size(),
      reference,
      reference + sizeof(reference) - 1,
      forcedDiagonalMotion,
      forcedHorizontalMotion,
      false,
      result);
  ASSERT_EQ(expected, result);
}

