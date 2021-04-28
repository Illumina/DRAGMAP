#include "gtest/gtest.h"

#include <algorithm>

#include "align/SmithWaterman.hpp"

TEST(SmithWaterman, Reset)
{
  using namespace dragenos::align;
  typedef SmithWatermanT<unsigned char, short, 48, 16, 11> SmithWaterman;
  constexpr auto                                           width      = SmithWaterman::width;
  const auto                                               NOT_A_BASE = SmithWaterman::NOT_A_QUERY_BASE;
  constexpr int                                            MATCH      = 2;
  constexpr int                                            MISMATCH   = -3;
  constexpr int                                            gapExtend  = 4;
  constexpr int                                            gapInit    = 5;
  const SimilarityScores                                   similarity(MATCH, MISMATCH);
  const auto random = std::string("ACGTTGAGGTTCCGTAGTATGACCTGTTTTAACGTTAGGCTGGAAAGTNNC");
  ASSERT_GE(random.size(), width);
  const Query    query(random.begin(), random.end());
  SmithWaterman  sw(similarity, gapInit, gapExtend);
  const Database database(random.begin(), random.end());
  sw.reset(
      query.data(), query.data() + query.size(), database.data(), database.data() + database.size(), false);
  const std::vector<unsigned char>& reversedQuery = sw.getReversedQuery();
  ASSERT_EQ(random.size() + 2 * width - 2, reversedQuery.size());
  ASSERT_EQ(reversedQuery.begin() - 1 + random.size(), sw.getQueryIt());
}

TEST(SmithWaterman, MoveDown)
{
  using namespace dragenos::align;
  typedef SmithWatermanT<unsigned char, short, 48, 16, 11> SmithWaterman;
  typedef SmithWaterman::Motion                            Motion;
  constexpr auto                                           width     = SmithWaterman::width;
  constexpr int                                            MATCH     = 3;
  constexpr int                                            MISMATCH  = -4;
  constexpr int                                            gapExtend = 2;
  constexpr int                                            gapInit   = 5;
  const SimilarityScores                                   similarity(MATCH, MISMATCH);
  const auto random = std::string("ACGTTGAGGTTCCGTAGTATGACCTGTTTTAACGTTAGGCTGGAAAGTNNC");
  //                               CNNTGAAAGGTCGGATTGCAATTTTGTCCAGTATGATGCCTTGGAGTTGCA
  ASSERT_GE(random.size(), width);
  const Query    query(random.begin(), random.end());
  const Database database(random.begin(), random.end());
  SmithWaterman  sw(similarity, gapInit, gapExtend);
  sw.reset(
      query.data(), query.data() + query.size(), database.data(), database.data() + database.size(), false);
  ASSERT_EQ(sw.getMotions().size(), 0);
  ASSERT_EQ(sw.getScores().size(), 0);
  sw.moveDown();
  ASSERT_EQ(sw.getMotions().size(), 1);
  ASSERT_EQ(sw.getScores().size(), 1);
  // first base matches, but first column is all 0 in HW
  ASSERT_EQ(sw.getScores().back()[0], 0);
  // all other have NOT_A_BASE on the query
  for (size_t i = 1; sw.getScores().back().size() > i; ++i)
    ASSERT_EQ(sw.getScores().back()[i], 0) << "i == " << i;
  ASSERT_EQ(sw.getMotions().back(), Motion::down);
  sw.moveDown();
  ASSERT_EQ(sw.getMotions().size(), 2);
  ASSERT_EQ(sw.getScores().size(), 2);
  // All mismatches (AC vs CA)
  for (size_t i = 0; sw.getScores().back().size() > i; ++i)
    ASSERT_EQ(sw.getScores().back()[i], 0) << "i == " << i;
  ASSERT_EQ(sw.getMotions().back(), Motion::down);
  sw.moveDown();
  ASSERT_EQ(sw.getMotions().size(), 3);
  ASSERT_EQ(sw.getScores().size(), 3);
  // ACG vs GCA: C matches
  ASSERT_EQ(sw.getScores().back()[0], 0);
  ASSERT_EQ(sw.getMotions().back(), Motion::down);
}

TEST(SmithWaterman, MoveRight)
{
  using namespace dragenos::align;
  typedef SmithWatermanT<unsigned char, short, 48, 16, 11> SmithWaterman;
  typedef SmithWaterman::Motion                            Motion;
  constexpr auto                                           width     = SmithWaterman::width;
  constexpr int                                            MATCH     = 3;
  constexpr int                                            MISMATCH  = -4;
  constexpr int                                            gapExtend = 2;
  constexpr int                                            gapInit   = 5;
  const SimilarityScores                                   similarity(MATCH, MISMATCH);
  const auto random = std::string("ACGTTGAGGTTCCGTAGTATGACCTGTTTTAACGTTAGGCTGGAAAGT") + "AATTCCAGNNC";
  ASSERT_GE(random.size(), width);
  const Query    query(random.begin(), random.begin() + width);
  const Database database(random.begin(), random.end());
  SmithWaterman  sw(similarity, gapInit, gapExtend);
  sw.reset(
      query.data(), query.data() + query.size(), database.data(), database.data() + database.size(), false);
  ASSERT_EQ(sw.getMotions().size(), 0);
  ASSERT_EQ(sw.getScores().size(), 0);
  // move down into the actual query string
  for (unsigned i = 0; width > i; ++i) sw.moveRight();
  ASSERT_EQ(sw.getMotions().size(), width);
  ASSERT_EQ(sw.getScores().size(), width);
}

TEST(SmithWaterman, PeakPosition)
{
  //FAIL() << "Not implemented";
}

TEST(SmithWaterman, NextMove)
{
  //FAIL() << "Not implemented";
}

TEST(SmithWaterman, NoSimilarity)
{
  using namespace dragenos::align;
  typedef SmithWatermanT<unsigned char, short, 48, 16, 11> SmithWaterman;
  constexpr auto                                           width     = SmithWaterman::width;
  constexpr int                                            MATCH     = 2;
  constexpr int                                            MISMATCH  = -3;
  constexpr int                                            gapExtend = 4;
  constexpr int                                            gapInit   = 5;
  const SimilarityScores                                   similarity(MATCH, MISMATCH);
  const auto                                               allC = std::vector<char>(48, 'c');
  const auto                                               allA = std::vector<char>(48, 'a');
  const Query                                              queryC(allC.begin(), allC.end());
  const Database                                           databaseA(allA.begin(), allA.end());
  SmithWaterman                                            sw(similarity, gapInit, gapExtend);
  sw.reset(
      queryC.data(),
      queryC.data() + queryC.size(),
      databaseA.data(),
      databaseA.data() + databaseA.size(),
      true);
  for (unsigned i = 0; width > i; ++i) {
    sw.moveRight();
  }

  const auto& s = sw.getSimilarities();
  for (unsigned i = 1; width > i; ++i) {
    ASSERT_EQ(MISMATCH, s[i]);
  }
}

TEST(SmithWaterman, AllSimilar)
{
  using namespace dragenos::align;
  typedef SmithWatermanT<unsigned char, short, 48, 16, 11> SmithWaterman;
  constexpr auto                                           width    = SmithWaterman::width;
  constexpr int                                            MATCH    = 2;
  constexpr int                                            MISMATCH = -3;
  const SimilarityScores                                   similarity(MATCH, MISMATCH);
  const auto random = std::string("ACGTTGAGGTTCCGTAGTATGACCTGTTTTAACGTTAGGCTGGAAAGT");
  ASSERT_EQ(random.size(), width);
  const Query    query(random.begin(), random.end());
  const Database database(random.begin(), random.end());
  SmithWaterman  sw(similarity, 5, 4);

  sw.reset(
      query.data(), query.data() + query.size(), database.data(), database.data() + database.size(), true);

  for (unsigned i = 0; width > i; ++i) {
    sw.moveRight();
  }

  const auto& s = sw.getSimilarities();
  for (unsigned i = 1; width > i; ++i) {
    ASSERT_EQ(MATCH, s[i]);
  }
}
