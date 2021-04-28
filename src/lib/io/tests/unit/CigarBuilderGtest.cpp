#include "gtest/gtest.h"
#include "inttypes.h"

#include "io/CigarBuilder.hpp"

TEST(CigarBuilderTest, Test1)
{
  dragenos::io::CigarBuilder cigar;
  for (int i = 0; i < 5; ++i) cigar.AddMatch();
  for (int i = 0; i < 15; ++i) cigar.AddInsertion();
  for (int i = 0; i < 25; ++i) cigar.AddDeletion();
  for (int i = 0; i < 35; ++i) cigar.AddSoftClip();
  cigar.ConsolidateRecords();

  std::ostringstream os;
  os << cigar;
  ASSERT_STREQ("5M15I25D35S", os.str().c_str());
}

TEST(CigarBuilderTest, LongSeries)
{
  dragenos::io::CigarBuilder cigar;
  for (int i = 0; i < 5; ++i) cigar.AddSoftClip();
  for (int i = 0; i < 130; ++i) cigar.AddInsertion();
  for (int i = 0; i < 130; ++i) cigar.AddMatch();
  for (int i = 0; i < 130; ++i) cigar.AddDeletion();
  for (int i = 0; i < 130; ++i) cigar.AddSoftClip();
  cigar.ConsolidateRecords();

  std::ostringstream os;
  os << cigar;
  ASSERT_STREQ("5S130I130M130D130S", os.str().c_str());
}

TEST(CigarBuilderTest, LongReverse)
{
  dragenos::io::CigarBuilder cigar;
  for (int i = 0; i < 5; ++i) cigar.AddSoftClip();
  for (int i = 0; i < 130; ++i) cigar.AddInsertion();
  for (int i = 0; i < 130; ++i) cigar.AddMatch();
  for (int i = 0; i < 130; ++i) cigar.AddDeletion();
  for (int i = 0; i < 130; ++i) cigar.AddSoftClip();

  cigar.Reverse();
  cigar.ConsolidateRecords();

  std::ostringstream os;
  os << cigar;
  ASSERT_STREQ("130S130D130M130I5S", os.str().c_str());
}
