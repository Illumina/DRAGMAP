#include "gtest/gtest.h"

#include "SeedPositionMocks.hpp"
#include "map/BestIntervalTracker.hpp"

using dragenos::map::BestIntervalTracker;
using dragenos::sequences::Read;
using dragenos::sequences::Seed;

TEST(BestIntervalTracker, isWorseThan)
{
  const Read     read;
  const unsigned primaryLength = 17;
  const unsigned readPosition  = 42;
  const Seed     seed(&read, readPosition, primaryLength);
  // less repeat, recursive
  BestIntervalTracker a(seed, 0, 16, 8);
  BestIntervalTracker b(seed, 5, 10, 10);
  ASSERT_EQ(b.isWorseThan(a), true);
  // more repeat, recursive
  BestIntervalTracker c(seed, 0, 1600, 8);
  BestIntervalTracker d(seed, 5, 1000, 10);
  ASSERT_EQ(d.isWorseThan(c), false);
  // less repeat, non-recursive
  BestIntervalTracker e(seed, 0, 16, 16);
  BestIntervalTracker f(seed, 5, 10, 10);
  ASSERT_EQ(f.isWorseThan(e), true);
  // more repeat, non-recursive
  BestIntervalTracker g(seed, 0, 1600, 16);
  BestIntervalTracker h(seed, 5, 1000, 10);
  ASSERT_EQ(h.isWorseThan(g), true);
}

TEST(BestIntervalTracker, isValidExtra)
{
  const Read     read;
  const auto     readLength    = read.getLength();
  const unsigned primaryLength = 17;
  const unsigned readPosition  = 42;
  ASSERT_GT(readLength, readPosition);
  const Seed seed(&read, readPosition, primaryLength);

  uint32_t longest_nonsample_seed_len = 0;
  uint32_t num_non_sample_seed_chains = 0;

  // too few num of chains activation
  BestIntervalTracker a(seed, 0, 16, 8);
  ASSERT_EQ(a.isValidExtra(num_non_sample_seed_chains, longest_nonsample_seed_len), true);

  // seed length activation
  num_non_sample_seed_chains = 10;
  BestIntervalTracker b(seed, 0, 16, 52);
  ASSERT_EQ(b.isValidExtra(num_non_sample_seed_chains, longest_nonsample_seed_len), true);

  // longest_nonsample_seed_len activation
  num_non_sample_seed_chains = 10;
  longest_nonsample_seed_len = 40;
  BestIntervalTracker c(seed, 0, 16, 52);
  ASSERT_EQ(c.isValidExtra(num_non_sample_seed_chains, longest_nonsample_seed_len), true);
}