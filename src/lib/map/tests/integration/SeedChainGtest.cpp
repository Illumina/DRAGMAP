#include "gtest/gtest.h"

#include <limits>

#include "map/SeedChain.hpp"

using dragenos::map::SeedChain;
using dragenos::map::SeedPosition;
using dragenos::sequences::Read;
using dragenos::sequences::Seed;
typedef SeedPosition::ReferencePosition ReferencePosition;
typedef Read::Name                      Name;
typedef Read::Bases                     Bases;
typedef Read::Qualities                 Qualities;

TEST(SeedChain, Constructor)
{
  SeedChain seedChain;
  // check that the seed chain is empty
  ASSERT_EQ(seedChain.begin(), seedChain.end());
  // check that it is not a reverse complement
  ASSERT_FALSE(seedChain.isReverseComplement());
  // check that it claims to have only random samples
  ASSERT_TRUE(seedChain.hasOnlyRandomSamples());
  Read           read;
  const uint64_t id       = 0;
  const unsigned position = 0;
  read.init(Name(), Bases(151), Qualities(151), id, position);
  const unsigned primaryLength = 17;
  const unsigned readPosition  = 42;
  const Seed     seed(&read, readPosition, primaryLength);
  // check that the seed chain accepts anything
  const unsigned halfExtension = 0;
  ASSERT_TRUE(seedChain.accepts(SeedPosition(seed, 0, halfExtension), false));
  ASSERT_TRUE(seedChain.accepts(SeedPosition(seed, 0, halfExtension), true));
  ASSERT_TRUE(seedChain.accepts(SeedPosition(seed, 100000, halfExtension), false));
  ASSERT_TRUE(seedChain.accepts(SeedPosition(seed, 100000, halfExtension), true));
  ASSERT_TRUE(seedChain.accepts(
      SeedPosition(seed, std::numeric_limits<ReferencePosition>::min(), halfExtension), false));
  ASSERT_TRUE(seedChain.accepts(
      SeedPosition(seed, std::numeric_limits<ReferencePosition>::min(), halfExtension), true));
  ASSERT_TRUE(seedChain.accepts(
      SeedPosition(seed, std::numeric_limits<ReferencePosition>::max(), halfExtension), false));
  ASSERT_TRUE(seedChain.accepts(
      SeedPosition(seed, std::numeric_limits<ReferencePosition>::max(), halfExtension), true));
}

TEST(SeedChain, SetReverseComplement)
{
  SeedChain seedChain;
  ASSERT_FALSE(seedChain.isReverseComplement());
  seedChain.setReverseComplement(true);
  ASSERT_TRUE(seedChain.isReverseComplement());
  seedChain.setReverseComplement(false);
  ASSERT_FALSE(seedChain.isReverseComplement());
  seedChain.setReverseComplement(true);
  ASSERT_TRUE(seedChain.isReverseComplement());
}

TEST(SeedChain, AddFirstSeed)
{
  SeedChain seedChain;
  // check that the seed chain is empty
  ASSERT_EQ(seedChain.begin(), seedChain.end());
  Read           read;
  const uint64_t id       = 0;
  const unsigned position = 0;
  read.init(Name(), Bases(151), Qualities(151), id, position);
  const unsigned primaryLength = 17;
  const unsigned readPosition  = 42;
  const Seed     seed(&read, readPosition, primaryLength);
  // make it a reverse complement and add a single non-random-sample hit
  seedChain.setReverseComplement(true);
  unsigned halfExtension = 0;
  seedChain.addSeedPosition({seed, 100000UL, halfExtension}, false);
  ASSERT_FALSE(seedChain.hasOnlyRandomSamples());
  ASSERT_TRUE(seedChain.isReverseComplement());
  seedChain.clear();
  ASSERT_FALSE(seedChain.isReverseComplement());
  ASSERT_TRUE(seedChain.hasOnlyRandomSamples());
  // make it a forward sequence and add a single random-sample
  seedChain.setReverseComplement(false);
  seedChain.addSeedPosition({seed, 100000UL, halfExtension}, true);
  ASSERT_FALSE(seedChain.isReverseComplement());
  ASSERT_TRUE(seedChain.hasOnlyRandomSamples());
}

TEST(SeedChain, AddSecondSeed)
{
  //FAIL() << "NOT IMPLEMENTED";
}

TEST(SeedChain, Clear)
{
  //FAIL() << "NOT IMPLEMENTED";
}

#if 0

TEST(SeedChainTest, TestConstructor)
{
  const uint64_t seq_pos     = 100000;
  const uint64_t ref_pos     = 35000;
  bool           rc          = false;
  bool           rs          = false;
  const uint8_t  ext_size    = 0;
  const int32_t  seed_length = 27;
  const uint32_t read_id     = 0x12345678;

  dragenos::align::SeedChain sc(seq_pos, ref_pos, seed_length, rc, rs, ext_size, read_id);
  ASSERT_EQ(1u, sc.size());
  ASSERT_EQ(ref_pos, sc.get_first_ref_pos());
  ASSERT_EQ(seq_pos, sc.get_first_seed_pos());
  ASSERT_EQ(ref_pos, sc.get_last_ref_pos());
  ASSERT_EQ(seq_pos, sc.get_last_seed_pos());
  ASSERT_EQ(read_id, sc.get_read_id());
}

TEST(SeedChainTest, RejectMismatchedRC)
{
  const uint64_t seq_pos     = 1;
  const uint64_t ref_pos     = 100001;
  const int32_t  seed_length = 27;
  bool           rc          = false;
  bool           rs          = false;
  const uint8_t  ext_size    = 0;
  const uint32_t read_id     = 0x12345678;

  dragenos::align::SeedChain sc(seq_pos, ref_pos, seed_length, rc, rs, ext_size, read_id);
  ASSERT_EQ(1u, sc.size());
  ASSERT_EQ(ref_pos, sc.get_first_ref_pos());
  ASSERT_EQ(seq_pos, sc.get_first_seed_pos());
  ASSERT_EQ(ref_pos, sc.get_last_ref_pos());
  ASSERT_EQ(seq_pos, sc.get_last_seed_pos());
  ASSERT_EQ(read_id, sc.get_read_id());

  ASSERT_FALSE(sc.accept(ref_pos, seq_pos + 2, !rc, rs, ext_size));
  ASSERT_EQ(1u, sc.size());
  ASSERT_EQ(ref_pos, sc.get_first_ref_pos());
  ASSERT_EQ(seq_pos, sc.get_first_seed_pos());
  ASSERT_EQ(ref_pos, sc.get_last_ref_pos());
  ASSERT_EQ(seq_pos, sc.get_last_seed_pos());
  ASSERT_EQ(read_id, sc.get_read_id());
}

TEST(SeedChainTest, SingleDeletion)
{
  uint64_t       seq_pos     = 1;
  uint64_t       ref_pos     = 100001;
  const int32_t  seed_length = 27;
  bool           rc          = false;
  bool           rs          = false;
  const uint8_t  ext_size    = 0;
  const uint32_t read_id     = 0x12345678;

  dragenos::align::SeedChain sc(seq_pos, ref_pos, seed_length, rc, rs, ext_size, read_id);
  ASSERT_EQ(1u, sc.size());
  ASSERT_EQ(ref_pos, sc.get_first_ref_pos());
  ASSERT_EQ(seq_pos, sc.get_first_seed_pos());
  ASSERT_EQ(ref_pos, sc.get_last_ref_pos());
  ASSERT_EQ(seq_pos, sc.get_last_seed_pos());
  ASSERT_EQ(read_id, sc.get_read_id());

  ref_pos += 3;
  seq_pos += 2;
  ASSERT_TRUE(sc.accept(seq_pos, ref_pos, rc, rs, ext_size));
  ASSERT_EQ(2u, sc.size());
  ASSERT_EQ(ref_pos, sc.get_last_ref_pos());
  ASSERT_EQ(seq_pos, sc.get_last_seed_pos());
  ASSERT_EQ(read_id, sc.get_read_id());
}

TEST(SeedChainTest, SingleInsertion)
{
  uint64_t       seq_pos     = 1;
  uint64_t       ref_pos     = 1000001;
  const int32_t  seed_length = 27;
  bool           rc          = false;
  bool           rs          = false;
  const uint8_t  ext_size    = 0;
  const uint32_t read_id     = 0x12345678;

  dragenos::align::SeedChain sc(seq_pos, ref_pos, seed_length, rc, rs, ext_size, read_id);
  ASSERT_EQ(1u, sc.size());

  ref_pos += 2;
  seq_pos += 3;
  ASSERT_TRUE(sc.accept(seq_pos, ref_pos, rc, rs, ext_size));
  ASSERT_EQ(2u, sc.size());
  ASSERT_EQ(ref_pos, sc.get_last_ref_pos());
  ASSERT_EQ(seq_pos, sc.get_last_seed_pos());
  ASSERT_EQ(read_id, sc.get_read_id());
}

TEST(SeedChainTest, LongStringOfInsertions)
{
  uint64_t       seq_pos        = 1;
  uint64_t       ref_pos        = 1000001;
  const int32_t  seed_length    = 27;
  const int      NUM_INSERTIONS = 40;
  bool           rc             = false;
  bool           rs             = false;
  const uint8_t  ext_size       = 0;
  const uint32_t read_id        = 0x12345678;

  dragenos::align::SeedChain sc(seq_pos, ref_pos, seed_length, rc, rs, ext_size, read_id);
  ASSERT_EQ(1u, sc.size());

  for (int i = 1; i < NUM_INSERTIONS; ++i) {
    seq_pos += 3;
    ref_pos += 2;

    ASSERT_TRUE(sc.accept(seq_pos, ref_pos, rc, rs, ext_size));
    ASSERT_EQ(static_cast<uint32_t>(i) + 1, sc.size());
    ASSERT_EQ(ref_pos, sc.get_last_ref_pos());
    ASSERT_EQ(seq_pos, sc.get_last_seed_pos());
    ASSERT_EQ(read_id, sc.get_read_id());
  }
}

TEST(SeedChainTest, RejectRCMismatch)
{
  const uint64_t seq_pos     = 1;
  const uint64_t ref_pos     = 1000001;
  const int32_t  seed_length = 27;
  bool           rc          = false;
  bool           rs          = false;
  const uint8_t  ext_size    = 0;
  const uint32_t read_id     = 0x12345678;

  dragenos::align::SeedChain sc(seq_pos, ref_pos, seed_length, rc, rs, ext_size, read_id);
  ASSERT_EQ(1u, sc.size());

  ASSERT_FALSE(sc.accept(seq_pos + 1, ref_pos + 1, !rc, rs, ext_size));
  ASSERT_EQ(1u, sc.size());
  ASSERT_EQ(ref_pos, sc.get_first_ref_pos());
  ASSERT_EQ(seq_pos, sc.get_first_seed_pos());
  ASSERT_EQ(ref_pos, sc.get_last_ref_pos());
  ASSERT_EQ(seq_pos, sc.get_last_seed_pos());
  ASSERT_EQ(read_id, sc.get_read_id());
}

// diameter test, none of the seeds being "too old"
TEST(SeedChainTest, TestDiameterLimitTooEarlyInReference)
{
  uint64_t       first_seq_pos = 1;
  uint64_t       first_ref_pos = 1000001;
  uint64_t       seq_pos       = first_seq_pos;
  uint64_t       ref_pos       = first_ref_pos;
  const int32_t  seed_length   = 27;
  bool           rc            = false;
  bool           rs            = false;
  const uint8_t  ext_size      = 0;
  const uint32_t read_id       = 0x12345678;

  dragenos::align::SeedChain sc(seq_pos, ref_pos, seed_length, rc, rs, ext_size, read_id);
  ASSERT_EQ(1u, sc.size());
  ASSERT_EQ(ref_pos, sc.get_first_ref_pos());
  ASSERT_EQ(seq_pos, sc.get_first_seed_pos());
  ASSERT_EQ(ref_pos, sc.get_last_ref_pos());
  ASSERT_EQ(seq_pos, sc.get_last_seed_pos());
  ASSERT_EQ(read_id, sc.get_read_id());

  seq_pos += 248;
  ref_pos += 200;
  ASSERT_TRUE(sc.accept(seq_pos, ref_pos, rc, rs, ext_size));
  ASSERT_EQ(2u, sc.size());
  ASSERT_EQ(ref_pos, sc.get_last_ref_pos());
  ASSERT_EQ(seq_pos, sc.get_last_seed_pos());
  ASSERT_EQ(read_id, sc.get_read_id());

  // Should be rejected due to the diagonal being too early
  // in the reference relative to previously accepted seeds
  ASSERT_FALSE(sc.accept(seq_pos + 49, ref_pos, rc, rs, ext_size));
  ASSERT_EQ(2u, sc.size());
  ASSERT_EQ(first_ref_pos, sc.get_first_ref_pos());
  ASSERT_EQ(first_seq_pos, sc.get_first_seed_pos());
  ASSERT_EQ(ref_pos, sc.get_last_ref_pos());
  ASSERT_EQ(seq_pos, sc.get_last_seed_pos());
  ASSERT_EQ(read_id, sc.get_read_id());
}

// diameter test, none of the seeds being "too old"
TEST(SeedChainTest, TestDiameterLimitTooLateInReference)
{
  const uint64_t seq_pos     = 10;
  const uint64_t ref_pos     = 3500007;
  const int32_t  seed_length = 27;
  bool           rc          = false;
  bool           rs          = false;
  const uint8_t  ext_size    = 0;
  const uint32_t read_id     = 0x12345678;

  dragenos::align::SeedChain sc(seq_pos, ref_pos, seed_length, rc, rs, ext_size, read_id);
  ASSERT_EQ(1u, sc.size());

  ASSERT_TRUE(sc.accept(seq_pos + 200, ref_pos + 250, rc, rs, ext_size));
  ASSERT_EQ(2u, sc.size());

  // Should be rejected due to the diagonal being too late
  // in the reference relative to previously accepted seeds
  ASSERT_FALSE(sc.accept(seq_pos + 200, ref_pos + 252, rc, rs, ext_size));
  ASSERT_EQ(2u, sc.size());
  ASSERT_EQ(seq_pos, sc.get_first_seed_pos());
  ASSERT_EQ(ref_pos, sc.get_first_ref_pos());
  ASSERT_EQ(seq_pos + 200, sc.get_last_seed_pos());
  ASSERT_EQ(ref_pos + 250, sc.get_last_ref_pos());
  ASSERT_EQ(read_id, sc.get_read_id());
}

// radius test.
TEST(SeedChainTest, TestRadiusLimit)
{
  const uint64_t seq_pos     = 100;
  const uint64_t ref_pos     = 1000000;
  const int32_t  seed_length = 27;
  bool           rc          = false;
  bool           rs          = false;
  const uint8_t  ext_size    = 0;
  const uint32_t read_id     = 0x12345678;

  dragenos::align::SeedChain sc(seq_pos, ref_pos, seed_length, rc, rs, ext_size, read_id);
  ASSERT_EQ(1u, sc.size());

  ASSERT_FALSE(sc.accept(seq_pos, ref_pos + 32, rc, rs, ext_size));  // fails the radius limit
  ASSERT_FALSE(sc.accept(seq_pos, ref_pos + 36, rc, rs, ext_size));  // fails the radius limit
  ASSERT_FALSE(sc.accept(seq_pos, ref_pos - 32, rc, rs, ext_size));  // fails the radius limit

  ASSERT_TRUE(sc.accept(seq_pos, ref_pos + 4, rc, rs, ext_size));  // passes the radius limit
  ASSERT_EQ(seq_pos, sc.get_last_seed_pos());
  ASSERT_EQ(ref_pos, sc.get_last_ref_pos());

  // ref_pos+36 now this passes the radius limit and the diameter limit!
  ASSERT_TRUE(sc.accept(seq_pos, ref_pos + 35, rc, rs, ext_size));
  ASSERT_EQ(seq_pos, sc.get_last_seed_pos());
  ASSERT_EQ(ref_pos, sc.get_last_ref_pos());
  ASSERT_EQ(seq_pos, sc.get_first_seed_pos());
  ASSERT_EQ(ref_pos, sc.get_first_ref_pos());
  ASSERT_EQ(read_id, sc.get_read_id());
}

TEST(SeedChainTest, TestTerminationOnEverythingAncient)
{
  const uint64_t seq_pos     = 10;
  const uint64_t ref_pos     = 3500007;
  const int32_t  seed_length = 27;
  bool           rc          = false;
  bool           rs          = false;
  const uint8_t  ext_size    = 0;
  const uint32_t read_id     = 0x12345678;

  dragenos::align::SeedChain sc(seq_pos, ref_pos, seed_length, rc, rs, ext_size, read_id);
  ASSERT_EQ(1u, sc.size());

  ASSERT_TRUE(sc.accept(seq_pos + 200, ref_pos + 250, rc, rs, ext_size));

  // If we skip way down the read without any more matches,
  // the seed chain should be terminated and reject all further matches.
  ASSERT_FALSE(sc.accept(seq_pos + 701, ref_pos + 951, rc, rs, ext_size));

  ASSERT_EQ(2u, sc.size());
  ASSERT_EQ(seq_pos, sc.get_first_seed_pos());
  ASSERT_EQ(ref_pos, sc.get_first_ref_pos());
  ASSERT_EQ(seq_pos + 200, sc.get_last_seed_pos());
  ASSERT_EQ(ref_pos + 250, sc.get_last_ref_pos());
  ASSERT_EQ(read_id, sc.get_read_id());
}

TEST(SeedChainTest, TestAcceptSlowDrift)
{
  uint64_t       seq_pos     = 10;
  uint64_t       ref_pos     = 1000005;
  int32_t        seed_length = 27;
  const int      NUM_SEEDS   = 2000;
  bool           rc          = false;
  bool           rs          = false;
  const uint8_t  ext_size    = 0;
  const uint32_t read_id     = 0x12345678;

  dragenos::align::SeedChain sc(seq_pos, ref_pos, seed_length, rc, rs, ext_size, read_id);
  ASSERT_EQ(1u, sc.size());

  for (int i = 1; i < NUM_SEEDS; ++i) {
    seq_pos += 28;
    ref_pos += 29;

    ASSERT_TRUE(sc.accept(seq_pos, ref_pos, rc, rs, ext_size));
    ASSERT_EQ(static_cast<uint32_t>(i) + 1, sc.size());
    ASSERT_EQ(ref_pos, sc.get_last_ref_pos());
    ASSERT_EQ(seq_pos, sc.get_last_seed_pos());
    ASSERT_EQ(read_id, sc.get_read_id());
  }
}

#endif
