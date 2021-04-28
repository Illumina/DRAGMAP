#include "gtest/gtest.h"

#include "SeedPositionMocks.hpp"

#include "map/SeedPosition.hpp"

using dragenos::map::SeedPosition;
using dragenos::sequences::Read;
using dragenos::sequences::Seed;
typedef SeedPosition::ReferencePosition ReferencePosition;

TEST(SeedPosition, getReferencePosition)
{
  const Read              read;
  const unsigned          primaryLength = 17;
  const unsigned          readPosition  = 42;
  const Seed              seed(&read, readPosition, primaryLength);
  const ReferencePosition referencePosition = 1000 + readPosition;
  // TODO: add tests for reverse seeds and extended seeds
  unsigned           halfExtension = 0;
  const SeedPosition seedPosition(seed, referencePosition, halfExtension);
  ASSERT_EQ(seedPosition.getReferencePosition(), referencePosition);
}

TEST(SeedPosition, getFirstProjection)
{
  const Read     read;
  const auto     readLength    = read.getLength();
  const unsigned primaryLength = 17;
  const unsigned readPosition  = 42;
  ASSERT_GT(readLength, readPosition);
  const Seed              seed(&read, readPosition, primaryLength);
  const ReferencePosition referencePosition = 1000 + readPosition;
  // TODO: add tests for reverse seeds and extended seeds
  unsigned           halfExtension = 0;
  const SeedPosition seedPosition(seed, referencePosition, halfExtension);
  ASSERT_EQ(seedPosition.getFirstProjection(false), referencePosition - readPosition);
  ASSERT_EQ(seedPosition.getFirstProjection(true), referencePosition - (readLength - readPosition - 1));
}

TEST(SeedPosition, getLastProjection)
{
  const Read     read;
  const auto     readLength    = read.getLength();
  const unsigned primaryLength = 17;
  const unsigned readPosition  = 42;
  ASSERT_GT(readLength, readPosition);
  const Seed              seed(&read, readPosition, primaryLength);
  const ReferencePosition referencePosition = 1000 + readPosition;
  // TODO: add tests for reverse seeds and extended seeds
  unsigned           halfExtension = 0;
  const SeedPosition seedPosition(seed, referencePosition, halfExtension);
  ASSERT_EQ(seedPosition.getLastProjection(false), seedPosition.getFirstProjection(false) + readLength - 1);
  ASSERT_EQ(seedPosition.getLastProjection(true), seedPosition.getFirstProjection(true) + readLength - 1);
}
