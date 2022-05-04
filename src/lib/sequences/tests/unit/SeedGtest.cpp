#include "SeedMocks.hpp"

#include "gtest/gtest.h"

#include <sstream>
#include <string>

#include "common/Exceptions.hpp"
#include "sequences/Seed.hpp"

namespace dragenos {
namespace common {

InvalidParameterException::InvalidParameterException(const std::string& message)
  : std::logic_error(message), ExceptionData(17, message)
{
}

PreConditionException::PreConditionException(const std::string& message)
  : std::logic_error(message), ExceptionData(17, message)
{
}

ExceptionData::ExceptionData(int e, std::string const&) : errorNumber_(e) {}

}  // namespace common
}  // namespace dragenos

TEST(SeedTest, getSeedOffsets)
{
  {
    const size_t readLength = 9;
    const size_t seedLength = 1;
    // every 4th seed, forcing the last 3
    // 0, 4, 6, 7, 8
    const uint32_t period     = 4;
    const uint32_t pattern    = 0x01;
    const uint8_t  forceLastN = 3;
    const auto     seedOffsets =
        dragenos::sequences::Seed::getSeedOffsets(readLength, seedLength, period, pattern, forceLastN);
    ASSERT_EQ(5u, seedOffsets.size());
    ASSERT_EQ(0u, seedOffsets[0]);
    ASSERT_EQ(4u, seedOffsets[1]);
    ASSERT_EQ(6u, seedOffsets[2]);
    ASSERT_EQ(7u, seedOffsets[3]);
    ASSERT_EQ(8u, seedOffsets[4]);
  }
  {
    const size_t readLength = 151;
    const size_t seedLength = 17;
    // every seed at even offset, forcing the last 3
    const uint32_t period     = 2;
    const uint32_t pattern    = 0x01;
    const uint8_t  forceLastN = 3;
    const auto     seedOffsets =
        dragenos::sequences::Seed::getSeedOffsets(readLength, seedLength, period, pattern, forceLastN);
    ASSERT_EQ(69u, seedOffsets.size());
    ASSERT_EQ(0u, seedOffsets[0]);
    ASSERT_EQ(2u, seedOffsets[1]);
    ASSERT_EQ(4u, seedOffsets[2]);
    ASSERT_EQ(130u, seedOffsets[65]);
    ASSERT_EQ(132u, seedOffsets[66]);
    ASSERT_EQ(133u, seedOffsets[67]);
    ASSERT_EQ(134u, seedOffsets[68]);
  }
}

TEST(SeedTest, generateReverseComplement)
{
  using dragenos::sequences::Seed;
  ASSERT_EQ(0x0u, Seed::generateReverseComplement(0x0UL, 0));
  ASSERT_EQ(0x0u, Seed::generateReverseComplement(0x3UL, 1));
  ASSERT_EQ(0x1u, Seed::generateReverseComplement(0x2UL, 1));
  ASSERT_EQ(0x2u, Seed::generateReverseComplement(0x1UL, 1));
  ASSERT_EQ(0x3u, Seed::generateReverseComplement(0x0UL, 1));
  // 2 bases that are reverse complement of each other
  ASSERT_EQ(0x3u, Seed::generateReverseComplement(0x3UL, 2));  // TA
  ASSERT_EQ(0x6u, Seed::generateReverseComplement(0x6UL, 2));  // GC
  ASSERT_EQ(0x9u, Seed::generateReverseComplement(0x9UL, 2));  // CG
  ASSERT_EQ(0xCu, Seed::generateReverseComplement(0xCUL, 2));  // AT
  // 2 bases that are the same
  ASSERT_EQ(0x0u, Seed::generateReverseComplement(0xFUL, 2));  // AA
  ASSERT_EQ(0x5u, Seed::generateReverseComplement(0xAUL, 2));  // CC
  ASSERT_EQ(0xAu, Seed::generateReverseComplement(0x5UL, 2));  // GG
  ASSERT_EQ(0xFu, Seed::generateReverseComplement(0x0UL, 2));  // TT
  // other 2 base combinations
  ASSERT_EQ(0x1u, Seed::generateReverseComplement(0xBUL, 2));  // CA
  ASSERT_EQ(0x2u, Seed::generateReverseComplement(0x7UL, 2));  // GA
  ASSERT_EQ(0xDu, Seed::generateReverseComplement(0x8UL, 2));  // CT
  ASSERT_EQ(0xEu, Seed::generateReverseComplement(0x4UL, 2));  // GT
  ASSERT_EQ(0x4u, Seed::generateReverseComplement(0xEUL, 2));  // AC
  ASSERT_EQ(0x7u, Seed::generateReverseComplement(0x2UL, 2));  // TC
  ASSERT_EQ(0x8u, Seed::generateReverseComplement(0xDUL, 2));  // AG
  ASSERT_EQ(0xBu, Seed::generateReverseComplement(0x1UL, 2));  // TG
  // some random sequences
  ASSERT_EQ(0x5F21u, Seed::generateReverseComplement(0xB70AUL, 8));       // CAGATTCC -> GGAATCTG
  ASSERT_EQ(0x35F21u, Seed::generateReverseComplement(0x2DC28UL, 9));     // CAGATTCCT -> AGGAATCTG
  ASSERT_EQ(0xE5F21u, Seed::generateReverseComplement(0xB70A4UL, 10));    // CAGATTCCGT -> ACGGAATCTG
  ASSERT_EQ(0x275F21u, Seed::generateReverseComplement(0x2DC289UL, 11));  // CAGATTCCTCG -> CGAGGAATCTG
  ASSERT_EQ(0xDE5F21u, Seed::generateReverseComplement(0xB70A48UL, 12));  // CAGATTCCGTCT -> AGACGGAATCTG
}

TEST(SeedTest, Constructor)
{
  using dragenos::sequences::Read;
  using dragenos::sequences::Seed;
  Read           read;
  const unsigned primaryLength = 21;
  const unsigned readPosition  = 12;
  Seed           seed(&read, readPosition, primaryLength);
  ASSERT_EQ(&read, seed.getRead());
  ASSERT_EQ(readPosition, seed.getReadPosition());
  ASSERT_EQ(primaryLength, seed.getPrimaryLength());
  ASSERT_NO_THROW(Seed(&read, 0, 32));
  ASSERT_THROW(Seed(&read, 0, 33), dragenos::common::InvalidParameterException);
}

dragenos::sequences::Seed::Data getBases(
    const unsigned begin, const unsigned end, const dragenos::sequences::Read& read)
{
  dragenos::sequences::Seed::Data data = 0;
  for (unsigned i = end; begin != i; --i) {
    data <<= 2;
    data |= (read.getBase2bpb(i - 1) & 3);
  }
  return data;
}

TEST(SeedTest, getPrimaryData)
{
  using dragenos::sequences::Read;
  using dragenos::sequences::Seed;
  Read read;
  // single base seeds
  ASSERT_EQ(151u, read.getLength());
  for (unsigned i = 0; read.getLength() > i; ++i) {
    ASSERT_EQ((unsigned long)read.getBase2bpb(i), Seed(&read, i, 1).getPrimaryData(false)) << "i: " << i;
    ASSERT_EQ((unsigned long)(~read.getBase2bpb(i)) & 3, Seed(&read, i, 1).getPrimaryData(true)) << "i: " << i;
  }
  ASSERT_THROW(
      Seed(&read, read.getLength(), 1).getPrimaryData(false), dragenos::common::PreConditionException);
  ASSERT_NO_THROW(Seed(&read, read.getLength() - 32, 32).getPrimaryData(false));
  ASSERT_EQ(
      getBases(read.getLength() - 32, read.getLength(), read),
      Seed(&read, read.getLength() - 32, 32).getPrimaryData(false));
  ASSERT_THROW(
      Seed(&read, read.getLength() - 32 + 1, 32).getPrimaryData(false),
      dragenos::common::PreConditionException);
  //typedef Seed::Data Data;
  //const Data bases0To20 = getBases(0, 21, read);
  // Check forward and reverse seed of length 21 bases at read positions 0, 10 and 130
  ASSERT_EQ(getBases(0, 21, read), Seed(&read, 0, 21).getPrimaryData(false));
  ASSERT_EQ(
      Seed::generateReverseComplement(getBases(0, 21, read), 21), Seed(&read, 0, 21).getPrimaryData(true));
  ASSERT_EQ(getBases(10, 31, read), Seed(&read, 10, 21).getPrimaryData(false));
  ASSERT_EQ(
      Seed::generateReverseComplement(getBases(10, 31, read), 21), Seed(&read, 10, 21).getPrimaryData(true));
  ASSERT_EQ(
      getBases(read.getLength() - 21, read.getLength(), read),
      Seed(&read, read.getLength() - 21, 21).getPrimaryData(false));
  ASSERT_EQ(
      Seed::generateReverseComplement(getBases(read.getLength() - 21, read.getLength(), read), 21),
      Seed(&read, read.getLength() - 21, 21).getPrimaryData(true));
}

TEST(SeedTest, getExtendedData)
{
  using dragenos::common::InvalidParameterException;
  using dragenos::common::PreConditionException;
  using dragenos::sequences::Read;
  using dragenos::sequences::Seed;
  typedef Seed::Data Data;
  Read               read;
  // single base wings
  ASSERT_EQ(151u, read.getLength());
  for (unsigned i = 1; read.getLength() > i + 1; ++i) {
    const Data l = read.getBase2bpb(i - 1);
    const Data r = read.getBase2bpb(i + 1);
    ASSERT_EQ(l | (r << 2), Seed(&read, i, 1).getExtendedData(0, 1, false)) << "i: " << i;
    ASSERT_EQ(((~r) & 3) | (((~l) & 3) << 2), Seed(&read, i, 1).getExtendedData(0, 1, true)) << "i: " << i;
  }
  // Failed precondition at the beginning
  ASSERT_NO_THROW(Seed(&read, 10, 1).getExtendedData(0, 10, false));
  ASSERT_THROW(Seed(&read, 10, 1).getExtendedData(0, 11, false), PreConditionException);
  ASSERT_THROW(Seed(&read, 10, 1).getExtendedData(10, 11, false), PreConditionException);
  // Failed precondition at the end
  ASSERT_NO_THROW(Seed(&read, read.getLength() - 11, 1).getExtendedData(0, 10, false));
  ASSERT_THROW(Seed(&read, read.getLength() - 11, 1).getExtendedData(0, 11, false), PreConditionException);
  ASSERT_THROW(Seed(&read, read.getLength() - 11, 1).getExtendedData(10, 11, false), PreConditionException);
  // from > to
  ASSERT_THROW(Seed(&read, 35, 1).getExtendedData(20, 11, false), InvalidParameterException);
  // wings too long
  ASSERT_NO_THROW(Seed(&read, 35, 1).getExtendedData(0, 16, false));
  ASSERT_THROW(Seed(&read, 35, 1).getExtendedData(0, 17, false), InvalidParameterException);
  ASSERT_NO_THROW(Seed(&read, 35, 1).getExtendedData(16, 32, false));
  ASSERT_THROW(Seed(&read, 35, 1).getExtendedData(16, 33, false), InvalidParameterException);
  // 2x6 bases wings
  const Seed seed(&read, 20, 21);
  {
    const Data l = getBases(12, 18, read);
    const Data r = getBases(43, 49, read);
    ASSERT_EQ(l | (r << 12), seed.getExtendedData(2, 8, false));
    // must use Reverse complement on each wing
    ASSERT_EQ(
        Seed::generateReverseComplement(r, 6) | (Seed::generateReverseComplement(l, 6) << 12),
        seed.getExtendedData(2, 8, true));
  }
}
