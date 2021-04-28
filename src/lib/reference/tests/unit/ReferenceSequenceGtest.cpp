#include "gtest/gtest.h"

#include <algorithm>
#include <cassert>

#include "reference/ReferenceSequence.hpp"

static std::string bases = "ACGTAACCGGTTAAACCCGGGTTT";
// the encoding is little endian. Bor byte at position B, 4 LSB are base (2xB) and 4 MSB are base (2xB + 1)
static std::vector<uint8_t> encoded{0x21, 0x84, 0x11, 0x22, 0x44, 0x88, 0x11, 0x21, 0x22, 0x44, 0x84, 0x88};
unsigned char               encodeBase(const char base);
char                        decodeBase(const unsigned char encodedBase);
std::vector<unsigned char>  generateSequence();

TEST(ReferenceSequence, decodeBase)
{
  using dragenos::reference::ReferenceSequence;
  ASSERT_EQ('P', ReferenceSequence::decodeBase(0x0));
  ASSERT_EQ('A', ReferenceSequence::decodeBase(0x1));
  ASSERT_EQ('C', ReferenceSequence::decodeBase(0x2));
  ASSERT_EQ('M', ReferenceSequence::decodeBase(0x3));
  ASSERT_EQ('G', ReferenceSequence::decodeBase(0x4));
  ASSERT_EQ('R', ReferenceSequence::decodeBase(0x5));
  ASSERT_EQ('S', ReferenceSequence::decodeBase(0x6));
  ASSERT_EQ('V', ReferenceSequence::decodeBase(0x7));
  ASSERT_EQ('T', ReferenceSequence::decodeBase(0x8));
  ASSERT_EQ('W', ReferenceSequence::decodeBase(0x9));
  ASSERT_EQ('Y', ReferenceSequence::decodeBase(0xA));
  ASSERT_EQ('H', ReferenceSequence::decodeBase(0xB));
  ASSERT_EQ('K', ReferenceSequence::decodeBase(0xC));
  ASSERT_EQ('D', ReferenceSequence::decodeBase(0xD));
  ASSERT_EQ('B', ReferenceSequence::decodeBase(0xE));
  ASSERT_EQ('N', ReferenceSequence::decodeBase(0xF));
  for (unsigned i = 0x10; 0x100 > i; ++i) {
    ASSERT_EQ(ReferenceSequence::decodeBase(i), ReferenceSequence::decodeBase(i & 0xF)) << "i: " << i;
  }
}

TEST(ReferenceSequence, translateTo2bpb)
{
  using dragenos::reference::ReferenceSequence;
  ASSERT_EQ(0, ReferenceSequence::translateTo2bpb(0x0));
  ASSERT_EQ(0, ReferenceSequence::translateTo2bpb(0x1));
  ASSERT_EQ(1, ReferenceSequence::translateTo2bpb(0x2));
  ASSERT_EQ(0, ReferenceSequence::translateTo2bpb(0x3));
  ASSERT_EQ(2, ReferenceSequence::translateTo2bpb(0x4));
  ASSERT_EQ(0, ReferenceSequence::translateTo2bpb(0x5));
  ASSERT_EQ(0, ReferenceSequence::translateTo2bpb(0x6));
  ASSERT_EQ(0, ReferenceSequence::translateTo2bpb(0x7));
  ASSERT_EQ(3, ReferenceSequence::translateTo2bpb(0x8));
  ASSERT_EQ(0, ReferenceSequence::translateTo2bpb(0x9));
  ASSERT_EQ(0, ReferenceSequence::translateTo2bpb(0xA));
  ASSERT_EQ(0, ReferenceSequence::translateTo2bpb(0xB));
  ASSERT_EQ(0, ReferenceSequence::translateTo2bpb(0xC));
  ASSERT_EQ(0, ReferenceSequence::translateTo2bpb(0xD));
  ASSERT_EQ(0, ReferenceSequence::translateTo2bpb(0xE));
  ASSERT_EQ(0, ReferenceSequence::translateTo2bpb(0xF));
  for (unsigned i = 0x10; 0x100 > i; ++i) {
    ASSERT_EQ(ReferenceSequence::translateTo2bpb(i), ReferenceSequence::translateTo2bpb(i & 0xF))
        << "i: " << i;
  }
}

TEST(ReferenceSequence, generateSequence)
{
  const auto generated = generateSequence();
  ASSERT_EQ(encoded.size(), generated.size());
  for (size_t i = 0; encoded.size() > i; ++i) {
    ASSERT_EQ(encoded[i], generated[i]) << "i: " << i;
  }
}

TEST(ReferenceSequence, getBase)
{
  const auto sequence = generateSequence();
  ASSERT_EQ(sequence.size() * 2, bases.size());
  using dragenos::reference::ReferenceSequence;
  typedef ReferenceSequence::Region Region;
  ReferenceSequence referenceSequence(std::vector<Region>(), sequence.data(), sequence.size());
  ASSERT_EQ(1, referenceSequence.getBase(0));
  ASSERT_EQ(2, referenceSequence.getBase(1));
  ASSERT_EQ(4, referenceSequence.getBase(2));
  ASSERT_EQ(8, referenceSequence.getBase(3));
  ASSERT_EQ(8, referenceSequence.getBase(bases.size() - 1));
  ASSERT_THROW(referenceSequence.getBase(bases.size()), dragenos::common::InvalidParameterException);
}

TEST(ReferenceSequence, getSequenceCopy) {}

TEST(ReferenceSequence, getBaseReference) {}

std::vector<unsigned char> generateSequence()
{
  assert(0 == (bases.size() % 2));
  std::vector<unsigned char> sequence;
  bool                       even = true;
  for (auto c : bases) {
    if (even) {
      sequence.push_back(encodeBase(c));
    } else {
      sequence.back() |= (encodeBase(c) << 4);
    }
    even = !even;
  }
  return sequence;
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

char decodeBase(const unsigned char encoded)
{
  static const std::string bases = "NACNGNNNTN";
  return bases[std::min(bases.size(), size_t(encoded))];
}

// Mocks
//
namespace dragenos {
namespace common {
ExceptionData::ExceptionData(int e, std::string const&) : errorNumber_(e) {}
InvalidParameterException::InvalidParameterException(std::string const&) : std::logic_error("dummy") {}
}  // namespace common
}  // namespace dragenos
