#include "gtest/gtest.h"

#include <algorithm>
#include <type_traits>

#include "align/Cigar.hpp"

TEST(Cigar, OperationNames)
{
  using dragenos::align::Cigar;
  ASSERT_EQ(Cigar::SEQUENCE_MISMATCH + 1, std::extent<decltype(Cigar::OPERATION_NAMES)>::value);
  ASSERT_EQ('M', Cigar::getOperationName(Cigar::ALIGNMENT_MATCH));
  ASSERT_EQ('I', Cigar::getOperationName(Cigar::INSERT));
  ASSERT_EQ('D', Cigar::getOperationName(Cigar::DELETE));
  ASSERT_EQ('N', Cigar::getOperationName(Cigar::SKIP));
  ASSERT_EQ('S', Cigar::getOperationName(Cigar::SOFT_CLIP));
  ASSERT_EQ('H', Cigar::getOperationName(Cigar::HARD_CLIP));
  ASSERT_EQ('P', Cigar::getOperationName(Cigar::PAD));
  ASSERT_EQ('=', Cigar::getOperationName(Cigar::SEQUENCE_MATCH));
  ASSERT_EQ('X', Cigar::getOperationName(Cigar::SEQUENCE_MISMATCH));
}

TEST(Cigar, Operations)
{
  using dragenos::align::Cigar;
  Cigar cigar;
  ASSERT_EQ(0, cigar.getNumberOfOperations());
  const Cigar::Operation opM17(Cigar::ALIGNMENT_MATCH, 17);
  const Cigar::Operation opI21(Cigar::INSERT, 21);
  const Cigar::Operation opD35(Cigar::DELETE, 35);
  cigar.emplace_back(opM17.first, opM17.second);
  ASSERT_EQ(1, cigar.getNumberOfOperations());
  cigar.push_back(opI21);
  ASSERT_EQ(2, cigar.getNumberOfOperations());
  cigar.push_back(Cigar::Operation(opD35.first, opD35.second));
  ASSERT_EQ(3, cigar.getNumberOfOperations());
  ASSERT_EQ(Cigar::ALIGNMENT_MATCH, cigar.getOperations()->first);
  ASSERT_EQ(Cigar::INSERT, (cigar.getOperations() + 1)->first);
  ASSERT_EQ(Cigar::DELETE, (cigar.getOperations() + 2)->first);
  ASSERT_EQ(17, cigar.getOperations()->second);
  ASSERT_EQ(21, (cigar.getOperations() + 1)->second);
  ASSERT_EQ(35, (cigar.getOperations() + 2)->second);
  cigar.clear();
  ASSERT_EQ(0, cigar.getNumberOfOperations());
}
