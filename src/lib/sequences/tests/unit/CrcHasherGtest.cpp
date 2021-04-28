#include "CrcHasherMocks.hpp"
#include "gtest/gtest.h"

#include <string>
#include "sequences/CrcHasher.hpp"

TEST(CrcHasher, KnownHashValues)
{
  using dragenos::sequences::CrcHasher;
  using dragenos::sequences::CrcPolynomial;
  // polynomial and value taken and verified from Hashtable V8
  const CrcPolynomial primaryPolynomial(54, std::string("2C991CE6A8DD55"));
  const CrcHasher     primaryHasher(primaryPolynomial);
  ASSERT_EQ(primaryHasher.getHash64(0x3543543543), 0x35de20855c4d8e);
  ASSERT_EQ(primaryHasher.getHash64(0x10d50d50d50), 0x2643a64efec778);
  ASSERT_EQ(primaryHasher.getHash64(0x14354354354), 0x24e8498f420b6f);
  ASSERT_EQ(primaryHasher.getHash64(0x150d50d50d5), 0x9287b35d36195);
  ASSERT_EQ(primaryHasher.getHash64(0xd50d50d50d), 0x22d3a73e8851c7);
}
