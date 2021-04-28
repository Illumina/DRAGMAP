#include "gtest/gtest.h"

#include "sequences/CrcPolynomial.hpp"

typedef dragenos::sequences::CrcPolynomial crc_polynomial_t;

TEST(CrcPolynomial, DefaultConstructor)
{
  crc_polynomial_t poly(0, nullptr);
  for (int i = 0; i < 16; ++i) {
    EXPECT_EQ(poly.getData()[i], 0);
  }
}

TEST(CrcPolynomial, OneByte)
{
  crc_polynomial_t poly(8, "b1");
  EXPECT_EQ(0xb1, poly.getData()[0]);
  for (int i = 1; i < 16; ++i) EXPECT_EQ(0, poly.getData()[i]);
}

TEST(CrcPolynomial, OneBit)
{
  crc_polynomial_t poly(1, "1");
  EXPECT_EQ(poly.getData()[0], 1);
  for (int i = 2; i < 16; ++i) EXPECT_EQ(0, poly.getData()[i]);
}

TEST(CrcPolynomial, LessThan8Bytes)
{
  crc_polynomial_t poly(56, "2C991CE6A8DD55");
  EXPECT_EQ(poly.getData()[0], 0x55);
  EXPECT_EQ(poly.getData()[1], 0xDD);
  EXPECT_EQ(poly.getData()[2], 0xA8);
  EXPECT_EQ(poly.getData()[3], 0xE6);
  EXPECT_EQ(poly.getData()[4], 0x1C);
  EXPECT_EQ(poly.getData()[5], 0x99);
  EXPECT_EQ(poly.getData()[6], 0x2C);
  for (int i = 7; i < 16; ++i) EXPECT_EQ(poly.getData()[i], 0);
}

TEST(CrcPolynomial, MoreThan8Bytes)
{
  crc_polynomial_t poly(96, "ABCDEF01022C991CE6A8DD55");
  EXPECT_EQ(poly.getData()[0], 0x55);
  EXPECT_EQ(poly.getData()[1], 0xDD);
  EXPECT_EQ(poly.getData()[2], 0xA8);
  EXPECT_EQ(poly.getData()[3], 0xE6);
  EXPECT_EQ(poly.getData()[4], 0x1C);
  EXPECT_EQ(poly.getData()[5], 0x99);
  EXPECT_EQ(poly.getData()[6], 0x2C);
  EXPECT_EQ(poly.getData()[7], 0x02);
  EXPECT_EQ(poly.getData()[8], 0x01);
  EXPECT_EQ(poly.getData()[9], 0xEF);
  EXPECT_EQ(poly.getData()[10], 0xCD);
  EXPECT_EQ(poly.getData()[11], 0xAB);

  for (int i = 12; i < 16; ++i) EXPECT_EQ(poly.getData()[i], 0);
}

TEST(CrcPolynomial, EqualityOperator)
{
  crc_polynomial_t poly1(96, "ABCDEF01022C991CE6A8DD55");
  crc_polynomial_t poly2(96, "ABCDEF01022C991CE6A8DD55");
  EXPECT_EQ(poly1, poly2);
}

TEST(CrcPolynomial, StringEqualityOperator)
{
  crc_polynomial_t poly(96, "ABCDEF01022C991CE6A8DD55");
  EXPECT_EQ(poly, "ABCDEF01022C991CE6A8DD55");
}
