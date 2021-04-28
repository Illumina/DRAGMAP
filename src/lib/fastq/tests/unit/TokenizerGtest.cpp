#include "gtest/gtest.h"

#include <algorithm>
#include <string>
#include <type_traits>

#include "fastq/Tokenizer.hpp"

const std::string ONE_RECORD(
    "@NB551322:14:HFVLLBGX9:4:11401:24054:1050 2:N:0:CGGCTATG+CCGTCGCC\n"
    "NTGTCGGGGCAGGCAGGGCTCCTCGGGCAGCGGCTCATGAGAGAAGACGGAATCCTCCCCTGAGGAGCACGTAGAGCTCCGGGTGTCGGGAAAGCTGGGGG\n"
    "+\n"
    "#AAAAEEEEEEEEEEEEEEEEEEEEEEE6EEEEEEEEEEEEEEEEEAEAEEEEEEEEEEEEAEEE/EEEEAEEEAEEEEEEEEEEEEEEEEEEEEE/EEEE\n");

using dragenos::fastq::Tokenizer;

TEST(Tokenizer, SmallBuffer)
{
  std::stringstream stm(ONE_RECORD);
  Tokenizer         t(stm, 102);

  ASSERT_EQ(true, t.token().empty());
  ASSERT_THROW(t.next(), std::logic_error);
  ASSERT_EQ(false, stm.eof());
}

TEST(Tokenizer, EnoughBuffer)
{
  std::stringstream stm(ONE_RECORD);
  Tokenizer         t(stm, 1024);

  ASSERT_EQ(true, t.token().empty());
  ASSERT_EQ(true, t.next());
  ASSERT_EQ(true, stm.eof());
  ASSERT_EQ(false, t.token().empty());
  ASSERT_EQ(true, t.token().valid());
}

TEST(Tokenizer, JustEnoughBuffer)
{
  std::stringstream stm(ONE_RECORD);
  Tokenizer         t(stm, ONE_RECORD.size());

  ASSERT_EQ(true, t.token().empty());
  ASSERT_EQ(true, t.next());
  ASSERT_EQ(false, stm.eof());
  ASSERT_EQ(false, t.token().empty());
  ASSERT_EQ(true, t.token().valid());

  ASSERT_EQ(false, t.next());
  ASSERT_EQ(true, stm.eof());
  ASSERT_EQ(true, t.token().empty());
  ASSERT_EQ(false, t.token().valid());
}

TEST(Tokenizer, AlmostNotEnoughBuffer)
{
  std::stringstream stm(ONE_RECORD);
  // read whole record minus last newline
  Tokenizer t(stm, ONE_RECORD.size() - 1);

  ASSERT_EQ(true, t.token().empty());
  ASSERT_EQ(true, t.next());
  ASSERT_EQ(false, stm.eof());
  ASSERT_EQ(false, t.token().empty());
  ASSERT_EQ(true, t.token().valid());

  ASSERT_EQ(false, t.next());
  ASSERT_EQ(true, stm.eof());
  ASSERT_EQ(true, t.token().empty());
  ASSERT_EQ(false, t.token().valid());
}

const std::string TWO_RECORDS(
    "@NB551322:14:HFVLLBGX9:4:11401:24054:1050 2:N:0:CGGCTATG+CCGTCGCC\n"
    "NTGTCGGGGCAGGCAGGGCTCCTCGGGCAGCGGCTCATGAGAGAAGACGGAATCCTCCCCTGAGGAGCACGTAGAGCTCCGGGTGTCGGGAAAGCTGGGGG\n"
    "+\n"
    "#AAAAEEEEEEEEEEEEEEEEEEEEEEE6EEEEEEEEEEEEEEEEEAEAEEEEEEEEEEEEAEEE/EEEEAEEEAEEEEEEEEEEEEEEEEEEEEE/EEEE\n"
    "@NB551322:14:HFVLLBGX9:4:11401:9125:1052 2:N:0:CGGCTATG+CCGTCGCC\n"
    "NTGCCTCGCACATAGCACCTGCCCCTGTGCAGCCCCCCATGATTCGGCG\n"
    "+\n"
    "#AAAAEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE");

TEST(Tokenizer, TooMuchBuffer)
{
  std::stringstream stm(TWO_RECORDS);
  // read whole record minus last newline
  Tokenizer t(stm, TWO_RECORDS.size() + 10);

  ASSERT_EQ(true, t.token().empty());

  ASSERT_EQ(true, t.next());
  ASSERT_EQ(true, stm.eof());
  ASSERT_EQ(false, t.token().empty());
  ASSERT_EQ(true, t.token().valid());

  ASSERT_EQ(true, t.next());
  ASSERT_EQ(true, stm.eof());
  ASSERT_EQ(false, t.token().empty());
  ASSERT_EQ(true, t.token().valid());

  ASSERT_EQ(false, t.next());
  ASSERT_EQ(true, stm.eof());
  ASSERT_EQ(true, t.token().empty());
  ASSERT_EQ(false, t.token().valid());
}

TEST(Tokenizer, EnoughBufferForOne)
{
  std::stringstream stm(TWO_RECORDS);
  // read whole record minus last newline
  Tokenizer t(stm, TWO_RECORDS.size() - 10);

  ASSERT_EQ(true, t.token().empty());
  ASSERT_EQ(true, t.next());
  ASSERT_EQ(false, stm.eof());
  ASSERT_EQ(false, t.token().empty());
  ASSERT_EQ(true, t.token().valid());

  ASSERT_EQ(true, t.next());
  ASSERT_EQ(true, stm.eof());
  ASSERT_EQ(false, t.token().empty());
  ASSERT_EQ(true, t.token().valid());

  ASSERT_EQ(false, t.next());
  ASSERT_EQ(true, stm.eof());
  ASSERT_EQ(true, t.token().empty());
  ASSERT_EQ(false, t.token().valid());
}

const std::string MIXED_NEWLINE(
    "@NB551322:14:HFVLLBGX9:4:11401:24054:1050 2:N:0:CGGCTATG+CCGTCGCC\n\r"
    "NTGTCGGGGCAGGCAGGGCTCCTCGGGCAGCGGCTCATGAGAGAAGACGGAATCCTCCCCTGAGGAGCACGTAGAGCTCCGGGTGTCGGGAAAGCTGGGGG\r\n"
    "+\n"
    "#AAAAEEEEEEEEEEEEEEEEEEEEEEE6EEEEEEEEEEEEEEEEEAEAEEEEEEEEEEEEAEEE/EEEEAEEEAEEEEEEEEEEEEEEEEEEEEE/EEEE\n"
    "@NB551322:14:HFVLLBGX9:4:11401:9125:1052 2:N:0:CGGCTATG+CCGTCGCC\r"
    "NTGCCTCGCACATAGCACCTGCCCCTGTGCAGCCCCCCATGATTCGGCG\n"
    "+\n\r"
    "#AAAAEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE");

TEST(Tokenizer, MixedNewline)
{
  std::stringstream stm(MIXED_NEWLINE);
  // buffer size does not matter here much
  Tokenizer t(stm, MIXED_NEWLINE.size() - 10, true);

  ASSERT_EQ(true, t.token().empty());
  ASSERT_EQ(true, t.next());
  ASSERT_EQ(false, stm.eof());
  ASSERT_EQ(false, t.token().empty());
  ASSERT_EQ(true, t.token().valid());

  ASSERT_EQ(true, t.next());
  ASSERT_EQ(true, stm.eof());
  ASSERT_EQ(false, t.token().empty());
  ASSERT_EQ(true, t.token().valid());

  ASSERT_EQ(false, t.next());
  ASSERT_EQ(true, stm.eof());
  ASSERT_EQ(true, t.token().empty());
  ASSERT_EQ(false, t.token().valid());
}

const std::string ZERO_LENGTH_READ(
    "@NB551322:14:HFVLLBGX9:4:11401:24054:1050 2:N:0:CGGCTATG+CCGTCGCC\n"
    "\n"
    "+\n"
    "\n"
    "@NB551322:14:HFVLLBGX9:4:11401:9125:1052 2:N:0:CGGCTATG+CCGTCGCC\n"
    "NTGCCTCGCACATAGCACCTGCCCCTGTGCAGCCCCCCATGATTCGGCG\n"
    "+\n"
    "#AAAAEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE");

TEST(Tokenizer, ZeroLength)
{
  std::stringstream stm(ZERO_LENGTH_READ);
  // buffer size does not matter here much
  Tokenizer t(stm, ZERO_LENGTH_READ.size() - 10);

  ASSERT_EQ(true, t.token().empty());
  ASSERT_EQ(true, t.next());
  ASSERT_EQ(false, stm.eof());
  ASSERT_EQ(false, t.token().empty());
  ASSERT_EQ(true, t.token().valid());
  ASSERT_EQ(0, t.token().readLength());

  ASSERT_EQ(true, t.next());
  ASSERT_EQ(true, stm.eof());
  ASSERT_EQ(false, t.token().empty());
  ASSERT_EQ(true, t.token().valid());
  ASSERT_NE(0, t.token().readLength());

  ASSERT_EQ(false, t.next());
  ASSERT_EQ(true, stm.eof());
  ASSERT_EQ(true, t.token().empty());
  ASSERT_EQ(false, t.token().valid());
}

TEST(Tokenizer, EmptyInput)
{
  std::stringstream stm("");
  // read whole record minus last newline
  Tokenizer t(stm, 10);

  ASSERT_EQ(true, t.token().empty());
  ASSERT_EQ(false, t.next());
  ASSERT_EQ(true, stm.eof());
  ASSERT_EQ(true, t.token().empty());
  ASSERT_EQ(false, t.token().valid());
  ASSERT_EQ(0, t.token().readLength());
}
