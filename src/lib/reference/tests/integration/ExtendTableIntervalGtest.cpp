#include "gtest/gtest.h"

#include "reference/ExtendTableInterval.hpp"

TEST(ExtendTableInterval, SetOfOne)
{
  // SLO
  // record 7bf9df58:f8030726 Start = 0x726 Length = 6
  // SLE(eflifts > 0, MSB=0)
  // record 7bf9df58:f9433726 exlifts = 0x43 Start = 0x26 Length = 0x37
  using dragenos::reference::ExtendTableInterval;
  using dragenos::reference::HashRecord;
  uint64_t s1       = 0x7bf9df58f8030726;
  uint64_t s2       = 0x7bf9df58f9433726;
  uint64_t ss[]     = {s1, s2};
  uint32_t esl[][3] = {{0, 0x726, 6}, {0x43, 0x26, 0x37}};
  for (unsigned i = 0; 2 > i; ++i) {
    const auto          begin = reinterpret_cast<const HashRecord*>(&(ss[i]));
    const auto          end   = begin + 1;
    ExtendTableInterval interval(begin, end);
    ASSERT_EQ(esl[i][0], interval.getExtraLiftoverMatchCount());
    ASSERT_EQ(esl[i][1], interval.getStart());
    ASSERT_EQ(esl[i][2], interval.getLength());
  }
}

TEST(ExtendTableInterval, SetOfTwo)
{
  // SL1, S
  // 0x7bf9df59f8123456, 0x7bf9fedcfa987654 exlifts = 0    start = 0x56987654 length = 0x1234
  // SLE(exlifts > 0, MSB=1), S
  // 0x7bf9df59f9123456, 0x7bf9fedcfa987654 exlifts = 0x12 start = 0x56987654 length = 0x34
  // SLE(exlifts = 0, MSB=0), L
  // 0x7bf9df58f9003456, 0x7bf9fedcfb987654 exlifts = 0    start = 0x56       length = 0x34987654
  // S, L
  // 0x7bf9df58fa123456, 0x7bf9fedcfb987654 exlifts = 0    start = 0x123456   length = 0x987654
  using dragenos::reference::ExtendTableInterval;
  using dragenos::reference::HashRecord;
  uint64_t  s1[]     = {0x7bf9df59f8123456, 0x7bf9fedcfa987654};
  uint64_t  s2[]     = {0x7bf9df59f9123456, 0x7bf9fedcfa987654};
  uint64_t  s3[]     = {0x7bf9df58f9003456, 0x7bf9fedcfb987654};
  uint64_t  s4[]     = {0x7bf9df58fa123456, 0x7bf9fedcfb987654};
  uint64_t* ss[]     = {s1, s2, s3, s4};
  uint32_t  esl[][3] = {
      {0, 0x56987654, 0x1234}, {0x12, 0x56987654, 0x34}, {0, 0x56, 0x34987654}, {0, 0x123456, 0x987654}};
  for (unsigned i = 0; 4 > i; ++i) {
    const auto begin = reinterpret_cast<const HashRecord*>(ss[i]);
    const auto end   = begin + 2;
    ASSERT_EQ(2, end - begin);
    ExtendTableInterval interval(begin, end);
    ASSERT_EQ(esl[i][0], interval.getExtraLiftoverMatchCount()) << "i: " << i;
    ASSERT_EQ(esl[i][1], interval.getStart()) << "i: " << i;
    ASSERT_EQ(esl[i][2], interval.getLength()) << "i: " << i;
  }
}

TEST(ExtendTableInterval, SetOfThree)
{
  // Only SLE + S + L but with different MSB and Carry flags
  // record: 7bf9df58:f2300726
  // record: 7bf9df59:f9000000: SLE exlift: 0 length: 0 start: 0
  // record: 7bf9df58:fa582838: S   start: 5777464
  // record: 7bf9df5a:fb0168c4: L   length: 92356
  //
  // Second example:
  // record: 044d5968:f2280559
  // record: 044d5969:f9000000: SLE exlift: 0 length: 0 start: 0
  // record: 044d5968:fa4092ac: S   start: 4231852
  // record: 044d596a:fb013a82: L   length: 80514
  //
  // Third example:
  // record: b33ac730:f2301041
  // record: b33ac731:f9000001: SLE exlift: 0 length: 0 start: 1
  // record: b33ac731:fa7c4fc2: S   start: 8146882
  // record: b33ac732:fb01a0a2: L   length: 106658
  //
  // Fourth example:
  // record: 9c50a1d8:f230116d
  // record: 9c50a1d9:f900001a: SLE exlift: 0 length: 0 start: 26
  // record: 9c50a1d9:fa46391e: S   start: 4602142
  // record: 9c50a1da:fb017940: L   length: 96576
  //
  // Fifth example:
  // record: c851b148:f2300595
  // record: c851b149:f900001a: SLE exlift: 0 length: 0 start: 26
  // record: c851b148:fabb10ae: S   start: 12259502
  // record: c851b14a:fb01f3b0: L   length: 127920
  //
  using dragenos::reference::ExtendTableInterval;
  using dragenos::reference::HashRecord;
  uint64_t  s1[]    = {0x7bf9df59f9000000, 0x7bf9df58fa582838, 0x7bf9df5afb0168c4};
  uint64_t  s2[]    = {0x044d5969f9000000, 0x044d5968fa4092ac, 0x044d596afb013a82};
  uint64_t  s3[]    = {0xb33ac731f9000001, 0xb33ac731fa7c4fc2, 0xb33ac732fb01a0a2};
  uint64_t  s4[]    = {0x9c50a1d9f900001a, 0x9c50a1d9fa46391e, 0x9c50a1dafb017940};
  uint64_t  s5[]    = {0xc851b149f900001a, 0xc851b148fabb10ae, 0xc851b14afb01f3b0};
  uint64_t* ss[]    = {s1, s2, s3, s4, s5};
  uint32_t  sl[][2] = {
      {5777464, 92356}, {4231852, 80514}, {0x27c4fc2, 106658}, {0x1b46391e, 96576}, {0x1abb10ae, 127920}};
  for (unsigned i = 0; 5 > i; ++i) {
    const auto begin = reinterpret_cast<const HashRecord*>(ss[i]);
    const auto end   = begin + 3;
    const auto type0 = begin->getType();
    const auto type1 = (begin + 1)->getType();
    const auto type2 = (begin + 2)->getType();
    ASSERT_EQ(3, end - begin);
    ASSERT_EQ(HashRecord::INTERVAL_SLE, type0);
    ASSERT_EQ(0, begin->getExlift());
    ASSERT_EQ(true, begin->isMsb());
    ASSERT_EQ(HashRecord::INTERVAL_S, type1);
    ASSERT_EQ(HashRecord::INTERVAL_L, type2);
    ExtendTableInterval interval(begin, end);
    ASSERT_EQ(sl[i][0], interval.getStart()) << "i: " << i;
    ASSERT_EQ(sl[i][1], interval.getLength()) << "i: " << i;
    ASSERT_EQ(0, interval.getExtraLiftoverMatchCount()) << "i: " << i;
  }
}
