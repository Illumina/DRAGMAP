#include "gtest/gtest.h"

#include <algorithm>

#include "align/Wavefront.hpp"

TEST(Wavefront, Constructor)
{
  using Antidiagonal = dragenos::align::Antidiagonal48;
  using Wavefront    = dragenos::align::Wavefront48;
  typedef Wavefront::History History;
  const Wavefront            w0;
  for (const auto v : w0.getLastScores()) {
    ASSERT_EQ(v, 0);
  }
  ASSERT_EQ(w0.next(), 0);
  ASSERT_EQ(w0.last(), 2);
  ASSERT_EQ(w0.penultimate(), 1);
  const auto uniformAntidiagonal = [](short init) {
    Antidiagonal tmp;
    for (auto& v : tmp.data) v = init;
    return tmp;
  };
  const Antidiagonal a0 = uniformAntidiagonal(0);
  const Antidiagonal a1 = uniformAntidiagonal(1);
  const Antidiagonal a2 = uniformAntidiagonal(2);
  const Antidiagonal a3 = uniformAntidiagonal(3);
  const Antidiagonal a4 = uniformAntidiagonal(4);
  const Antidiagonal a5 = uniformAntidiagonal(5);
  const History      h2{a0, a1, a2};
  const History      h3{a1, a2, a3};
  const History      h4{a2, a3, a4};
  const History      h5{a3, a4, a5};
  const Wavefront    w1(h2, h3, h4, 1, Wavefront::right);
  ASSERT_EQ(w1.next(), 1);
  ASSERT_EQ(w1.last(), 0);
  ASSERT_EQ(w1.penultimate(), 2);
  const Wavefront w2(h3, h4, h5, 2, Wavefront::down);
  ASSERT_EQ(w2.next(), 2);
  ASSERT_EQ(w2.last(), 1);
  ASSERT_EQ(w2.penultimate(), 0);
}

TEST(Wavefront, moveRightE)
{
  using Wavefront    = dragenos::align::Wavefront48;
  using Antidiagonal = dragenos::align::Antidiagonal48;
  typedef Wavefront::History History;
  constexpr auto             width               = Antidiagonal::width;
  constexpr short            gapInit             = 5;
  constexpr short            gapExtend           = 3;
  const auto                 uniformAntidiagonal = [](short init) {
    Antidiagonal tmp;
    for (auto& v : tmp.data) v = init;
    return tmp;
  };
  const auto hatAntidiagonal = []() {
    Antidiagonal tmp;
    short        i         = 0;
    short        increment = 1;
    for (auto& v : tmp.data) {
      if (1 + Antidiagonal::width / 2 == i) increment = -1;
      v = i;
      i += increment;
    };
    return tmp;
  };
  const Antidiagonal a0 = uniformAntidiagonal(20);
  const Antidiagonal a1 = hatAntidiagonal();
  const History      e{a0};
  const History      f{};
  const History      h{a1};
  Wavefront          w(e, f, h, 1, Wavefront::right);
  const auto         a = w.moveRightE(gapInit, gapExtend);
  // should be 17 == 20 - gapExtend when coming from E -
  for (unsigned int i = 0; 23 > i; ++i) ASSERT_EQ(a[i], 17);
  // should be H - gapInit when coming from H
  for (unsigned int i = 23; 29 > i; ++i) ASSERT_EQ(a[i], h[0][i] - gapInit);
  // should be 17 == 20 - gapExtend when coming from E -
  for (unsigned int i = 29; width > i; ++i) ASSERT_EQ(a[i], 17);
}

TEST(Wavefront, moveDownE)
{
  // same tests as moveRightE but with a +1 shift relative to the last antidiagonal
  using Wavefront    = dragenos::align::Wavefront48;
  using Antidiagonal = dragenos::align::Antidiagonal48;
  typedef Wavefront::History History;
  constexpr auto             width               = Antidiagonal::width;
  constexpr short            gapInit             = 5;
  constexpr short            gapExtend           = 3;
  const auto                 uniformAntidiagonal = [](short init) {
    Antidiagonal tmp;
    for (auto& v : tmp.data) v = init;
    return tmp;
  };
  const auto hatAntidiagonal = []() {
    Antidiagonal tmp;
    short        i         = 0;
    short        increment = 1;
    for (auto& v : tmp.data) {
      if (1 + Antidiagonal::width / 2 == i) increment = -1;
      v = i;
      i += increment;
    };
    return tmp;
  };
  const Antidiagonal a0 = uniformAntidiagonal(20);
  const Antidiagonal a1 = hatAntidiagonal();
  const History      e{a0};
  const History      f{};
  const History      h{a1};
  Wavefront          w(e, f, h, 1, Wavefront::right);
  const auto         a = w.moveDownE(gapInit, gapExtend);
  // because of the shift, the first element should be 0
  ASSERT_EQ(a.front(), -1);
  // should be 17 == 20 - gapExtend when coming from E -
  for (unsigned int i = 1; 24 > i; ++i) ASSERT_EQ(a[i], 17) << "a[i==" << i << "]";
  // should be h[i-1] - gapInit when coming from H -
  for (unsigned int i = 24; 29 > i; ++i) ASSERT_EQ(a[i], h[0][i - 1] - gapInit) << "a[i==" << i << "]";
  // should be 17 == 20 - gapExtend when coming from E -
  for (unsigned int i = 29; width > i; ++i) ASSERT_EQ(a[i], 17) << "a[i==" << i << "]";
}

TEST(Wavefront, moveRightF)
{
  // same tests as moveRightE but using F and
  // with a -1 shift relative to the last antidiagonal
  using Wavefront    = dragenos::align::Wavefront48;
  using Antidiagonal = dragenos::align::Antidiagonal48;
  typedef Wavefront::History History;
  constexpr auto             width               = Antidiagonal::width;
  constexpr short            gapInit             = 5;
  constexpr short            gapExtend           = 3;
  const auto                 uniformAntidiagonal = [](short init) {
    Antidiagonal tmp;
    for (auto& v : tmp.data) v = init;
    return tmp;
  };
  const auto hatAntidiagonal = []() {
    Antidiagonal tmp;
    short        i         = 0;
    short        increment = 1;
    for (auto& v : tmp.data) {
      if (1 + Antidiagonal::width / 2 == i) increment = -1;
      v = i;
      i += increment;
    };
    return tmp;
  };
  const Antidiagonal a0 = uniformAntidiagonal(20);
  const Antidiagonal a1 = hatAntidiagonal();
  const History      e{};
  const History      f{a0};
  const History      h{a1};
  Wavefront          w(e, f, h, 1, Wavefront::right);
  const auto         a = w.moveRightF(gapInit, gapExtend);
  // should be 17 == 20 - gapExtend when coming from E -
  for (unsigned int i = 0; 22 > i; ++i) ASSERT_EQ(a[i], 17) << "a[i==" << i << "]";
  // should be H - gapInit when coming from H
  for (unsigned int i = 22; 28 > i; ++i) ASSERT_EQ(a[i], h[0][i + 1] - gapInit) << "a[i==" << i << "]";
  // should be 17 == 20 - gapExtend when coming from E -
  for (unsigned int i = 28; width > i + 1; ++i) ASSERT_EQ(a[i], 17) << "a[i==" << i << "]";
  // because of the shift the first element should be 0
  ASSERT_EQ(a.back(), -1);
}

TEST(Wavefront, moveDownF)
{
  // same tests as moveRightF but without the gap
  using Wavefront    = dragenos::align::Wavefront48;
  using Antidiagonal = dragenos::align::Antidiagonal48;
  typedef Wavefront::History History;
  constexpr auto             width               = Antidiagonal::width;
  constexpr short            gapInit             = 5;
  constexpr short            gapExtend           = 3;
  const auto                 uniformAntidiagonal = [](short init) {
    Antidiagonal tmp;
    for (auto& v : tmp.data) v = init;
    return tmp;
  };
  const auto hatAntidiagonal = []() {
    Antidiagonal tmp;
    short        i         = 0;
    short        increment = 1;
    for (auto& v : tmp.data) {
      if (1 + Antidiagonal::width / 2 == i) increment = -1;
      v = i;
      i += increment;
    };
    return tmp;
  };
  const Antidiagonal a0 = uniformAntidiagonal(20);
  const Antidiagonal a1 = hatAntidiagonal();
  const History      e{};
  const History      f{a0};
  const History      h{a1};
  Wavefront          w(e, f, h, 1, Wavefront::right);
  const auto         a = w.moveDownF(gapInit, gapExtend);
  // should be 17 == 20 - gapExtend when coming from E -
  for (unsigned int i = 0; 23 > i; ++i) ASSERT_EQ(a[i], 17);
  // should be H - gapInit when coming from H
  for (unsigned int i = 23; 29 > i; ++i) ASSERT_EQ(a[i], h[0][i] - gapInit);
  // should be 17 == 20 - gapExtend when coming from E -
  for (unsigned int i = 29; width > i; ++i) ASSERT_EQ(a[i], 17);
}

TEST(Wavefront, moveRightH)
{
  // requires an history of 2 H, a query and a database
  using Wavefront    = dragenos::align::Wavefront48;
  using Antidiagonal = dragenos::align::Antidiagonal48;
  typedef Wavefront::History History;
  constexpr int              MATCH        = 2;
  constexpr short            MISMATCH     = -3;
  const Antidiagonal         similarities = []() {
    Antidiagonal tmp;
    for (size_t i = 0; tmp.size() > i; ++i) tmp[i] = MATCH * (0 == (i % 2)) + MISMATCH * (1 == (i % 2));
    return tmp;
  }();
  // slowly incrementing from 0 to ensure that some values are smaller than the mismatch penalty
  const auto antidiagonal = []() {
    Antidiagonal tmp;
    for (size_t i = 0; tmp.data.size() > i; ++i) {
      tmp.data[i] = i / 8;
    };
    return tmp;
  };
  const Antidiagonal a1 = antidiagonal();
  const History      e{};
  const History      f{};
  const History      h{a1};
  {
    Wavefront  fromRight(e, f, h, 2, Wavefront::right);
    const auto a = fromRight.moveRightH(similarities);
    for (unsigned i = 0; a.size() > i + 1; ++i) {
      ASSERT_EQ(a[i], a1[i + 1] + similarities[i]) << "a[i==" << i << "]";
    }
    ASSERT_EQ(a.back(), MISMATCH);  // odd positions are MISMATCH
  }
  {
    Wavefront  fromDown(e, f, h, 2, Wavefront::down);
    const auto a = fromDown.moveRightH(similarities);
    for (unsigned i = 0; a.size() > i; ++i) {
      ASSERT_EQ(a[i], a1[i] + similarities[i]) << "a[i==" << i << "]";
    }
  }
}

TEST(Wavefront, moveDownH)
{
  // requires an history of 2 H, a query and a database
  using Wavefront    = dragenos::align::Wavefront48;
  using Antidiagonal = dragenos::align::Antidiagonal48;
  typedef Wavefront::History History;
  constexpr int              MATCH        = 2;
  constexpr short            MISMATCH     = -3;
  const Antidiagonal         similarities = []() {
    Antidiagonal tmp;
    for (size_t i = 0; tmp.size() > i; ++i) tmp[i] = MATCH * (0 == (i % 2)) + MISMATCH * (1 == (i % 2));
    return tmp;
  }();
  // slowly incrementing from 0 to ensure that some values are smaller than the mismatch penalty
  const auto antidiagonal = []() {
    Antidiagonal tmp;
    for (size_t i = 0; tmp.data.size() > i; ++i) {
      tmp.data[i] = i / 8;
    };
    return tmp;
  };
  const Antidiagonal a1 = antidiagonal();
  const History      e{};
  const History      f{};
  const History      h{a1};
  {
    Wavefront  fromRight(e, f, h, 2, Wavefront::right);
    const auto a = fromRight.moveDownH(similarities);
    for (unsigned i = 0; a.size() > i; ++i) {
      ASSERT_EQ(a[i], a1[i] + similarities[i]) << "a[i==" << i << "]";
    }
  }
  {
    Wavefront  fromDown(e, f, h, 2, Wavefront::down);
    const auto a = fromDown.moveDownH(similarities);
    ASSERT_EQ(a.front(), MATCH);  // even positions are MATCH
    for (unsigned i = 1; a.size() > i; ++i) {
      ASSERT_EQ(a[i], a1[i - 1] + similarities[i]) << " with i=" << i;
    }
  }
}

TEST(Wavefront, moveRight)
{
  using Wavefront    = dragenos::align::Wavefront48;
  using Antidiagonal = dragenos::align::Antidiagonal48;
  typedef Wavefront::History History;
  constexpr short            MISMATCH  = -3;
  constexpr short            gapInit   = 5;
  constexpr short            gapExtend = 3;
  const Antidiagonal         similarities{MISMATCH,
                                  MISMATCH,
                                  MISMATCH,
                                  MISMATCH,
                                  MISMATCH,
                                  MISMATCH,
                                  MISMATCH,
                                  MISMATCH,
                                  MISMATCH,
                                  MISMATCH,
                                  MISMATCH,
                                  MISMATCH,
                                  MISMATCH,
                                  MISMATCH,
                                  MISMATCH};
  // E is unshifted
  const Antidiagonal ae{
      gapExtend - 1, gapExtend + 1, gapExtend + 2, gapExtend + 3, gapExtend + 4, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  // F is shifted by one to the right (nextF[i] = lastF[i-1]
  const Antidiagonal af{-2, -2, -2, gapExtend + 5, gapExtend + 6, 0, gapExtend + 7, 0, 0, 0, 0};
  const Antidiagonal a0{};
  const History      e{a0, ae};
  const History      f{a0, af};
  {
    // H is shifted by one to the left (nextH[i] = penultimateH[i+1])
    const Antidiagonal ah{0, 0, 0, 0, 0, 8 - MISMATCH, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    const History      h{ah};
    Wavefront          fromRight(e, f, h, 2, Wavefront::right);
    const auto         a = fromRight.moveRight(similarities, gapInit, gapExtend);
    ASSERT_EQ(a[0], 0);  // max(0, nextE==-1)
    ASSERT_EQ(a[1], 1);  // max(0, nextE==+1)
    ASSERT_EQ(a[2], 5);  // max(max(0, nextE==+2), nextF==5)
    ASSERT_EQ(a[3], 6);  // max(max(0, nextE==+3), nextF==6)
    ASSERT_EQ(a[4], 8);  // max(max(max(0, nextE==+4), nextF==0),nextH==8)
  }
  {
    // H is unshifted
    const Antidiagonal ah{0, 0, 0, 0, 8 - MISMATCH, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    const History      h{ah};
    Wavefront          fromDown(e, f, h, 2, Wavefront::down);
    const auto         a = fromDown.moveRight(similarities, gapInit, gapExtend);
    ASSERT_EQ(a[0], 0);  // max(0, nextE==-1)
    ASSERT_EQ(a[1], 1);  // max(0, nextE==+1)
    ASSERT_EQ(a[2], 5);  // max(max(0, nextE==+2), nextF==5)
    ASSERT_EQ(a[3], 6);  // max(max(0, nextE==+3), nextF==6)
    ASSERT_EQ(a[4], 8);  // max(max(max(0, nextE==+4), nextF==0),nextH==8)
  }
}

TEST(Wavefront, moveDown)
{
  using Wavefront    = dragenos::align::Wavefront48;
  using Antidiagonal = dragenos::align::Antidiagonal48;
  typedef Wavefront::History History;
  constexpr short            MISMATCH  = -3;
  constexpr short            gapInit   = 5;
  constexpr short            gapExtend = 3;
  const Antidiagonal         similarities{MISMATCH,
                                  MISMATCH,
                                  MISMATCH,
                                  MISMATCH,
                                  MISMATCH,
                                  MISMATCH,
                                  MISMATCH,
                                  MISMATCH,
                                  MISMATCH,
                                  MISMATCH,
                                  MISMATCH,
                                  MISMATCH,
                                  MISMATCH,
                                  MISMATCH,
                                  MISMATCH};
  // E is shifted by one yo the right (nextE[i] = lastE[i-1])
  const Antidiagonal ae{
      gapExtend - 1, gapExtend + 1, gapExtend + 2, gapExtend + 3, gapExtend + 4, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  // F is unshifted
  const Antidiagonal af{0, 0, 0, gapExtend + 5, gapExtend + 6, 0, gapExtend + 7, 0, 0, 0, 0};
  const Antidiagonal a0{};
  const History      e{a0, ae};
  const History      f{a0, af};
  {
    // H is unshifted
    const Antidiagonal ah{0, 0, 0, 0, 0, 8 - MISMATCH, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    const History      h{ah};
    Wavefront          fromRight(e, f, h, 2, Wavefront::right);
    const auto         a = fromRight.moveDown(similarities, gapInit, gapExtend);
    ASSERT_EQ(a[0], 0);  // max(0, nextE==0)
    ASSERT_EQ(a[1], 0);  // max(0, nextE==-1)
    ASSERT_EQ(a[2], 1);  // max(0, nextE==+1)
    ASSERT_EQ(a[3], 5);  // max(max(0, nextE==+2), nextF==5)
    ASSERT_EQ(a[4], 6);  // max(max(0, nextE==+3), nextF==6)
    ASSERT_EQ(a[5], 8);  // max(max(max(0, nextE==+4), nextF==0),nextH==8)
  }
  {
    // H is shifted by 1 to the right (nextH[i]=penultimateH[i-1])
    const Antidiagonal ah{0, 0, 0, 0, 8 - MISMATCH, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    const History      h{ah};
    Wavefront          fromDown(e, f, h, 2, Wavefront::down);
    const auto         a = fromDown.moveDown(similarities, gapInit, gapExtend);
    ASSERT_EQ(a[0], 0);  // max(0, nextE==0)
    ASSERT_EQ(a[0], 0);  // max(0, nextE==-1)
    ASSERT_EQ(a[2], 1);  // max(0, nextE==+1)
    ASSERT_EQ(a[3], 5);  // max(max(0, nextE==+2), nextF==1)
    ASSERT_EQ(a[4], 6);  // max(max(0, nextE==+2), nextF==3)
    ASSERT_EQ(a[5], 8);  // max(max(max(0, nextE==+2), nextF==3),nextH==4)
  }
}
