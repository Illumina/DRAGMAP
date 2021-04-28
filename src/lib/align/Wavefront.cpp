/**
 ** DRAGEN Open Source Software
 ** Copyright (c) 2019-2020 Illumina, Inc.
 ** All rights reserved.
 **
 ** This software is provided under the terms and conditions of the
 ** GNU GENERAL PUBLIC LICENSE Version 3
 **
 ** You should have received a copy of the GNU GENERAL PUBLIC LICENSE Version 3
 ** along with this program. If not, see
 ** <https://github.com/illumina/licenses/>.
 **
 **/

#include "align/Wavefront.hpp"
#include <assert.h>

namespace dragenos {
namespace align {

template <typename T, int WIDTH, int ALIGN>
void WavefrontT<T, WIDTH, ALIGN>::reset()
{
  next_  = 0;
  moved_ = Motion::right;  // doesn't matter
  e_     = decltype(e_)();
  f_     = decltype(f_)();
  h_     = decltype(h_)();
}

template <typename T, int WIDTH, int ALIGN>
const typename WavefrontT<T, WIDTH, ALIGN>::Antidiagonal& WavefrontT<T, WIDTH, ALIGN>::moveRight(
    const Antidiagonal& similarities, const Int gapInit, const Int gapExtend)
{
  resetBs();
  moveRightE(gapInit, gapExtend);
  moveRightF(gapInit, gapExtend);
  moveRightH(similarities);
  moved_ = right;
  return setNextToMax();
}

template <typename T, int WIDTH, int ALIGN>
void WavefrontT<T, WIDTH, ALIGN>::resetBs()
{
  for (unsigned i = 0; WIDTH > i; ++i) {
    bs_[i] = none;
  }
}

template <typename T, int WIDTH, int ALIGN>
const typename WavefrontT<T, WIDTH, ALIGN>::Antidiagonal& WavefrontT<T, WIDTH, ALIGN>::moveDown(
    const Antidiagonal& similarities, const Int gapInit, const Int gapExtend)
{
  resetBs();
  moveDownE(gapInit, gapExtend);
  moveDownF(gapInit, gapExtend);
  moveDownH(similarities);
  moved_ = down;
  return setNextToMax();
}

template <typename T, int WIDTH, int ALIGN>
const typename WavefrontT<T, WIDTH, ALIGN>::Antidiagonal& WavefrontT<T, WIDTH, ALIGN>::setNextToMax()
{
  auto const& nextE = e_[next_];
  auto const& nextF = f_[next_];
  auto&       nextH = h_[next_];

  Antidiagonal h;

  for (unsigned i = 0; WIDTH > i; ++i) {
    h[i] = std::max<T>(0, std::max(std::max(nextE[i], nextF[i]), nextH[i]));

    if (h[i] == nextE[i]) {
      bs_[i] = Backstep(bs_[i] | horz);
    }

    if (h[i] == nextF[i]) {
      bs_[i] = Backstep(bs_[i] | vert);
    }

    if (h[i] == nextH[i]) {
      bs_[i] = Backstep(bs_[i] | diag);
    }

    nextH[i] = h[i];
    if (!nextH[i]) {
      bs_[i] = none;
      ;
    }
  }

  next_ = (next_ + 1) % SIZE;
  return nextH;
}

template <typename T, int WIDTH, int ALIGN>
bool WavefrontT<T, WIDTH, ALIGN>::selectBest(const T extend, const T open, T& ret)
{
  if (open > extend) {
    ret = open;
  } else {
    ret = extend;
    return true;
  }

  return false;
}

// Ei,j = max(Ei,j-1 - gapExtend, Hi,j-1 - gapInit): tracks deletions
// when moving right (constant i), the offset is 0 between the successive antidiagonals
template <typename T, int WIDTH, int ALIGN>
const typename WavefrontT<T, WIDTH, ALIGN>::Antidiagonal& WavefrontT<T, WIDTH, ALIGN>::moveRightE(
    const Int gapInit, const Int gapExtend)
{
  auto&       nextE = e_[next_];
  auto const& lastE = e_[last()];
  auto const& lastH = h_[last()];
  for (unsigned i = 0; nextE.size() > i; ++i) {
    if (selectBest(lastE[i] - gapExtend, lastH[i] - gapInit, nextE[i])) {
      bs_[i] = Backstep(bs_[i] | extHFlag);
    }
  }
  return nextE;
}

// Ei,j = max(Ei,j-1 - gapExtend, Hi,j-1 - gapInit): tracks deletions
// when moving down (increasing i), the offset is 1 between the successive antidiagonals
// Bottom-left element initialized from non-value. Others are shifted
template <typename T, int WIDTH, int ALIGN>
const typename WavefrontT<T, WIDTH, ALIGN>::Antidiagonal& WavefrontT<T, WIDTH, ALIGN>::moveDownE(
    const Int gapInit, const Int gapExtend)
{
  auto&       nextE = e_[next_];
  auto const& lastE = e_[last()];
  auto const& lastH = h_[last()];
  nextE.front()     = -1;
  for (unsigned i = 1; nextE.size() > i; ++i) {
    if (selectBest(lastE[i - 1] - gapExtend, lastH[i - 1] - gapInit, nextE[i])) {
      bs_[i] = Backstep(bs_[i] | extHFlag);
    }
  }
  return nextE;
}

// Fi,j = max(Fi-1,j - gapExtend, Hi-1,j - gapInit): tracks insertions
// when moving right (increasing j), the offset is 1 between the successive antidiagonals
// Top-right element initialized from non-value. Others are shifted
template <typename T, int WIDTH, int ALIGN>
const typename WavefrontT<T, WIDTH, ALIGN>::Antidiagonal& WavefrontT<T, WIDTH, ALIGN>::moveRightF(
    const Int gapInit, const Int gapExtend)
{
  auto&       nextF = f_[next_];
  auto const& lastF = f_[last()];
  auto const& lastH = h_[last()];
  for (unsigned i = 0; nextF.size() > i + 1; ++i) {
    if (selectBest(lastF[i + 1] - gapExtend, lastH[i + 1] - gapInit, nextF[i])) {
      bs_[i] = Backstep(bs_[i] | extVFlag);
    }
  }
  nextF.back() = -1;
  return nextF;
}

// Fi,j = max(Fi-1,j - gapExtend, Hi-1,j - gapInit): tracks insertions
// when moving down (increasing i), the offset is 1 between the successive antidiagonals
template <typename T, int WIDTH, int ALIGN>
const typename WavefrontT<T, WIDTH, ALIGN>::Antidiagonal& WavefrontT<T, WIDTH, ALIGN>::moveDownF(
    const Int gapInit, const Int gapExtend)
{
  auto&       nextF = f_[next_];
  auto const& lastF = f_[last()];
  auto const& lastH = h_[last()];
  for (unsigned i = 0; nextF.size() > i; ++i) {
    if (selectBest(lastF[i] - gapExtend, lastH[i] - gapInit, nextF[i])) {
      bs_[i] = Backstep(bs_[i] | extVFlag);
    }
  }
  return nextF;
}

// Hi,j = Hi-1,j-1 + similarity(query, database): tracks mismatches
// Hi-1,j-1 is the penultimate antidiagonal and the previous motion (moved_) is needed
// when moved_==down, the overall motion was diagonal and there is no offset
// when moved_==right, the overall motion is right and there is an offset of 1
template <typename T, int WIDTH, int ALIGN>
const typename WavefrontT<T, WIDTH, ALIGN>::Antidiagonal& WavefrontT<T, WIDTH, ALIGN>::moveRightH(
    const Antidiagonal& similarities)
{
  auto& nextH              = h_[next_];
  nextH                    = similarities;
  auto const& penultimateH = h_[penultimate()];
  if (down == moved_)  // down then right
  {
    for (unsigned i = 0; nextH.size() > i; ++i) {
      nextH[i] += penultimateH[i];
    };
  } else  // two consecutive rights - shift by -1
  {
    for (unsigned i = 0; nextH.size() > i + 1; ++i) {
      nextH[i] += penultimateH[i + 1];
    }
    // nextH.back() keeps the appropriate similarity value
  }
  return nextH;
}

// Hi,j = Hi-1,j-1 + similarity(query, database): tracks mismatches
// Hi-1,j-1 is the penultimate antidiagonal and the previous motion (moved_) is needed
// when moved_==right, the overall motion was diagonal and there is no offset
// when moved_==down, the overall motion is down and there is an offset of 1
template <typename T, int WIDTH, int ALIGN>
const typename WavefrontT<T, WIDTH, ALIGN>::Antidiagonal& WavefrontT<T, WIDTH, ALIGN>::moveDownH(
    const Antidiagonal& similarities)
{
  auto& nextH              = h_[next_];
  nextH                    = similarities;
  auto const& penultimateH = h_[penultimate()];
  if (right == moved_)  // right then down
  {
    for (unsigned i = 0; nextH.size() > i; ++i) {
      nextH[i] += penultimateH[i];
    };
  } else  // two consecutive down - shift by +1
  {
    // nextH[0] keeps the appropriate similarity value
    for (unsigned i = 1; nextH.size() > i; ++i) {
      nextH[i] += penultimateH[i - 1];
    }
  }
  return nextH;
}

template class WavefrontT<short, 48, 16>;
// for tests
template class WavefrontT<short, 8, 16>;
}  // namespace align
}  // namespace dragenos
