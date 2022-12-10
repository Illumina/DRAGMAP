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

#include "align/SmithWaterman.hpp"
#include "common/DragenLogger.hpp"

#include <bitset>
#include <cassert>
#include <numeric>

namespace dragenos {
namespace align {

/**
 *   \param forcedHorizontalMotion distance between bottom-left of antidiagonal and database begin
 */
template <typename C, typename T, int WIDTH, int ALIGN, unsigned STEERING_DELAY>
void SmithWatermanT<C, T, WIDTH, ALIGN, STEERING_DELAY>::buildWavefronts(
    const size_t forcedVerticalMotion, const size_t forcedDiagonalMotion, const size_t forcedHorizontalMotion)
{
  DRAGEN_S_W_STEER_LOG << "--------" << std::endl;
  forcedHorizontalMotion_ = forcedHorizontalMotion;
  forcedVerticalMotion_   = forcedVerticalMotion;
  forcedDiagonalMotion_   = forcedDiagonalMotion;
  autoSteerEnabled_       = false;

  // adaptive steering until a move that would go past the end of query or database
  drift_ = 0;  // initialized to be centered on the expected diagonal
  while (true) {
    const auto move = getNextMove();
    if (Motion::right == move) {
      if (getDatabaseOffset() >= getDatabaseSize() - 1) {
        break;
      } else {
        moveRight();
      }
    } else  // move down
    {
      if (queryIt_ <= query_.begin() + width - 1) {
        break;
      } else {
        assert(queryIt_ > query_.begin() + width - 1);
        moveDown();
      }
    }
  }

  // dragen thing...
  setDatabaseSize(std::min(getDatabaseOffset() + 1 + 2 * width + 1, getDatabaseSize()));
  while (getDatabaseOffset() < getDatabaseSize() - 1) {
    moveRight();
  }

  // move down all the way past the end of the actual query unless we hit the bottom
  while (queryIt_ != query_.begin() && (getQueryOffset() < getQuerySize() - 1))
  //    || (getQueryOffset() + 1 - getQuerySize()) < (getDatabaseSize() - getDatabaseOffset() - 1)))
  {
    moveDown();
  }
}

template <typename C, typename T, int WIDTH, int ALIGN, unsigned STEERING_DELAY>
T SmithWatermanT<C, T, WIDTH, ALIGN, STEERING_DELAY>::buildCigar(
    const std::size_t querySize, std::string& cigar)
{
  //    std::cerr << "backtracking from global max score " << *maxIndices_.at(globalMax_) << " globalMax_=" << globalMax_ << scores_.at(globalMax_) << std::endl;
  return buildCigar(querySize, globalMax_, globalMaxOffset_, cigar);
}

/**
 * \param i       current antidiagonal to start backtracking
 * \param offset  offset in the current antidiagonal
 */
template <typename C, typename T, int WIDTH, int ALIGN, unsigned STEERING_DELAY>
T SmithWatermanT<C, T, WIDTH, ALIGN, STEERING_DELAY>::buildCigar(
    const std::size_t querySize, int i, int offset, std::string& ret)
{
  ret.clear();
  static const std::array<char, 8> operations = {
      '?',
      'D',  // deletion
      'I',  // insertion
      '?',
      'M',  // map
      '?',
      '?',
      '?',
  };

  ///  std::cerr << "i " << i << std::endl;
  const int rightDrift = std::count(motions_.begin(), motions_.begin() + i + 1, Motion::right);
  //  std::cerr << "rightDrift " << rightDrift << std::endl;
  // building CIGAR backwards. Compute end offset
  const int hardClipEnd = std::distance(databaseBeginIt_, databaseEndIt_) - rightDrift + width - offset - 1;
  //  std::cerr << "hardClipEnd " << hardClipEnd << std::endl;
  const int downDrift = width + std::count(motions_.begin(), motions_.begin() + i + 1, Motion::down);
  //  std::cerr << "downDrift " << downDrift << std::endl;
  //  std::cerr << "offset " << offset << std::endl;
  //  std::cerr << "querySize " << querySize << std::endl;
  int queryPos = downDrift - offset;
  //  std::cerr << "queryPos " << queryPos << std::endl;
  const int softClipEnd = querySize - queryPos;
  //  std::cerr << "softClipEnd " << softClipEnd << std::endl;
  assert(0 <= softClipEnd);

  Backstep extFlag = Backstep::none;
  for (; 0 < i && 0 <= offset && queryPos && WIDTH > offset; --i) {
    //     std::cerr << "i:" << i << " offset:" << offset << " score:" << scores_.at(i)[offset] << std::endl;
    //    std::cerr << "bs.size()=" << bs_.size() << " offset=" << offset << " bs_.at(i).size()=" << bs_.at(i).size() << " bs_.at(i).at(offset)=" << bs_.at(i).at(offset) << std::endl;
    Backstep bs = bs_.at(i).at(offset);
    if (extFlag & Backstep::extHFlag) {
      extFlag = Backstep(bs & Backstep::extHFlag);
      bs      = Backstep::horz;
    } else if (extFlag & Backstep::extVFlag) {
      extFlag = Backstep(bs & Backstep::extVFlag);
      bs      = Backstep::vert;
    } else if (bs & Backstep::diag) {
      bs = Backstep::diag;
    } else if (bs & Backstep::vert) {
      extFlag = Backstep(bs & Backstep::extVFlag);
      bs      = Backstep::vert;
    } else if (bs & Backstep::horz) {
      extFlag = Backstep(bs & Backstep::extHFlag);
      bs      = Backstep::horz;
    } else {
      // nowhere to go. Done backtracking
      //      std::cerr << "nowhere to go at i=" << i << " offset=" << offset << " ad: " << scores_.at(i) << std::endl;
      break;
    }

    ret += operations.at(bs);
    //    std::cerr << ret.back() << std::endl;

    // offset on the next antidiagonal if we transition as deletion. Insertion offset is delOffset + 1 ...
    offset -= (Motion::down == motions_.at(i));

    if (Backstep::horz != bs) {
      ++offset;
      --queryPos;
    }

    if (Backstep::diag == bs) {
      // next offset on the diagonal if we're match or mismatch. i never equals 0, so this is safe
      offset -= (Motion::down == motions_.at(i - 1));
      // when moving diagonally, skip one
      --i;
    }
  }

  // calculate clipping
  if (!i && width - 1 == offset) {
    // reached the corner.
    ret += operations.at(Backstep::diag);
  }

  //   std::cerr << "ret " << ret << std::endl;
  //   std::cerr << "i " << i << std::endl;
  //   std::cerr << "offset " << offset << std::endl;

  const std::size_t maps       = std::count(ret.begin(), ret.end(), 'M');
  const std::size_t insertions = std::count(ret.begin(), ret.end(), 'I');
  const std::size_t deletions  = std::count(ret.begin(), ret.end(), 'D');
  //   std::cerr << "maps " << maps << " insertions " << insertions << " deletions " << deletions << std::endl;

  const std::size_t coveredBases = maps + insertions;
  //  const int softClipEnd = querySize - coveredBases - softClipStart;
  const int softClipStart = querySize - coveredBases - softClipEnd;
  //  std::cerr << "softClipStart " << softClipStart << std::endl;

  //  std::cerr << "coveredBases " << coveredBases << std::endl;
  //  std::cerr << "softClipEnd " << softClipEnd << std::endl;
  assert(0 <= softClipStart);

  //  ret = std::string(hardClipEnd, 'N') +
  //    std::string(softClipEnd, 'S') + ret +
  //    std::string(softClipStart, 'S');
  const int hardClipStart = std::distance(databaseBeginIt_, databaseEndIt_) - maps - deletions - hardClipEnd;
  std::reverse(ret.begin(), ret.end());
  ret =
      std::string(hardClipStart, 'N') + std::string(softClipStart, 'S') + ret + std::string(softClipEnd, 'S');

  return getMaxScore() - (softClipStart ? 0 : unclipScore_) - (softClipEnd ? 0 : unclipScore_);
}

template <typename C, typename T, int WIDTH, int ALIGN, unsigned STEERING_DELAY>
T SmithWatermanT<C, T, WIDTH, ALIGN, STEERING_DELAY>::align(
    const C*     queryBegin,
    const C*     queryEnd,
    const C*     databaseBegin,
    const C*     databaseEnd,
    const int    forcedDiagonalMotion,
    const int    forcedHorizontalMotion,
    bool         reverseQuery,
    std::string& cigar)
{
  const int querySize = std::distance(queryBegin, queryEnd);
  const int dbSize    = std::distance(databaseBegin, databaseEnd);

  reset(queryBegin, queryEnd, databaseBegin, databaseEnd, reverseQuery);

  const int w = width;  // gcc 8.4 has trouble linking
  const int v = 0;
  const int h = std::min(forcedHorizontalMotion + w, dbSize - 1);
  const int d = std::max(0, std::min(dbSize - h, std::min(forcedDiagonalMotion, querySize - width - v)));
  //  std::cerr << "v:" << v << " h:" << h << " d:" << d << std::endl;
  buildWavefronts(v, d, h);

  DRAGEN_SMITH_WATERMAN_LOG << *this << std::endl;
#ifdef TRACE_SMITH_WATERMAN
  printHistoryNice<3>(std::cerr, *this, scores_);
#endif  // TRACE_SMITH_WATERMAN

  return buildCigar(querySize, cigar);
  // TODO: store position, score and cigar into alignment
}

template <typename C, typename T, int WIDTH, int ALIGN, unsigned STEERING_DELAY>
void SmithWatermanT<C, T, WIDTH, ALIGN, STEERING_DELAY>::setQuery(
    const C* queryBegin, const C* queryEnd, bool reverseQuery)
{
  const std::size_t querySize = std::distance(queryBegin, queryEnd);
  query_.resize(querySize + 2 * width - 2);
  std::fill_n(query_.begin(), width - 1, Query::value_type(NOT_A_QUERY_BASE));
  if (!reverseQuery) {
    // since db is revesed for gap left shifting, only reverse queries need to be reversed
    const std::reverse_iterator<const C*> rbegin(queryEnd);
    const std::reverse_iterator<const C*> rend(queryBegin);
    std::copy(rbegin, rend, query_.begin() + width - 1);
  } else {
    // no need to reverse query if it is forward
    std::copy(queryBegin, queryEnd, query_.begin() + width - 1);
  }
  std::fill_n(query_.begin() + querySize + width - 1, width - 1, Query::value_type(NOT_A_QUERY_BASE));
  queryIt_ = query_.begin() + querySize - 1;
}

template <typename C, typename T, int WIDTH, int ALIGN, unsigned STEERING_DELAY>
const char SmithWatermanT<C, T, WIDTH, ALIGN, STEERING_DELAY>::NOT_A_QUERY_BASE = 'b';

template <typename C, typename T, int WIDTH, int ALIGN, unsigned STEERING_DELAY>
const char SmithWatermanT<C, T, WIDTH, ALIGN, STEERING_DELAY>::NOT_A_DB_BASE = 'B';

template <typename C, typename T, int WIDTH, int ALIGN, unsigned STEERING_DELAY>
void SmithWatermanT<C, T, WIDTH, ALIGN, STEERING_DELAY>::setDatabase(
    const C* databaseBegin, const C* databaseEnd)
{
  const std::size_t dbSize = std::distance(databaseBegin, databaseEnd);
  // ensure that getSimilarities does not find matches outside db sequence
  reversedRef_.clear();
  // comparison vector instructions will go outside the dbSize
  // db is always reversed to achieve left shifting of gaps
  //  const std::reverse_iterator<const C*> rbegin(databaseEnd);
  //  const std::reverse_iterator<const C*> rend(databaseBegin);
  // reversedRef_.insert(reversedRef_.end(), rbegin, rend);
  reversedRef_.insert(reversedRef_.end(), databaseBegin, databaseEnd);
  reversedRef_.resize(dbSize + width, NOT_A_DB_BASE);
  reversedRef_.resize(dbSize);
  //  reversedRef_.at(0) = NOT_A_BASE;
  databaseBeginIt_ = reversedRef_.begin();
  databaseEndIt_   = reversedRef_.end();
  dbSize_          = reversedRef_.size();
}

template <typename C, typename T, int WIDTH, int ALIGN, unsigned STEERING_DELAY>
void SmithWatermanT<C, T, WIDTH, ALIGN, STEERING_DELAY>::reset(
    const C* queryBegin, const C* queryEnd, const C* databaseBegin, const C* databaseEnd, bool reverseQuery)
{
  wavefront_.reset();
  setQuery(queryBegin, queryEnd, reverseQuery);
  setDatabase(databaseBegin, databaseEnd);

  scores_.clear();
  const std::size_t querySize = std::distance(queryBegin, queryEnd);
  // this is required to avoid reallocation of score vectors when scores_ runs out of capacity for adding new
  scores_.reserve(querySize * std::distance(databaseBegin, databaseEnd) + width);
  bs_.clear();
  bs_.reserve(querySize * std::distance(databaseBegin, databaseEnd) + width);

  motions_.clear();
  globalMax_       = 0;
  globalMaxScore_  = 0;
  globalMaxOffset_ = -1;
  //  finalMax_ = 0;
  //  finalMaxScore_ = 0;
  //  finalMaxOffset_ = -1;
  queryOffset_ = width - 1;
  // first move right sets tail to the first ref base, which is 0
  dbOffset_ = -width;
  // Apart from avoiding memory reallocations, this one is important as iterators
  // into members of scores_ are kept in maxIndices_.
  // We would not want them to get invalidated when scores_ gets reallocated.
  motions_.reserve(scores_.capacity());

  hyst_shifter_              = 0;
  steer_hyst_                = -HYST_STAGES_C;
  cell_gte_last_next_.first  = 0;
  cell_gte_last_next_.second = 0;
}

template <typename C, typename T, int WIDTH, int ALIGN, unsigned STEERING_DELAY>
void SmithWatermanT<C, T, WIDTH, ALIGN, STEERING_DELAY>::updateHyst(const bool move_ver_v)
{
  const int hyst_add_v = move_ver_v ? 1 : -1;
  const int hyst_sub_v = (hyst_shifter_ & (1 << (HYST_STAGES_C - 1))) ? 1 : -1;
  steer_hyst_          = steer_hyst_ + hyst_add_v - hyst_sub_v;
  hyst_shifter_ <<= 1;
  hyst_shifter_ |= move_ver_v;

  //  std::cerr << "qpos=" << std::max<int>(0, getQueryOffset() - int(WIDTH) + 1) <<
  //  ",rpos=" << (getDatabaseOffset() + int(WIDTH) - 1 + std::min<int>(0, getQueryOffset() - int(WIDTH) + 1)) << ",steer_hyst_=" << int(steer_hyst_) << "," <<
  //    std::bitset<HYST_STAGES_C>(hyst_shifter_) << ",hyst_add=" << hyst_add_v << ",hyst_sub=" << hyst_sub_v << ",move_ver_v=" << move_ver_v << std::endl;
}

template <typename C, typename T, int WIDTH, int ALIGN, unsigned STEERING_DELAY>
void SmithWatermanT<C, T, WIDTH, ALIGN, STEERING_DELAY>::moveDown()
{
  //  std::cerr << "down" << std::endl;
  --queryIt_;
  ++queryOffset_;

  updateHyst(1);
  auto similarities = getSimilarities();
  if (!motions_.empty() && motions_.back() == Motion::down) {
    similarities.front() = 0;
  }

  assert(scores_.capacity() > scores_.size());
  scores_.push_back(wavefront_.moveDown(similarities, gapInit_, gapExtend_));
  bs_.push_back(wavefront_.getLastBs());
  updateMotions(Motion::down);
  updateMax();

  // std::cerr << "qpos=" << std::max<int>(0, getQueryOffset() - int(WIDTH) + 1) <<
  //  ",rpos=" << (getDatabaseOffset() + int(WIDTH) - 1 + std::min<int>(0, getQueryOffset() - int(WIDTH) + 1)) << ",steer_hyst_=" << int(steer_hyst_) << "," <<
  //    std::bitset<HYST_STAGES_C>(hyst_shifter_) << ",hyst_add=" << hyst_add_v << ",hyst_sub=" << hyst_sub_v <<
  //  ",move_ver_v=1, scores:" << scores_.back() << std::endl;
}

template <typename C, typename T, int WIDTH, int ALIGN, unsigned STEERING_DELAY>
void SmithWatermanT<C, T, WIDTH, ALIGN, STEERING_DELAY>::moveRight()
{
  //    std::cerr << "right" << std::endl;
  ++dbOffset_;
  updateHyst(0);

  auto similarities = getSimilarities();
  if (getQueryOffset() >= WIDTH && !motions_.empty() && motions_.back() == Motion::right) {
    // zero to all double right moved ends unless we're the top row
    similarities.back() = 0;
  }

  assert(scores_.capacity() > scores_.size());
  scores_.push_back(wavefront_.moveRight(similarities, gapInit_, gapExtend_));

  bs_.push_back(wavefront_.getLastBs());
  updateMotions(Motion::right);
  updateMax();

  // std::cerr << "qpos=" << std::max<int>(0, getQueryOffset() - int(WIDTH) + 1) <<
  //  ",rpos=" << (getDatabaseOffset() + int(WIDTH) - 1 + std::min<int>(0, getQueryOffset() - int(WIDTH) + 1)) << ",steer_hyst_=" << int(steer_hyst_) << "," <<
  //    std::bitset<HYST_STAGES_C>(hyst_shifter_) << ",hyst_add=" << hyst_add_v << ",hyst_sub=" << hyst_sub_v <<
  //  ",move_ver_v=0, scores:" << scores_.back() << std::endl;
}

template <typename C, typename T, int WIDTH, int ALIGN, unsigned STEERING_DELAY>
void SmithWatermanT<C, T, WIDTH, ALIGN, STEERING_DELAY>::updateMotions(const Motion move)
{
  motions_.push_back(move);
  if (motions_.size() > STEERING_DELAY)  // drop oldest move
  {
    if (move != motions_[motions_.size() - STEERING_DELAY - 1]) {
      drift_ += (move == Motion::right) * 2 - (move == Motion::down) * 2;
    }
  } else  // only accumulate new move
  {
    drift_ += (move == Motion::right) - (move == Motion::down);
  }
}

template <typename C, typename T, int WIDTH, int ALIGN, unsigned STEERING_DELAY>
void SmithWatermanT<C, T, WIDTH, ALIGN, STEERING_DELAY>::updateMax()
{
  auto& s = scores_.back();
  // diagonal can start past the end of the query for short queries
  const int startDeadCells =
      std::max(getQueryOffset() - std::min<int>(getQuerySize() - 1, getQueryOffset()), -getDatabaseOffset());
  // std::cerr << "startDeadCells:" << startDeadCells << std::endl;

  const int dbCovered = std::max(0, getDatabaseSize() - getDatabaseOffset() - startDeadCells);
  // std::cerr << "dbCovered:" << dbCovered << std::endl;

  const auto searchStart = s.begin() + startDeadCells;
  const auto searchEnd =
      s.begin() + std::min<int>(WIDTH, startDeadCells + std::min(dbCovered, getQuerySize()));
  if (searchEnd != searchStart) {
    if (getQueryOffset() >= getQuerySize() - 1) {
      // bonus to bottom row
      s[getQueryOffset() - getQuerySize() + 1] =
          std::max<T>(s[getQueryOffset() - getQuerySize() + 1] + unclipScore_, unclipScore_);
      //    std::cerr << " qo=" << getQueryOffset() << ",do=" << getDatabaseOffset() << ",so=" << (getQueryOffset() - getQuerySize() + 1) << ",bonus=" << scores_.back()[getQueryOffset() - getQuerySize() + 1];
    }

    const auto maxElement = std::max_element(searchStart, searchEnd);
    //  std::cerr << "maxElement:" << maxIndices_.size() << "(" << *maxElement <<  ")," << std::distance(s.begin(), maxElement) << std::endl;

    if (1 == scores_.size() || globalMaxScore_ < *maxElement) {
      globalMaxScore_  = *maxElement;
      globalMax_       = scores_.size() - 1;
      globalMaxOffset_ = std::distance(s.begin(), maxElement);
    }
  }

  //  // offset of the antidiagonal element that contains score of last query base
  //  const int offset = getQueryOffset() - (getQuerySize() - 1);
  //  if (0 <= offset && int(width) > offset)
  //  {
  //    if (finalMaxScore_ <= s.at(offset) + unclipScore_)
  //    {
  //      finalMax_ = scores_.size() - 1;
  //      finalMaxScore_ = s.at(offset) + unclipScore_;
  //      finalMaxOffset_ = offset;
  //    }
  //  }
  //  std::cerr << "getQueryOffset():" << getQueryOffset() << " offset:" << offset  << " finalMax_:" << finalMax_ << " finalMaxScore_:" << finalMaxScore_ << " finalMaxOffset_:" << finalMaxOffset_ << std::endl;
}

template <typename C, typename T, int WIDTH, int ALIGN, unsigned STEERING_DELAY>
std::pair<int, int> SmithWatermanT<C, T, WIDTH, ALIGN, STEERING_DELAY>::getPeakPosition(
    const Antidiagonal& scores, const T minValue) const
{
  auto i     = scores.begin();
  auto first = scores.begin();
  for (; i != scores.end(); ++i) {
    if (*i >= minValue) {
      first = i;
      break;
    }
  }

  auto last = scores.end();
  for (auto i = scores.end() - 1; i != first; --i) {
    if (*i >= minValue) {
      last = i;
      break;
    }
  }
  return std::make_pair((scores.end() - last - 1), (scores.end() - first - 1));
}

template <typename C, typename T, int WIDTH, int ALIGN, unsigned STEERING_DELAY>
typename SmithWatermanT<C, T, WIDTH, ALIGN, STEERING_DELAY>::Motion
SmithWatermanT<C, T, WIDTH, ALIGN, STEERING_DELAY>::getNextMove()
{
  const int index = scores_.size() - 1 - STEERING_DELAY;

  bool steer_ver = false;

  // dragen autosteer takes 4 cycles to determine highest score in a wavefront.
  // By then the wave front is replaced by another 4 times
  // then it takes another cycle to find the bounds of the cells above the threshold
  // and another cycle to compute the steer target.
  // The point is that we need to have a history of all those variables in order to emulate
  // the steering decisions that HW makes
  const int CYCLES_AFTER_PEAK = 4;
  if (0 <= index) {
    static const int ALN_CFG_STEER_DELTA = 12;
    const T          steer_score_v       = *std::max_element(scores_[index].begin(), scores_[index].end());
    assert(CYCLES_AFTER_PEAK <= int(STEERING_DELAY));

    if (!forcedHorizontalMotion_ && !forcedVerticalMotion_ && !forcedDiagonalMotion_ && !autoSteerEnabled_) {
      autoSteerEnabled_ = true;
      // std::cerr << "--------------autosteer---------------" << std::endl;
    }

    //  const int steer_target = steer_lohi.first + steer_lohi.second;
    int steer_target = std::min(cell_gte_last_.first, cell_gte_last_prev_.first) +
                       std::max(cell_gte_last_.second, cell_gte_last_prev_.second);
    //    if (!steer_target)
    //    {
    //      steer_lohi.second = width - 1;
    //      steer_target = steer_lohi.second;
    //    }

    const int steer_baseline = WIDTH + steer_hyst_;
    steer_ver                = steer_target >= steer_baseline;
    const bool move_ver_v    = steer_ver;
    const int  hyst_add_v    = move_ver_v ? 1 : -1;
    const int  hyst_sub_v    = (hyst_shifter_ & (1 << (HYST_STAGES_C - 1))) ? 1 : -1;

    if (autoSteerEnabled_) {
      DRAGEN_S_W_STEER_LOG << "as:"
                           << "steer_score_v=" << steer_score_v << " " << cell_gte_last_.first << ","
                           << cell_gte_last_.second << " "
                           << (scores_.size() % 2 ? cell_gte_last_prev_.first : cell_gte_last_.first) << ","
                           << (scores_.size() % 2 ? cell_gte_last_prev_.second : cell_gte_last_.second) <<
          //        " " << (cell_gte_last_prev_.first) << "," << (cell_gte_last_prev_.second) <<
          ",steer_target=" << steer_target << ",steer_hyst=" << steer_hyst_
                           << ",steer_baseline=" << steer_baseline << ",steer_ver=" << steer_ver << ","
                           << std::bitset<7>(hyst_shifter_) << ",hyst_add=" << hyst_add_v
                           << ",hyst_sub=" << hyst_sub_v << ",move_ver_v=" << steer_ver << ",index=" << index
                           << ",scores_.size()=" << scores_.size()
                           << ",scores=" << scores_[index + CYCLES_AFTER_PEAK - 1]
                           << ",qpos=" << std::max<int>(0, getQueryOffset() - int(WIDTH) + 1) << ",rpos="
                           << (getDatabaseOffset() + int(WIDTH) - 1 +
                               std::min<int>(0, getQueryOffset() - int(WIDTH) + 1))
                           << std::endl;
    }

    cell_gte_last_prev_ = cell_gte_last_;
    cell_gte_last_      = cell_gte_last_next_;
    cell_gte_last_next_ =
        getPeakPosition(scores_[index + CYCLES_AFTER_PEAK], steer_score_v - ALN_CFG_STEER_DELTA);
  }

  if (forcedHorizontalMotion_) {
    --forcedHorizontalMotion_;
    steer_ver = false;
  } else if (forcedVerticalMotion_) {
    --forcedVerticalMotion_;
    steer_ver = true;
  } else if (forcedDiagonalMotion_) {
    forcedDiagonalMotion_ -= !(scores_.size() % 2);
    // query of db too short to have enough scores. Just go diagonally
    steer_ver = !(scores_.size() % 2) ? false : true;
    //  std::cerr << "--diag--" << std::endl;
  } else if (scores_.size() <= STEERING_DELAY - 1) {
    // query of db too short to have enough scores. Just go diagonally
    steer_ver = scores_.size() % 2 ? false : true;
  }

  return steer_ver ? Motion::down : Motion::right;
}

template <typename C, typename T, int WIDTH, int ALIGN, unsigned STEERING_DELAY>
typename SmithWatermanT<C, T, WIDTH, ALIGN, STEERING_DELAY>::Antidiagonal
SmithWatermanT<C, T, WIDTH, ALIGN, STEERING_DELAY>::getSimilarities() const
{
  Antidiagonal similarities;
  const auto   binOp = [this](char q, char d) { return similarity_(q, d); };
  // NOTE, this will go over the databaseEndIt_, but the
  // reversedRef_ is padded to width with NOT_A_BASE chars to achieve consistent results
  //  assert(
  //      std::distance(databaseBeginIt_, databaseIt_) <
  //      (std::distance(databaseBeginIt_, databaseEndIt_) + width));
  //  const int dbCovered = std::min<int>(width, std::distance(databaseBeginIt_, databaseIt_) + 1);
  const int dbCovered = std::min<int>(width, getDatabaseOffset() + width);
  const int qryOffset = width - dbCovered;
  std::fill(similarities.begin(), similarities.end(), similarity_(1, 2));
  //  std::transform(queryIt_ + qryOffset, queryIt_ + width, databaseIt_ - dbCovered + 1, similarities.begin() + qryOffset, binOp);
  std::transform(
      queryIt_ + qryOffset,
      queryIt_ + width,
      databaseBeginIt_ + (getDatabaseOffset() + width - dbCovered),
      similarities.begin() + qryOffset,
      binOp);
  //  std::cerr << "qo:" << getQueryOffset() << "dbCovered:" << dbCovered << ",qryOffset:" << qryOffset << "similarities:" << similarities << std::endl;

  if (getQueryOffset() < WIDTH) {
    // bonus to top row
    similarities[getQueryOffset()] += unclipScore_;
  }

  if (getDatabaseOffset() <= 0) {
    similarities[-getDatabaseOffset()] = 0;
  }

  //  std::cerr << std::endl;
  return similarities;
}

template class SmithWatermanT<unsigned char, short, 48, 16, 9>;
template class SmithWatermanT<unsigned char, short, 48, 16, 11>;

// unit test instances
template class SmithWatermanT<char, short, 48, 16, 9>;
template class SmithWatermanT<char, short, 48, 16, 10>;
template class SmithWatermanT<char, short, 48, 16, 11>;
template class SmithWatermanT<char, short, 8, 16, 10>;
template class SmithWatermanT<char, short, 8, 16, 4>;
}  // namespace align
}  // namespace dragenos
