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

#ifndef ALIGN_SMITH_WATERMAN_HPP
#define ALIGN_SMITH_WATERMAN_HPP

#include <algorithm>
#include <array>
#include <iomanip>
#include <vector>

#include "align/Alignment.hpp"
#include "align/Antidiagonal.hpp"
#include "align/Database.hpp"
#include "align/Query.hpp"
#include "align/SimilarityScores.hpp"
#include "align/Wavefront.hpp"

namespace dragenos {
namespace align {

/**
 ** \brief a stripped Smith-Waterman implementation with dynamic steering of the wavefront
 **
 ** Only a fixed number of cells are calculated on each antidiagonal. These
 ** cells are alled the Wavefront. The wavefront position for antidiagonal i+1
 ** is either directly right or directly down of wavefront for antidiagonal i.
 ** The direction is predetermined and forced at the beginning and at the end
 ** of the calculation. Otherwise it is dynamically calculated using a simple
 ** rule: keeping the "peak" score in the center of the wavefront.
 **
 ** Note: the steering is based on the estimate of the maximum several cycles
 ** earlier to mimic the implementation on FPGA.
 **
 ** Note: the wavefront is oriented from botom-left to top-right  - this is
 ** the same direction as the database and the direction opposite to the query.
 **/
template <typename C, typename T, int WIDTH, int ALIGN, unsigned STEERING_DELAY>
class SmithWatermanT {
public:
  static const int                       SW_CELLS = WIDTH;
  typedef AntidiagonalT<T, WIDTH, ALIGN> Antidiagonal;
  static const int                       width = Antidiagonal::width;
  static const char                      NOT_A_QUERY_BASE;
  static const char                      NOT_A_DB_BASE;
  typedef std::vector<char>              Database;
  typedef WavefrontT<T, WIDTH, ALIGN>    Wavefront;
  typedef typename Wavefront::Motion     Motion;
  typedef typename Wavefront::Backstep   Backstep;
  typedef typename Wavefront::Backsteps  Backsteps;
  // time required on the FPGA to resolve the steering.
  static constexpr short MAX_RANGE = 3;
  /**
   ** \brief Constructor, registering the appropriate scores
   **/
  SmithWatermanT(
      const SimilarityScores& similarity, const int gapInit, const int gapExtend, const int unclipScore = 0)
    : similarity_(similarity),
      gapInit_(gapInit),
      gapExtend_(gapExtend),
      unclipScore_(unclipScore),
      globalMax_(0),
      cell_gte_last_next_({std::numeric_limits<int>::max(), std::numeric_limits<int>::min()}),
      //    finalMax_(0), finalMaxScore_(0), finalMaxOffset_(-1),
      databaseBeginIt_(0),
      databaseEndIt_(0),
      drift_(0)
  {
  }
  /**
   ** \brief aligns the query against the database
   **
   ** \param query the read to align
   ** \param database the reference to use, with the database expected to start
   ** with "width+2" bases before the expected starting point of the alignment
   ** if any (padding with "not-a-base" if this would go beyond the beginning
   ** of the reference sequence).
   ** \param output parameter to store the alignment information (cigar, score
   ** and position relative to the beginning  of the database)
   **
   **/
  T align(
      const C*     queryBegin,
      const C*     queryEnd,
      const C*     databaseBegin,
      const C*     databaseEnd,
      const int    forcedDiagonalMotion,
      const int    forcedHorizontalMotion,
      bool         reverseQuery,
      std::string& cigar);
  void reset(
      const C* queryBegin,
      const C* queryEnd,
      const C* databaseBegin,
      const C* databaseEnd,
      bool     forwardQuery);
  Antidiagonal getSimilarities() const;
  // identify the position of the "peak" of the scores
  std::pair<int, int> getPeakPosition(const Antidiagonal& scores, const T minValue) const;
  /// decide between right or down based on the peak position and the current relative offset
  Motion getNextMove();
  void   moveRight();
  void   moveDown();
  void   updateMax();
  void   updateMotions(Motion move);
  T      getMaxScore() const { return globalMaxScore_; };

private:
  const SimilarityScores    similarity_;
  const T                   gapInit_;
  const T                   gapExtend_;
  const T                   unclipScore_;
  Wavefront                 wavefront_;
  std::vector<Antidiagonal> scores_;
  std::vector<Backsteps>    bs_;
  // position of the max score in each antidiagonal
  int globalMaxOffset_ = -1;
  // index of the antidiagonal containing the global max
  size_t globalMax_      = -1;
  T      globalMaxScore_ = 0;

  //  // index of the antidiagonal containing the maximum score of last query base
  //  size_t finalMax_;
  //  T finalMaxScore_;
  //  int finalMaxOffset_;
  static const int STAGE_AUTO_STEER_VLD_C = STEERING_DELAY;
  static const int HYST_STAGES_C          = 7;

  static const int SHIFTER_MASK_ = (1 << HYST_STAGES_C) - 1;

  int                 hyst_shifter_ = 0;               //SHIFTER_MASK_;
  int                 steer_hyst_   = -HYST_STAGES_C;  //0;
  std::pair<int, int> cell_gte_last_next_;
  std::pair<int, int> cell_gte_last_;
  std::pair<int, int> cell_gte_last_prev_;
  // keeping track of motions
  std::vector<typename Wavefront::Motion> motions_;
  /**
   ** \brief copy of the query in reverse memory order and padding
   **
   ** To calculate the similarities, the query and the database must be
   ** accessed in opposite directions, which prevents the compiler from
   ** using SSE instructions to load the sequence read against memory
   ** order. Also, for the initialization and finalization of alignment,
   ** some form of "not-a-base" padding is needed either on the query or
   ** on the database to calculate the scores in the top-left and bottom-
   ** right half-corners.
   ** To resolve both concerns without adding tests in the inner loops of
   ** the implementation, the idea is to make a copy of the query with the
   ** elements stored in reverse memory order and with suddicient padding
   ** at both ends;
   **
   ** Note: the SSE 4.2 instruction to load a string (ant mem-alignment)
   ** into an SIMD register is faster than shifting the register and
   ** inserting a value at the end.
   **
   ** Extra note: Since reference needs always to be reversed in order to
   ** obtain left-shifted gaps, query may or may not be reversed from
   ** original input depending on reverseQuery parameter of the align method
   **/
  std::vector<C> query_;
  std::vector<C> reversedRef_;
  /**
   ** \brief forward iterator on the reverseQuery_
   **
   ** Use "widh" positions starting at queryIt_ to fill the Antidiagonal.
   ** Decrease when moving the wavefront down.
   **/
  typename std::vector<C>::const_iterator queryIt_;
  // offset of the bottom left element of antidiagonal in query sequence
  int queryOffset_ = -1;
  /**
   ** \brief forward iterator on the database
   **
   ** Use "width" positions starting at databasEIt_ to fill the Antidiagonal.
   ** Increase when moving right.
   **/
  typename std::vector<C>::const_iterator databaseBeginIt_;
  // offset of bottom left element of antidiagonal in db sequence
  int dbOffset_ = 0;
  //  const C *databaseBeginIt_;
  //  const C *databaseIt_;
  /// End of the database
  typename std::vector<C>::const_iterator databaseEndIt_;
  int                                     dbSize_ = 0;
  //  const C *databaseEndIt_;
  /// drift accumulated over the last STEERING_DELAY cycles (+1/-1 for right/down resp)
  ssize_t drift_;

  int  forcedHorizontalMotion_ = 0;
  int  forcedVerticalMotion_   = 0;
  int  forcedDiagonalMotion_   = 0;
  bool autoSteerEnabled_       = false;

  void setQuery(const C* queryBegin, const C* queryEnd, bool leftShiftIndels);
  // Note: the database is not copied and must remain valid (including same iterators) while being used
  void setDatabase(const C* databaseBegin, const C* databaseEnd);
  void setDatabaseSize(int dbSize)
  {
    assert(dbSize_ >= dbSize);
    dbSize_ = dbSize;
  }
  void buildWavefronts(
      const size_t forcedVerticalMotion,
      const size_t forcedDiagonalMotion,
      const size_t forcedHorizontalMotion);
  T buildCigar(const std::size_t querySize, std::string& cigar);
  T buildCigar(const std::size_t querySize, int antidiagonal, int maxOffset, std::string& ret);

  int getQuerySize() const { return query_.size() - 2 * width + 2; }
  // offset of the bottom left element of antidiagonal in query sequence
  int getQueryOffset() const { return queryOffset_; }

  int getDatabaseSize() const { return dbSize_; }
  // offset of bottom left element of antidiagonal in db sequence
  int getDatabaseOffset() const { return dbOffset_; }

  void updateHyst(const bool move_ver_v);

public:
  auto             getScores() const -> const decltype(scores_)& { return scores_; }
  auto             getGlobalMax() const -> decltype(globalMax_) { return globalMax_; }
  auto             getMotions() const -> const decltype(motions_)& { return motions_; }
  auto             getDrift() const -> decltype(drift_) { return drift_; }
  auto             getReversedQuery() const -> const decltype(query_)& { return query_; }
  auto             getQueryIt() const -> decltype(queryIt_) { return queryIt_; }
  const Wavefront& getWavefront() const { return wavefront_; }

  //  template <typename P, typename V> P getPrintable(V v)
  //  {
  //    return v;
  //  }

  static T getPrintable(T v) { return v; }

  static const char* getPrintable(Backstep bs)
  {
    static const char* BS[] = {
        "x",    "<",     "^",    "<^",   "\\",    "\\<", "\\^", "\\<^", "-x",   "-<",    "-^",
        "-<^",  "-\\",   "-\\<", "-\\^", "-\\<^", "|x",  "|<",  "|^",   "|<^",  "|\\",   "|\\<",
        "|\\^", "|\\<^", "+x",   "+<",   "+^",    "+<^", "+\\", "+\\<", "+\\^", "+\\<^",
    };

    assert(size_t(bs) < sizeof(BS) / sizeof(BS[0]));

    return BS[bs];
  }

  template <int JUSTIFY, typename HistoryT>
  static void printHistoryNice(std::ostream& os, const SmithWatermanT& sw, const HistoryT& history)
  {
    //    static char ACGT[256] = {'A', 'C', 'G', 'T'};
    static char ACGT[256]       = {'N', 'A', 'C', '?', 'G', '?', '?', '?', 'T'};
    ACGT[int('A')]              = 'A';
    ACGT[int('C')]              = 'C';
    ACGT[int('G')]              = 'G';
    ACGT[int('T')]              = 'T';
    ACGT[int('a')]              = 'A';
    ACGT[int('c')]              = 'C';
    ACGT[int('g')]              = 'G';
    ACGT[int('t')]              = 'T';
    ACGT[int('n')]              = 'N';
    ACGT[int('N')]              = 'N';
    ACGT[int(NOT_A_QUERY_BASE)] = 'x';
    ACGT[int(NOT_A_DB_BASE)]    = 'X';

    os << std::setfill(' ') << std::setw(JUSTIFY) << " ";
    for (auto dbit = sw.databaseBeginIt_; sw.databaseEndIt_ != dbit; ++dbit) {
      os << std::setw(JUSTIFY) << ACGT[int(*dbit)];
    }
    os << std::endl;

    for (std::size_t q = 0; q < sw.query_.size() - 2 * width + 2; ++q) {
      os << std::setw(JUSTIFY) << ACGT[int(*(sw.query_.rbegin() + (q + width - 1)))];
      int offset = -q - 1 + width;
      for (std::size_t i = 0; history.size() > i; ++i) {
        offset += (sw.motions_.at(i) == Motion::down);
        if (q <= i) {
          if (offset >= 0 && offset < int(Antidiagonal::width)) {
            os << std::setw(JUSTIFY) << getPrintable(*(history.at(i).begin() + offset));
          } else if (offset < 0) {
            os << std::setw(JUSTIFY) << " ";
          }
        }
      }
      os << std::endl;
    }
    os << std::endl;
  };

  template <typename HistoryT>
  static void printHistory(std::ostream& os, const SmithWatermanT& sw, const HistoryT& history)
  {
    int qpos = width - 1;
    int rpos = -width;

    for (std::size_t i = 0; history.size() > i; ++i) {
      qpos += (sw.motions_.at(i) == Motion::down);
      rpos += (sw.motions_.at(i) == Motion::right);
      os << "qpos=" << std::max(0, qpos - int(sw.width) + 1)
         << ", rpos=" << rpos + int(sw.width) - 1 + std::min(0, qpos - int(sw.width) + 1) << ", scores=";
      const auto& h = history.at(i);
      for (int offset = h.size() - 1; 0 <= offset; --offset) {
        if ((qpos - offset) < 0) {
        } else if ((rpos + offset) >= sw.getDatabaseSize()) {
          os << "-,";
        } else if ((qpos - offset) >= sw.getQuerySize()) {
          os << "+,";
        } else {
          os << h.at(offset) << ",";
        }
      }
      os << std::endl;
    }
  };

  friend std::ostream& operator<<(std::ostream& os, const SmithWatermanT& sw)
  {
    assert(sw.motions_.size() == sw.scores_.size());
    //    for (const auto& s : sw.scores_)
    //    {
    //      os << s << std::endl;
    //    }
    //
    //    for (const auto& m : sw.motions_)
    //    {
    //      os << (Motion::down == m ? "down" : "right") << std::endl;
    //    }

    //    printHistory(os, sw, sw.e_);
    //    printHistory(os, sw, sw.f_);
    printHistory(os, sw, sw.scores_);
    //    printHistoryNice<3>(std::cerr, sw, sw.scores_);
    //    printHistoryNice<5>(os, sw, sw.bs_);

    return os;
  }
};

typedef SmithWatermanT<unsigned char, short, 48, 16, 9> SmithWaterman;

template <typename C, typename T, int WIDTH, int ALIGN, unsigned STEERING_DELAY>
const int SmithWatermanT<C, T, WIDTH, ALIGN, STEERING_DELAY>::width;

}  // namespace align
}  // namespace dragenos

#endif  // #ifndef ALIGN_SMITH_WATERMAN_HPP
