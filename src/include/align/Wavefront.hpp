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

#ifndef ALIGN_WAVEFRONT_HPP
#define ALIGN_WAVEFRONT_HPP

#include <array>

#include "align/Antidiagonal.hpp"

namespace dragenos {
namespace align {

/**
 ** \brief enables incremental calculation of the Smith-Waterman scoring matrix
 **
 ** The wavefront is simply a partial antidiagonal in the Smith-Waterman
 ** scoring matrix. The assumption being that the gaps are smaller than the
 ** width of the wavefront and that the scores on the antidiagonal outside the
 ** wavefront are all null. This assumption enables calculating the scoring
 ** matrices for longer reads without the prohibitive quatratic cost.
 **
 ** The component provides the two elementary methods to move the wavefront
 ** right and down. Each of these methods is further split into the detail
 ** calculation of the three underlying components used to keep track of the
 ** operations related to opening and extending insertions and deletions.
 **
 ** The responsibility of this component is strictly restricted to the short
 ** term management of the wavefrons and only keeps track of the information
 ** strictly necessary to move the wavefront one step right or one step down.
 ** It is the responsibility of the client code to manage the steering of the
 ** wavefront and to store the history of scores and travel pathh.
 **
 ** As this is a very low level component, additional functionalities are
 ** provided to allow initialization and inspection of the internals. Also,
 ** the high-level methos moveRight and moveDown are deconstructed to enable
 ** fine grained control of this component in diverse contexts.
 **
 ** Note that the indices in the antidiagonal are from top-right to
 ** botom-left.
 **/
template <typename T = short, int WIDTH = 48, int ALIGN = 16>
class WavefrontT {
  /// size of the history stored internally: penultimate, last and next
  static constexpr size_t SIZE = 3;

public:
  typedef AntidiagonalT<T, WIDTH, ALIGN>    Antidiagonal;
  typedef typename Antidiagonal::value_type Int;
  typedef std::array<Antidiagonal, SIZE>    History;
  enum Motion { right = 0, down = 1 };

  enum Backstep : char {
    none     = 0,
    horz     = 1,
    vert     = 2,
    horzVert = 3,
    diag     = 4,
    horzDiag = 5,
    vertDiag = 6,
    all      = 7,

    extHFlag  = 8,
    horzExtH  = horz | extHFlag,
    vertH     = vert | extHFlag,
    horzVertH = horzVert | extHFlag,
    diagH     = diag | extHFlag,
    horzDiagH = horzDiag | extHFlag,
    vertDiagH = vertDiag | extHFlag,
    allH      = all,

    extVFlag  = 16,
    horzV     = horz | extVFlag,
    vertV     = vert | extVFlag,
    horzVertV = horz | extVFlag,
    diagV     = diag | extVFlag,
    horzDiagV = horzDiag | extVFlag,
    vertDiagV = vertDiag | extVFlag,
    allV      = all | extVFlag,

    flags     = extVFlag | extHFlag,
    horzA     = horz | flags,
    vertA     = vert | flags,
    horzVertA = horz | flags,
    diagA     = diag | flags,
    horzDiagA = horzDiag | flags,
    vertDiagA = vertDiag | flags,
    allVA     = all | flags,

    everything = flags | all
  };
  typedef std::array<Backstep, WIDTH> Backsteps;

  WavefrontT() : e_(), f_(), h_(), next_(0), moved_(right) {}
  WavefrontT(History e, History f, History h, size_t next, Motion moved)
    : e_(e), f_(f), h_(h), next_(next), moved_(moved)
  {
  }
  WavefrontT(const WavefrontT& w) : e_(w.e_), f_(w.f_), h_(w.h_), next_(w.next_), moved_(w.moved_) {}
  WavefrontT& operator=(const WavefrontT w)
  {
    e_     = w.e_;
    f_     = w.f_;
    h_     = w.h_;
    next_  = w.next_;
    moved_ = w.moved_;
    return *this;
  }
  const History& getE() const { return e_; }
  const History& getF() const { return f_; }
  const History& getH() const { return h_; }
  /// reset all scores to 0 to start a new alignment
  void reset();
  /// scores along the last calculated antidiagonal - 0 is top right
  const Antidiagonal& getLastScores() const { return h_[last()]; }
  const Antidiagonal& getLastE() const { return e_[last()]; }
  const Antidiagonal& getLastF() const { return f_[last()]; }
  const Backsteps&    getLastBs() const { return bs_; }
  /**
   ** \brief move the wavefront one cell to the right
   **
   ** \return the calculated scores
   **
   ** \param query has a begin() methods that provides an iterator to the begining
   **        of the query sequence relative to the wavefront
   ** \param query has an end() method that provides an iterator to the end of the
   **        reference database relative to the wavefront
   ** \param similarity the scores provided are consistenten with the encoding of
   **        bases in query and database
   **
   ** Both query and database have at least as many bases as the width of the
   ** wavefront
   **
   ** These assumptions will enable the client code to meaningfully handle all
   ** situations where the scoring matrix calculation goes past the beginning
   ** or past the end of the query or database sequences.
   ** It is the responsibility of the client code to manage the query and the
   ** database accordingly.
   **/
  const Antidiagonal& moveRight(const Antidiagonal& similarities, const Int gapInit, const Int gapExtend);
  /**
   ** \brief move the wavefront one cell downwards
   **
   ** see "moveRight" for parameters and return value
   **/
  const Antidiagonal& moveDown(const Antidiagonal& similarities, const Int gapInit, const Int gapExtend);
  /// deconstruction of moveRight, only for the E component
  const Antidiagonal& moveRightE(const Int gapInit, const Int gapExtend);
  /// deconstruction of moveDown, only for the E component
  const Antidiagonal& moveDownE(const Int gapInit, const Int gapExtend);
  /// deconstruction of moveRight, only for the F component
  const Antidiagonal& moveRightF(const Int gapInit, const Int gapExtend);
  /// deconstruction of moveDown, only for the F component
  const Antidiagonal& moveDownF(const Int gapInit, const Int gapExtend);
  /// deconstruction of moveRight, only for the H component, assuming E anf F have already been calculated
  const Antidiagonal& moveRightH(const Antidiagonal& similarities);
  /// deconstruction of moveDown, only for the F component, assuming E anf F have already been calculated
  const Antidiagonal& moveDownH(const Antidiagonal& similarities);
  /// set the next antidiagonal values to the max between 0, Eij, Fij and Hij and advance the wavefront to the
  /// next position
  const Antidiagonal& setNextToMax();
  size_t              next() const { return next_; }
  size_t              last() const { return (next_ + SIZE - 1) % SIZE; }
  size_t              penultimate() const { return (next_ + SIZE - 2) % SIZE; }

private:
  // tracking deletions - we need only 2 E, but it we need 3 H
  std::array<Antidiagonal, SIZE> e_;
  // tracking insertions - we need only 2 F, but it we need 3 H
  std::array<Antidiagonal, SIZE> f_;
  // actual similarity scores - needs the penultimate for cells positioned diagonally
  std::array<Antidiagonal, SIZE> h_;
  // the next offset to use in e_, f_ and h_

  // directions towards the start of the alignment
  Backsteps bs_;

  size_t next_;
  // needed for the relative position of h_
  Motion moved_;

  void resetBs();
  bool selectBest(const T extend, const T open, T& ret);
};

// for testing
typedef WavefrontT<short, 48, 16> Wavefront48;

}  // namespace align
}  // namespace dragenos

#endif  // #ifndef ALIGN_WAVEFRONT_HPP
