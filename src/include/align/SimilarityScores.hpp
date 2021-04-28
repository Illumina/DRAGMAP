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

#ifndef ALIGN_SIMILARITY_SCORES_HPP
#define ALIGN_SIMILARITY_SCORES_HPP

#include "align/Wavefront.hpp"

namespace dragenos {
namespace align {

/**
 ** \brief similarity scores used in the Smith-Waterman algorithm
 **
 **/
class SimilarityScores {
public:
  SimilarityScores(const short match, const short mismatch, const short nScore = -1)
    : match_(match), mismatch_(mismatch), nScore_(nScore)
  {
  }
  short operator()(char queryBase, char databaseBase) const
  {
    if ((0xF == queryBase) || (0xF == databaseBase) || (0 == queryBase) || (0 == databaseBase)) {
      return nScore_;
    }
    return match_ * (queryBase == databaseBase) + mismatch_ * (queryBase != databaseBase);
  }

  short getSnpCost() const { return match_ - mismatch_; }

  const short match_;
  const short mismatch_;
  const short nScore_;
};

}  // namespace align
}  // namespace dragenos

#endif  // #ifndef ALIGN_SIMILARITY_SCORES_HPP
