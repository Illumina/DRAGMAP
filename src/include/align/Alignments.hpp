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

#ifndef ALIGN_ALIGNMENTS_HPP
#define ALIGN_ALIGNMENTS_HPP

#include "align/Alignment.hpp"
#include "align/SimilarityScores.hpp"
#include "align/SmithWaterman.hpp"
#include "align/VectorSmithWaterman.hpp"

namespace dragenos {
namespace align {

/**
 ** \brief container for alignments
 **
 ** The container exposes a range of valid alignments via the begin() and end()
 ** methods. Otherwise, management of the underlying collection of Alignment
 ** instances is hidden to enable localized caching.
 **/
class Alignments : public std::vector<Alignment> {
public:
  Alignments()
  {
    // TODO: fix the resizing issue in pair builder
    reserve(100000);
  }
  typedef std::vector<Alignment>::iterator       iterator;
  typedef std::vector<Alignment>::const_iterator const_iterator;
  Alignment&                                     addAlignment()
  {
    resize(size() + 1);
    return back();
  }

  void append(const Alignment& a) { push_back(a); }
};

}  // namespace align
}  // namespace dragenos

#endif  // #ifndef ALIGN_ALIGNMENTS_HPP
