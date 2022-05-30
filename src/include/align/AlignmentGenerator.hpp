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

#ifndef ALIGN_ALIGNMENT_GENERATOR_HPP
#define ALIGN_ALIGNMENT_GENERATOR_HPP

#include "align/Alignments.hpp"
#include "align/SmithWaterman.hpp"
#include "align/VectorSmithWaterman.hpp"
#include "map/ChainBuilder.hpp"
#include "reference/ReferenceDir.hpp"
#include "sequences/Read.hpp"

namespace dragenos {
namespace align {

class AlignmentGenerator {
public:
  typedef sequences::Read   Read;
  typedef map::ChainBuilder ChainBuilder;
  AlignmentGenerator(
      const reference::ReferenceSequence& refSeq,
      const reference::HashtableConfig&   htConfig,
      SmithWaterman&                      smithWaterman,
      VectorSmithWaterman&                vectorSmithWaterman,
      bool                                vectorizedSW)
    : refSeq_(refSeq),
      htConfig_(htConfig),
      smithWaterman_(smithWaterman),
      vectorSmithWaterman_(vectorSmithWaterman),
      vectorizedSW_(vectorizedSW)
  {
  }
  /// delegate for the Aligner generateAlignments method
  void generateAlignments(
      const Read& read, const map::ChainBuilder& chainBuilder, Alignments& alignments, const int readIdx);
  bool generateAlignment(
      const ScoreType alnMinScore,
      const Read&     read,
      map::SeedChain  seedChain,
      Alignment&      alignment,
      const int       readIdx);

private:
  const reference::ReferenceSequence& refSeq_;
  const reference::HashtableConfig&   htConfig_;
  SmithWaterman&                      smithWaterman_;
  VectorSmithWaterman&                vectorSmithWaterman_;
  const bool                          vectorizedSW_;
  void updateFetchChain(const Read& read, map::SeedChain& seedChain, Alignment& alignment);
};  // class AlignmentGenerator

}  // namespace align
}  // namespace dragenos

#endif  // #ifndef ALIGN_ALIGNMENT_GENERATOR_HPP
