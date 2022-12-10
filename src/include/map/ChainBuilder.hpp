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

#ifndef MAP_CHAIN_BUILDER_HPP
#define MAP_CHAIN_BUILDER_HPP

#include <algorithm>
#include <cstdint>
#include <utility>
#include <vector>

#include "map/SeedChain.hpp"
#include "map/SeedPosition.hpp"

namespace dragenos {
namespace map {

/**
 ** \brief a component dedicated to constructing chains of related mapping seed positions
 **
 ** A chain is a sequence defined by a pair of iterator (begin, end]. Each element of
 ** the chain is the pair (seed, position).
 **
 ** This component serves two purposes, first encapsulating the logic of chain creation
 ** into an independent class, and second enabling the client code to reuse allocated
 ** internal buffers for positions and chains.
 **
 **/
class ChainBuilder {
public:
  ChainBuilder(double chainFilterRatio) : seedChainCount_(0), chainFilterRatio_(chainFilterRatio) {}
  void                                   clear() { seedChainCount_ = 0; }
  std::vector<SeedChain>::const_iterator begin() const { return seedChains_.begin(); }
  std::vector<SeedChain>::const_iterator end() const { return seedChains_.begin() + seedChainCount_; }

  std::vector<SeedChain>::iterator begin() { return seedChains_.begin(); }
  std::vector<SeedChain>::iterator end() { return seedChains_.begin() + seedChainCount_; }

  SeedChain&       at(std::size_t i) { return seedChains_.at(i); }
  const SeedChain& at(std::size_t i) const { return seedChains_.at(i); }
  const SeedChain& back() const { return seedChains_.at(seedChainCount_ - 1); }

  void   reserve(std::size_t max) { seedChains_.reserve(max); }
  size_t size() const { return seedChainCount_; }
  void   addSeedPosition(const SeedPosition& seedPosition, bool reverseComplement, bool randomSample);
  void   addSeedChain(const SeedChain& seedChain);
  void   filterChains();
  template <class F>
  void sort(F compare)
  {
    std::sort(seedChains_.begin(), seedChains_.begin() + seedChainCount_, compare);
  }

  friend std::ostream& operator<<(std::ostream& os, const ChainBuilder& chains)
  {
    for (const SeedChain& chain : chains) {
      os << "\n" << chain;
    }
    return os;
  }

  void setFilterConstant(int a) { chainFilterConstant_ = 1 - a; }

private:
  std::vector<SeedChain> seedChains_;
  size_t                 seedChainCount_;
  double                 chainFilterRatio_    = 2.0;
  double                 chainFilterConstant_ = 0.0;
  // int                    numRandomSampleHits_ = 0;
};

}  // namespace map
}  // namespace dragenos

#endif  // #ifndef MAP_CHAIN_BUILDER_HPP
