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

#include <assert.h>

#include "common/DragenLogger.hpp"
#include "map/ChainBuilder.hpp"

namespace dragenos {
namespace map {

void ChainBuilder::addSeedPosition(
    const SeedPosition& seedPosition, const bool reverseComplement, const bool randomSample)
{
///////////
#ifdef TRACE_SEED_CHAINS
  std::cerr << "  ChainBuilder::addSeedPosition: seedPosition: " << std::hex << "0x" << std::uppercase
            << std::setw(9) << std::setfill('0') << seedPosition.getReferencePosition() << std::dec
            << std::endl;
#endif
  ///////////

  bool accepted = false;
  for (auto& seedChain : *this) {
    if (seedChain.accepts(seedPosition, reverseComplement)) {
///////////
#ifdef TRACE_SEED_CHAINS
      std::cerr << seedPosition << "  Accepted by: " << seedChain << std::endl;
#endif
      ///////////
      seedChain.addSeedPosition(seedPosition, randomSample);
      accepted = true;
      // TODO: verify that we don't have to add the seed possition to all accepting seed chains
      //return;
    }
/////////
#ifdef TRACE_SEED_CHAINS
    else {
      //      std::cerr << seedPosition << "  Rejected by: " << seedChain << std::endl;
    }
#endif
    /////////
  }
  if (accepted) {
    return;
  }
  // none of the existing chains accept the seed at this position. add new chain with this seed
/////////
#ifdef TRACE_SEED_CHAINS
  std::cerr << seedPosition << "  None Accepted, formed its own chain " << std::endl;
#endif
  ///////////
  if (seedChains_.size() == seedChainCount_) {
    seedChains_.resize(seedChainCount_ + 1);
  }
  seedChains_[seedChainCount_].clear();
  seedChains_[seedChainCount_].setReverseComplement(reverseComplement);
  seedChains_[seedChainCount_].addSeedPosition(seedPosition, randomSample);
  ++seedChainCount_;
}

void ChainBuilder::addSeedChain(const SeedChain& seedChain)
{
  if (seedChains_.size() == seedChainCount_) {
    seedChains_.push_back(seedChain);
  } else {
    seedChains_[seedChainCount_] = seedChain;
  }
  ++seedChainCount_;
}

void ChainBuilder::filterChains()
{
#ifdef TRACE_SEED_CHAINS
  std::cerr << "\n------------------------ ChainBuilder::filterChains()" << std::endl;
#endif

  SeedChain* inferior = nullptr;

  int maxcov_beg   = std::numeric_limits<int>::max();
  int maxcov_end   = std::numeric_limits<int>::min();
  int max_coverage = 0;
  for (int i = 0; i < seedChainCount_; ++i) {
    if (seedChains_[i].getReadCovLength() > max_coverage) {
      max_coverage = seedChains_[i].getReadCovLength();
      maxcov_beg   = seedChains_[i].firstReadBase();
      maxcov_end   = seedChains_[i].lastReadBase();
    }
  }

  assert(!seedChainCount_ || maxcov_beg != std::numeric_limits<int>::max());
  assert(!seedChainCount_ || maxcov_end != std::numeric_limits<int>::min());

  for (int j = 0; j < seedChainCount_; ++j) {
    inferior = &seedChains_[j];

    int chain_filt_ratio = int(chainFilterRatio_ * (1 << 8));
    int cov_thresh       = (chain_filt_ratio * inferior->getReadCovLength()) >> 8;  // + chainFilterConstant_;
    {
      int len_div4_v = inferior->getReadCovLength() >> 2;
      if (maxcov_beg <= inferior->firstReadBase() + len_div4_v &&
          maxcov_end >= inferior->lastReadBase() - len_div4_v && max_coverage >= cov_thresh &&
          !inferior->isExtra() && !inferior->hasOnlyRandomSamples()) {
        inferior->setFiltered(true);
#ifdef TRACE_SEED_CHAINS
        std::cerr <<  // "superior:" << *superior << "\n"
            "inferior:" << *inferior << std::endl;
        std::cerr << "chain_length=" << inferior->getReadCovLength() << ",len_div4_v=" << len_div4_v
                  << ",max_coverage=" << max_coverage << ",cov_thresh=" << cov_thresh
                  << ",maxcov_beg=" << maxcov_beg << ",maxcov_end=" << maxcov_end
                  << ",map_cfg.chain_filt_ratio=" << chain_filt_ratio
                  << " chainFilterConstant_:" << chainFilterConstant_ << std::endl;
#endif
      }
    }
  }

  for (const SeedChain& chain : *this) {
    DRAGEN_CHAIN_LOG << chain << std::endl;
  }
}

}  // namespace map
}  // namespace dragenos
