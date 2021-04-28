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

#include "map/ChainBuilder.hpp"
#include "common/DragenLogger.hpp"

namespace dragenos {
namespace map {

void ChainBuilder::addSeedPosition(
    const SeedPosition& seedPosition, const bool reverseComplement, const bool randomSample)
{
///////////
#ifdef DEBUG_LOG
  std::cerr << "  ChainBuilder::addSeedPosition: seedPosition: " << std::hex << "0x" << std::uppercase
            << std::setw(9) << std::setfill('0') << seedPosition.getReferencePosition() << std::dec
            << std::endl;
#endif
  ///////////

  bool accepted = false;
  for (auto& seedChain : *this) {
    if (seedChain.accepts(seedPosition, reverseComplement)) {
///////////
#ifdef DEBUG_LOG
      std::cerr << "  Accepted by: " << seedChain << std::endl;
#endif
      ///////////
      seedChain.addSeedPosition(seedPosition, randomSample);
      accepted = true;
      // TODO: verify that we don't have to add the seed possition to all accepting seed chains
      //return;
    }
/////////
#ifdef DEBUG_LOG
    else {
      std::cerr << "  Rejected by: " << seedChain << std::endl;
    }
#endif
    /////////
  }
  if (accepted) {
    return;
  }
  // none of the existing chains accept the seed at this position. add new chain with this seed
/////////
#ifdef DEBUG_LOG
  std::cerr << "  None Accepted, formed its own chain " << std::endl;
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
  /*
#ifdef DEBUG_LOG
  std::cerr << "\n------------------------ ChainBuilder::filterChains()"<<std::endl;
#endif

  SeedChain* superior = nullptr;
  SeedChain* inferior = nullptr;

  for (int i = 0; i < seedChainCount_; ++i)
  {
    if(seedChains_[i].isFiltered()) continue;
    for (int j = i+1; j < seedChainCount_; ++j)
    {
      if(seedChains_[j].isFiltered()) continue;
      if(seedChains_[i].getReadSpanLength() > seedChains_[j].getReadSpanLength())
      {
        superior = &seedChains_[i];
        inferior = &seedChains_[j];
      }
      else
      {
        superior = &seedChains_[j];
        inferior = &seedChains_[i];
      }

      double filterThresh = chainFilterRatio_ * static_cast<float>(inferior->getReadSpanLength()) + chainFilterConstant_;
#ifdef DEBUG_LOG
      std::cerr << "superior:"<<*superior<<"superior->getReadSpanLength():"<<superior->getReadSpanLength()<<"\n";
      std::cerr << "inferior:"<<*inferior<<"inferior->getReadSpanLength():"<<inferior->getReadSpanLength()<<std::endl;
      std::cerr << "superior->size():"<<superior->size()<<"\n";
      std::cerr << "chainFilterRatio_:"<<chainFilterRatio_<<"\n";
      std::cerr << "inferior->getReadSpanLength():"<<inferior->getReadSpanLength()<<"\n";
      std::cerr << "chainFilterConstant_:"<<chainFilterConstant_<<"\n";
      std::cerr << "(chainFilterRatio_ * (inferior->getReadSpanLength())) + chainFilterConstant_):"<<filterThresh<<std::endl;
#endif
      if (superior->size() >= filterThresh) {
        double overlapThresh = inferior->getReadSpanLength() * 0.25;
        if(superior->firstReadBase() - overlapThresh <= inferior->firstReadBase()
           and superior->lastReadBase() + overlapThresh >= inferior->lastReadBase())
          inferior->setFiltered(true);
      }
    }
  }
*/
  for (const SeedChain& chain : *this) {
    DRAGEN_CHAIN_LOG << chain << std::endl;
  }
}

}  // namespace map
}  // namespace dragenos
