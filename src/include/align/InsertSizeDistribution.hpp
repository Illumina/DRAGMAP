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

#ifndef ALIGN_INSERT_SIZE_DISTRIBUTION_HPP
#define ALIGN_INSERT_SIZE_DISTRIBUTION_HPP

#include "align/InsertSizeParameters.hpp"

#include "host/dragen_api/sampling/readgroup_insert_stats.hpp"

namespace dragenos {
namespace align {

class InsertSizeDistribution {
  bool                 samplingEnabled_;
  ReadGroupInsertStats dragenInsertStats_;
  StatsInterval        fixedStats_;

public:
  InsertSizeDistribution(
      bool          samplingEnabled,
      int           alignerPeq25,
      int           alignerPeq50,
      int           alignerPeq75,
      double        alignerPeMeanInsert,
      double        alignerPeStddevInsert,
      double        alignerPeMeanReadLen,
      uint64_t      peStatsIntervalSize,
      uint64_t      peStatsSampleSize,
      uint8_t       peStatsIntervalMemory,
      uint8_t       peStatsIntervalDelay,
      bool          peStatsContinuousUpdate,
      bool          peStatsUpdateLogOnly,
      uint32_t      alignerMapqMax,
      uint32_t      alignerPeOrientation,
      double        alignerResqueSigmas,
      double        alignerResqueCeilFactor,
      uint32_t      alignerResqueMinIns,
      uint32_t      alignerResqueMaxIns,
      std::ostream& logStream);
  //  const InsertSizeParameters& getInsertSizeParameters() const { return insertSizeParameters_; }
  InsertSizeParameters getInsertSizeParameters(std::size_t r1ReadLen);
  void add(const align::SerializedAlignment& aln, const dragenos::sequences::SerializedRead& read);
  bool notGoingToBlock() { return !dragenInsertStats_.justSentAllInitRecords(); }
  void forceInitDoneSending();

  friend std::ostream& operator<<(std::ostream& os, InsertSizeDistribution& d)
  {
    if (d.samplingEnabled_) {
      d.dragenInsertStats_.finalLog();
    } else {
      os << "Paired-end configuration parameters";
      d.fixedStats_.printDescriptiveLog(os);
    }
    return os;
  }
};

}  // namespace align
}  // namespace dragenos

#endif  // #ifndef ALIGN_INSERT_SIZE_DISTRIBUTION_HPP
