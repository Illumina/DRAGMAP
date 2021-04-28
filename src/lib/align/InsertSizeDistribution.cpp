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

#include "align/InsertSizeDistribution.hpp"

namespace dragenos {
namespace align {

InsertSizeDistribution::InsertSizeDistribution(
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
    std::ostream& logStream)
  : samplingEnabled_(samplingEnabled),
    dragenInsertStats_(
        0,                        // read group index
        false,                    // true for RNA, false for DNA
        peStatsContinuousUpdate,  // whether to update stats continuously
        peStatsUpdateLogOnly,     // whether updates should just be logged, or also applied to input reads
        peStatsIntervalSize,      // how many reads to accumulate per interval
        peStatsSampleSize,        // how many inserts to include in each calculation
        peStatsIntervalMemory,    // max # of intervals to include in each calculation
        peStatsIntervalDelay,     // lag between calculating and using stats
        std::min(20u, alignerMapqMax),  // minimum MAPQ required for inclusion in stats
        alignerPeOrientation,           // expected orientation of the pairs
        alignerResqueSigmas,            // num std deviations spread to allow rescue alignments
        alignerResqueCeilFactor,        // multiple of read len to compare to rescue sigmas
        alignerResqueMinIns,            // override of rescue min insert length
        alignerResqueMaxIns,            // override of rescue max insert length
        logStream),
    fixedStats_(
        alignerPeMeanInsert,
        alignerPeStddevInsert,
        alignerPeq25,
        alignerPeq50,
        alignerPeq75,
        alignerPeOrientation,
        alignerResqueSigmas,
        alignerResqueCeilFactor,
        alignerResqueMinIns,
        alignerResqueMaxIns,
        alignerPeMeanReadLen,
        false)  // rna mode
{
  if (!samplingEnabled_) {
    std::cerr << "Automatic detection of paired-end insert stats is disabled." << std::endl;
  }
}

InsertSizeParameters InsertSizeDistribution::getInsertSizeParameters(std::size_t r1ReadLen)
{
  InputDbamRecord idr1(r1ReadLen);
  InputDbamRecord idr2;
  if (samplingEnabled_) {
    dragenInsertStats_.waitForValidInterval();
    dragenInsertStats_.saveForRemapping(idr1);
    // this one turns out to collect the total number of bases in order to compute
    // the average read length. This process completely ignores read 2. Don't know why
    dragenInsertStats_.fillPeInsertStats(idr1);
    // There is a tricky interplay between saveForRemapping and fillPeInsertStats
    // saveForRemapping sets state to WAITING, but fillPeInsertStats fails assertion
    // if state is WAITING. So, a dummy R2 record is needed to trigger WAITING state
    // without calling fillPeInsertStats
    dragenInsertStats_.saveForRemapping(idr2);
  } else {
    fixedStats_.fillPeInsertStats(idr1);
  }
  const InsertSizeParameters ret(
      idr1.getInsertStats()->peMinInsert,
      idr1.getInsertStats()->peMaxInsert,
      idr1.getInsertStats()->peMeanInsert,
      idr1.getInsertStats()->rescueMinInsert,
      idr1.getInsertStats()->rescueMaxInsert,
      idr1.getInsertStats()->insertSigmaFactor,
      static_cast<InsertSizeParameters::Orientation>(idr1.getInsertStats()->peOrientation));

  return ret;
}

void InsertSizeDistribution::add(const align::SerializedAlignment& aln)
{
  if (samplingEnabled_) {
    const DbamHeader dbh(aln);
    dragenInsertStats_.sample(&dbh);
  }
}

void InsertSizeDistribution::forceInitDoneSending()
{
  dragenInsertStats_.setInitDoneSending();
  dragenInsertStats_.checkForInitComplete();
}

}  // namespace align
}  // namespace dragenos
