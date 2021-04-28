// Copyright 2016 Edico Genome Corporation. All rights reserved.
//
// This file contains confidential and proprietary information of the Edico Genome
// Corporation and is protected under the U.S. and international copyright and other
// intellectual property laws.
//
// $Id$
// $Author$
// $Change$
// $DateTime$
//

#include "stats_interval.hpp"
#include "dragen_run_log.hpp"
#include "input_dbam_record.hpp"
#include "kernel_density.hpp"
#include "output_dbam_header.hpp"
#include "run_stats.hpp"

//--------------------------------------------------------------------------------adamb
// Instantiation/initializations of static constant members:
const int      StatsInterval::DEFAULT_STDDEV            = 10000;
const float    StatsInterval::OUTLIER_BOUND             = 2.0;
const double   StatsInterval::MAPPING_BOUND             = 3.0;
const double   StatsInterval::MAX_STDDEV                = 4.0;
const double   StatsInterval::DEFAULT_RESCUE_SIGMAS     = 2.5;
const uint32_t StatsInterval::DEFAULT_RESCUE_MIN_INSERT = 1;
const uint32_t StatsInterval::DEFAULT_RESCUE_MAX_INSERT = 1000;

//--------------------------------------------------------------------------------adamb
// constructor
//
StatsInterval::StatsInterval(
    const uint16_t rgId,     // which readgroup
    const uint8_t  id,       // which stats interval, within that read group
    const bool     rnaMode,  // true for RNA, false for DNA
    const uint32_t peOrientation,
    const double   rescueSigmas,      // num std deviations spread to allow rescue alignments
    const double   rescueCeilFactor,  // multiple of read len to compare to rescue sigmas
    const uint32_t rescueMinInsert,   // override of rescue min insert length
    const uint32_t rescueMaxInsert)   // override of rescue max insert length
  : m_rgId(rgId),
    m_id(id),
    m_enableRescueAlignments(true),
    m_rnaMode(rnaMode),
    m_peOrientation(peOrientation),
    m_valid(false),
    m_numSent(0),
    m_doneSending(false),
    m_numReceived(0),
    m_q25(-1),
    m_q50(0),
    m_q75(0),
    m_s50(-1),
    m_low(0),
    m_high(0),
    m_numInsertsInMean(0),
    m_mean(0),
    m_stddev(0),
    m_minInsert(1),
    m_maxInsert(-1),
    m_insertSigmaFactor(0),
    m_rescueSigmasOverridden(rescueSigmas != 0),
    m_rescueSigmas(m_rescueSigmasOverridden ? rescueSigmas : DEFAULT_RESCUE_SIGMAS),
    m_rescueCeilFactor(rescueCeilFactor),
    m_rescueRadius(0),
    m_rescueMinOverridden(rescueMinInsert != 0),
    m_rescueMinInsert(m_rescueMinOverridden ? rescueMinInsert : DEFAULT_RESCUE_MIN_INSERT),
    m_rescueMaxOverridden(rescueMaxInsert != 0),
    m_rescueMaxInsert(m_rescueMaxOverridden ? rescueMaxInsert : DEFAULT_RESCUE_MAX_INSERT),
    m_warningMode(WARN_OK)
{
  if (m_rnaMode) m_rescueMaxInsert = 0;
}

//--------------------------------------------------------------------------------adamb
// Version of the constructor to use when forcing stats based on config file
//
StatsInterval::StatsInterval(
    const double   mean,              // mean value for insert-size
    const double   stddev,            // std-deviation for insert-size
    const uint32_t q25,               // first quartile
    const uint32_t q50,               // second quartile
    const uint32_t q75,               // third quartile
    const uint32_t peOrientation,     // expected orientations of two ends of pairs
    const double   rescueSigmas,      // num std deviations spread to allow rescue alignments
    const double   rescueCeilFactor,  // multiple of read len to compare to rescue sigmas
    const uint32_t rescueMinInsert,   // override of rescue min insert length
    const uint32_t rescueMaxInsert,   // override of rescue max insert length
    const float    meanReadLen,       // average bases per read
    const bool     rnaMode)               // true for RNA, false for DNA
  : m_rgId(255),
    m_id(255),
    m_enableRescueAlignments(true),
    m_rnaMode(rnaMode),
    m_peOrientation(peOrientation),
    m_valid(true),
    m_numSent(0),
    m_doneSending(false),
    m_numReceived(0),
    m_q25(q25),
    m_q50(q50),
    m_q75(q75),
    m_s50(q75 - q25),
    m_low(0),
    m_high(0),
    m_numInsertsInMean(0),
    m_mean(mean),
    m_stddev(stddev),
    m_minInsert(getMinInsert()),
    m_maxInsert(getMaxInsert()),
    m_insertSigmaFactor(getSigmaFactor()),
    m_rescueSigmasOverridden(rescueSigmas != 0),
    m_rescueSigmas(m_rescueSigmasOverridden ? rescueSigmas : DEFAULT_RESCUE_SIGMAS),
    m_rescueCeilFactor(rescueCeilFactor),
    m_rescueRadius(0),
    m_rescueMinOverridden(rescueMinInsert != 0),
    m_rescueMinInsert(m_rescueMinOverridden ? rescueMinInsert : 0),
    m_rescueMaxOverridden(rescueMaxInsert != 0),
    m_rescueMaxInsert((m_rescueMaxOverridden && !m_rnaMode) ? rescueMaxInsert : 0),
    m_warningMode(WARN_OK)
{
  updateOutlierFilter();
  updateRescueInsertLimits(meanReadLen);
}

//--------------------------------------------------------------------------------adamb
// Copy constructor
//
StatsInterval::StatsInterval(const StatsInterval& other)
  : m_rgId(other.m_rgId),
    m_id(other.m_id),
    m_enableRescueAlignments(other.m_enableRescueAlignments),
    m_rnaMode(other.m_rnaMode),
    m_peOrientation(other.m_peOrientation),
    m_valid(other.m_valid),
    m_numSent(other.m_numSent),
    m_doneSending(other.m_doneSending),
    m_numReceived(other.m_numReceived),
    m_warningMode(other.m_warningMode)
{
  copyStats(other);
}

//--------------------------------------------------------------------------------adamb
// when we start re-filling an interval, first thing we need to do is set it back to
// an "invalid" state and clear its stats.
//
void StatsInterval::reset()
{
  m_valid             = false;
  m_numSent           = 0;
  m_doneSending       = false;
  m_numReceived       = 0;
  m_q25               = -1;
  m_q50               = 0;
  m_q75               = 0;
  m_s50               = -1;
  m_low               = 0;
  m_high              = 0;
  m_numInsertsInMean  = 0;
  m_mean              = 0;
  m_stddev            = 0;
  m_minInsert         = 1;
  m_maxInsert         = -1;
  m_insertSigmaFactor = 0;
  //  m_rescueCeilFactor = 1.0; don't change rescue ceil factor from constructor value
  m_rescueRadius = 0;

  if (!m_rescueSigmasOverridden) {
    m_rescueSigmas = DEFAULT_RESCUE_SIGMAS;
  }

  if (!m_rescueMinOverridden) {
    m_rescueMinInsert = DEFAULT_RESCUE_MIN_INSERT;
  }

  if (m_rnaMode)
    m_rescueMaxInsert = 0;
  else if (!m_rescueMaxOverridden) {
    m_rescueMaxInsert = DEFAULT_RESCUE_MAX_INSERT;
  }
  m_warningMode = WARN_OK;
}

//--------------------------------------------------------------------------------adamb
// Copy the stats results from a different interval
//
void StatsInterval::copyStats(const StatsInterval& other)
{
  m_valid                  = other.m_valid;
  m_q25                    = other.m_q25;
  m_q50                    = other.m_q50;
  m_q75                    = other.m_q75;
  m_s50                    = other.m_s50;
  m_low                    = other.m_low;
  m_high                   = other.m_high;
  m_numInsertsInMean       = other.m_numInsertsInMean;
  m_mean                   = other.m_mean;
  m_stddev                 = other.m_stddev;
  m_minInsert              = other.m_minInsert;
  m_maxInsert              = other.m_maxInsert;
  m_insertSigmaFactor      = other.m_insertSigmaFactor;
  m_rescueSigmasOverridden = other.m_rescueSigmasOverridden;
  m_rescueSigmas           = other.m_rescueSigmas;
  m_rescueCeilFactor       = other.m_rescueCeilFactor;
  m_rescueRadius           = other.m_rescueRadius;
  m_rescueMinOverridden    = other.m_rescueMinOverridden;
  m_rescueMinInsert        = other.m_rescueMinInsert;
  m_rescueMaxOverridden    = other.m_rescueMaxOverridden;
  m_rescueMaxInsert        = other.m_rescueMaxInsert;
}

//--------------------------------------------------------------------------------swhitmore/adamb
// For an interval that doesn't include any inserts that fall between the outlier bounds,
// we need to just set default values.  In version 1.0, this led to an error message about
// there being no high-quzality pairs
//
void StatsInterval::setDefaultStats()
{
  m_mean      = 0;
  m_q25       = 0;
  m_q50       = 0;
  m_q75       = 0;
  m_low       = 0;
  m_high      = 0;
  m_stddev    = DEFAULT_STDDEV;
  m_minInsert = getMinInsert();
  m_maxInsert = getMaxInsert();

  if (!m_rescueMinOverridden) {
    m_rescueMinInsert = DEFAULT_RESCUE_MIN_INSERT;
  }
  if (m_rnaMode)
    m_rescueMaxInsert = 0;
  else if (!m_rescueMaxOverridden) {
    m_rescueMaxInsert = DEFAULT_RESCUE_MAX_INSERT;
  }
  m_warningMode = WARN_NONE;
}

//--------------------------------------------------------------------------------adamb
// Output current stats to the log file
//
void StatsInterval::log(std::ostream& os) const
{
  os << m_rgId << "\t";
  os << (int)m_id << "\t";
  os << m_q25 << "\t" << m_q50 << "\t" << m_q75 << "\t";
  os << m_mean << "\t";
  os << m_stddev << "\t";
  os << m_minInsert << "\t";
  os << m_maxInsert << "\t";
  os << m_rescueMinInsert << "\t";
  os << m_rescueMaxInsert << "\t";
  os << m_numInsertsInMean << "\t";
  os << m_inserts.size() << "\t";
}

//--------------------------------------------------------------------------------theoh
// update final insertSize stats to Runstats ( one element per read group )

void StatsInterval::updateRunStats()
{
  RunStats::RGinsertSizeStats rgStats;

  switch (m_warningMode) {
  case (WARN_NONE):
    rgStats.m_errMsg = "WARNING: limited reads to estimate insert stats. Use defaults.";
    break;
  case (WARN_THREE):
    rgStats.m_errMsg = "WARNING: limited reads to estimate insert stats. Use default standard deviation.";
    break;
  case (WARN_TWENTYEIGHT):
    rgStats.m_errMsg = "WARNING: limited reads to estimate insert stats. Use small samples formula.";
    break;
  default:
    rgStats.m_errMsg = "";
    break;
  }

  rgStats.m_numInsertsInMean = m_numInsertsInMean;
  rgStats.m_q25              = m_q25;
  rgStats.m_q50              = m_q50;
  rgStats.m_q75              = m_q75;
  rgStats.m_mean             = m_mean;
  rgStats.m_stddev           = m_stddev;
  rgStats.m_low              = m_low;
  rgStats.m_high             = m_high;
  rgStats.m_minInsert        = m_minInsert;
  rgStats.m_maxInsert        = m_maxInsert;
  RunStats::Instance()->addInsertSizeRG(rgStats);
}

//--------------------------------------------------------------------------------adamb
// Output stats in the "final" format, rather than the "iterative" format used
// in log().
//
void StatsInterval::printDescriptiveLog(std::ostream& os) const
{
  os << ", based on " << m_numInsertsInMean << " high quality pairs for FR orientation\n";
  os << "        Quartiles (25 50 75) = " << m_q25 << " " << m_q50 << " " << m_q75 << "\n";
  os << "        Mean = " << m_mean << "\n";
  os << "        Standard deviation = " << m_stddev << "\n";

  // Indicate to the user the rescue radius and rescue sigmas that was used, since the host software
  // may override the default based on read length. For RNA mode, rescue scans are disabled, so we never
  // print it for RNA mode.
  if (!m_rnaMode) {
    const double effectiveRescueSigmas = m_rescueRadius / m_stddev;
    os << "        Rescue radius = " << m_rescueRadius << "\n";
    os << "        Effective rescue sigmas = " << effectiveRescueSigmas << "\n";
    if (effectiveRescueSigmas < m_rescueSigmas) {
      os << "        WARNING: Default rescue sigmas value of " << DEFAULT_RESCUE_SIGMAS
         << " was overridden by host software!\n";
      os << "        The user may wish to set a rescue sigmas value explicitly with --Aligner.rescue-sigmas\n";
    }
  }

  os << "        Boundaries for mean and standard deviation: low = " << m_low << ", high = " << m_high
     << "\n";
  os << "        Boundaries for proper pairs: low = " << m_minInsert << ", high = " << m_maxInsert << "\n";
  os << "        NOTE: DRAGEN's insert estimates include corrections for clipping "
     << "(so they are not identical to TLEN)\n";
}

//--------------------------------------------------------------------------------adamb
// If there weren't enough samples to get good stats, issue a warning.
//
void StatsInterval::warnIfTooLittleData(std::ostream& os) const
{
  switch (m_warningMode) {
  case (WARN_NONE):
    os << "WARNING: No high quality pairs found for calculating paired-end "
       << "statistics -- used defaults\n";
    break;
  case (WARN_THREE):
    os << "WARNING: Less than 3 high quality pairs found - standard deviation "
       << "set to " << DEFAULT_STDDEV << "\n";

    LOG_RUN("Number of inserts in mean is less than 3, stddev was ", m_stddev);
    break;
  case (WARN_TWENTYEIGHT):
    os << "WARNING: Less than 28 high quality pairs found - standard deviation "
       << "is calculated from the small samples formula\n";
    LOG_RUN("Number of inserts in mean is less than 28, stddev was ", m_stddev);
    break;
  default:
    break;
  }
}

//--------------------------------------------------------------------------------adamb
// Print detailed info to help diagnose an error condition
//
void StatsInterval::logError() const
{
  std::cerr << "  Num sent: " << m_numSent << std::endl;
  std::cerr << "  Num received: " << m_numReceived << std::endl;
  std::cerr << "  Done sending? " << (m_doneSending ? "yes" : "no") << std::endl;
}

//--------------------------------------------------------------------------------adamb
// Compute a low and high threshold after setting the quartiles, and prior to computing
// mean and standard deviation. Insert sizes not within the low-high threshold will
// be discarded from mean and standard deviation computation
//
void StatsInterval::updateOutlierFilter()
{
  m_low = static_cast<int>(m_q25 - (OUTLIER_BOUND * m_s50) + .499);
  if (m_low < 1) {
    m_low = 1;
  }
  m_high = static_cast<int>(m_q75 + (OUTLIER_BOUND * m_s50) + .499);
  if (m_high < 1) {
    m_high = std::numeric_limits<int>::max();
  }
}

//--------------------------------------------------------------------------------adamb
// update rescue-{min,max}-insert
//
void StatsInterval::updateRescueInsertLimits(const float meanReadLen)
{
  if (m_rnaMode) {
    m_rescueMinInsert = 0;
    m_rescueMaxInsert = 0;
  } else {
    if (m_rescueMinOverridden && m_rescueMaxOverridden) return;

    // Rescues alignments should be run plus-or-minus some number of standard
    // deviations from the mean read length. This is specified by the rescue-sigmas
    // parameter.
    if (m_rescueSigmasOverridden) {
      m_rescueRadius = m_rescueSigmas * m_stddev;
    } else {
      // Under default conditions, the host software further limits the rescue scan
      // radius as a function of the read length.
      const double rescueFloorFactor = 1.0;
      m_rescueRadius                 = std::max(rescueFloorFactor * meanReadLen, m_rescueSigmas * m_stddev);
      m_rescueRadius                 = std::min(m_rescueCeilFactor * meanReadLen, m_rescueRadius);
    }

    if (m_rescueRadius / m_stddev < 0.5 || !m_enableRescueAlignments /*|| m_rnaMode*/) {
      // But not for RNA mode, or if there is no variance in insert size, or if
      // rescues are entirely turned off.
      if (!m_rescueMinOverridden) m_rescueMinInsert = 0;
      if (!m_rescueMaxOverridden) m_rescueMaxInsert = 0;

      // If rescue scans are disabled, then update m_rescueRadius to indicate it as such.
      if (!m_rescueSigmasOverridden and m_rescueMinInsert == 0 and m_rescueMaxInsert == 0) {
        m_rescueRadius = 0;
      }
    } else {
      if (!m_rescueMinOverridden) {
        m_rescueMinInsert = std::max(round(m_mean - m_rescueRadius), static_cast<double>(m_minInsert));
      }

      if (!m_rescueMaxOverridden) {
        m_rescueMaxInsert = std::min(round(m_mean + m_rescueRadius), static_cast<double>(m_maxInsert));
      }
    }
  }
}

//-------------------------------------------------------------------------------swhitmore
// calculate - Once sampling is finished, calculate paired-end insert statistics
// from insert sizes of high quality pairs.
//
void StatsInterval::calculate(Inserts_c& inserts, const double meanReadLen)
{
  if (inserts.empty()) {
    setDefaultStats();
    m_valid = true;
    return;
  }

  std::sort(inserts.begin(), inserts.end());

  // For RNA, perform kernel density estimation to alter the inserts list
  if (m_rnaMode and inserts.size() > 100) {
    calculateStatsRna(inserts);
  }

  // Compute the quartiles
  const size_t len = inserts.size();
  if (len > 1) {
    m_q25 = inserts[static_cast<int>((.25 * len) + .499)];
    m_q50 = inserts[static_cast<int>((.50 * len) + .499)];
    m_q75 = inserts[static_cast<int>((.75 * len) + .499)];
  } else {
    // There is only one insert in the list
    m_q25 = inserts[0];
    m_q50 = m_q25;
    m_q75 = m_q25;
  }

  // Compute the center 50 span
  m_s50 = m_q75 - m_q25;

  updateOutlierFilter();

  // Compute the mean
  for (size_t i = 0; i < inserts.size(); ++i) {
    if ((inserts[i] >= (unsigned)m_low) && inserts[i] <= (unsigned)m_high) {
      m_mean += inserts[i];
      ++m_numInsertsInMean;
    } else {
      //      ++m_numInsertsFiltered;
      if (inserts[i] > (unsigned)m_high) {
        // We have a sorted list - all other elements in the list are going to be larger
        break;
      }
    }
  }

  if (!m_numInsertsInMean) {
    // No non-outlier stats.  We are done.
    setDefaultStats();
    m_valid = true;
    return;
  }

  m_mean = m_mean / m_numInsertsInMean;

  // Compute the standard deviation
  for (size_t i = 0; i < inserts.size(); ++i) {
    if ((inserts[i] >= (unsigned)m_low) && inserts[i] <= (unsigned)m_high) {
      m_stddev += (inserts[i] - m_mean) * (inserts[i] - m_mean);
    }
  }
  m_stddev = sqrt(m_stddev / m_numInsertsInMean);

  // Handle case with a small number of samples (and stddev is not reliable)
  if (m_numInsertsInMean < 3) {
    m_stddev = DEFAULT_STDDEV;
  } else if (m_numInsertsInMean < 28) {
    m_stddev = (25 * (m_stddev + 1)) / (m_numInsertsInMean - 2);
  }
  m_stddev = std::max(12.0, m_stddev);

  // Compute the min and max insert for proper pairs
  m_minInsert         = getMinInsert();
  m_maxInsert         = getMaxInsert();
  m_insertSigmaFactor = getSigmaFactor();

  updateRescueInsertLimits(meanReadLen);

  m_valid = true;
}

//--------------------------------------------------------------------------------adamb
// Copy this interval's stats into the destination mapper input record
//
void StatsInterval::fillPeInsertStats(InputDbamRecord& dbr)
{
  ASSERT(
      m_valid,
      "m_id: ",
      (int)m_id,
      ", m_doneSending: ",
      (m_doneSending ? "true" : "false"),
      ", m_numSent: ",
      m_numSent,
      ", m_numReceived: ",
      m_numReceived);
  DbamInsertStats* stats = dbr.getInsertStats();
  if (!stats) return;

  stats->peMeanInsert      = m_mean;
  stats->peMinInsert       = m_minInsert;
  stats->peMaxInsert       = m_maxInsert;
  stats->rescueMinInsert   = m_rescueMinInsert;
  stats->rescueMaxInsert   = m_rescueMaxInsert;
  stats->insertSigmaFactor = m_insertSigmaFactor;
  stats->unused0           = 0;
  stats->unused1[0] = stats->unused1[1] = stats->unused1[2] = 0;
  stats->peOrientation                                      = m_peOrientation;
}

void print_intlist(StatsInterval::Inserts_c& ints)
{
  auto it = ints.begin();
  while (true) {
    for (uint32_t i = 0; i < 10; ++i) {
      if (it == ints.end()) return;
      std::cerr << (*it) << "\n";
      ++it;
    }
  }
}

//--------------------------------------------------------------------------------adamb
// Return the most common value in the list of inserts
//
uint32_t StatsInterval::getHistogramPeak(const Inserts_c& inserts)
{
  uint32_t peak    = 0;
  size_t   peak_n  = 0;
  size_t   n       = 0;
  uint32_t lastVal = 0;

  for (auto it = inserts.begin(); it != inserts.end(); ++it) {
    if ((*it) == lastVal) {
      // We're in the middle of a string of identical values.  Increment the count.
      ++n;
    } else {
      // We've reached the end of a string of identical values.  Check if we had
      // more of them than the previous winner:
      if (n > peak_n) {
        peak_n = n;
        peak   = (*it);
      }

      // Move on to the new value, and restart counting from zero
      lastVal = (*it);
      n       = 0;
    }
  }

  return peak;
}

//-------------------------------------------------------------------------------seant
// CalculateStatsRna - Function to adjust the inserts list for RNA processing by
// using a kernel density estimator to smooth out the samples and determine the peak
// value.
//
void StatsInterval::calculateStatsRna(
    Inserts_c& inserts)  // the insert sizes to use for calculating the distribution
{
  // This function performs the following steps to get a modified sample set:
  //  1. Use the kernel density estimator ONLY to find the smoothed-peak position (insert size).
  //     The estimator is used in a 2-pass method. The first pass to find the peak P1. Then in
  //     in the second pass we keep only samples that are < 2*P1 and run the estimator again.
  //     The peak value of this iteration P is used going forward.
  //     This should be ≥ 1, considering all the samples were positive integers.
  //  2. Round the peak to an integer value P ≥ 1.
  //  3. Discard insert samples > P.
  //  4. For every insert sample X < P, add an artificial sample P + (P – X).
  //  5. Then proceed with legacy sample processing.  All samples are positive, between 1 and 2P–1, and the
  //  mean and median are P.
  //

  //  std::cerr << "calculateStatsRna input:\n";
  //  print_intlist(inserts);

  // Kernel Density Estimator for first pass
  KernelDensity* kd = new KernelDensity();

  // Add the samples
  for (size_t i = 0; i < inserts.size(); ++i) {
    kd->AddSample(static_cast<float>(inserts[i]));
  }

  // Calculate the density and determine the peak
  double density_max_x = 0;
  double density_max_y = 0;

  double   min_x          = kd->GetMin(0);
  uint32_t histogram_peak = getHistogramPeak(inserts);
  double   max_x          = 3.0 * histogram_peak;

  // Bail if there turned out to be no variation in the stats
  const double TOL = 1e-3;
  if (fabs(min_x - max_x) < TOL) {
    std::cerr << "No variation in RNA stats." << std::endl;
    delete kd;
    return;
  }

  double x_increment = (max_x - min_x) / 1000.0;
  for (double x = min_x; x <= max_x; x += x_increment) {
    double y = kd->Pdf(x);
    if (y > density_max_y) {
      density_max_y = y;
      density_max_x = x;
    }
  }

  // Finished with first pass
  delete kd;

  // Second pass
  kd = new KernelDensity();

  // Add the samples, but only samples that are < 2*peak from first pass
  uint32_t threshold = 2 * static_cast<uint32_t>(density_max_x + .499);
  for (size_t i = 0; i < inserts.size(); ++i) {
    // The list is already sorted, so we can break the loop once we pass the threshold
    if (inserts[i] >= threshold) {
      break;
    }
    kd->AddSample(static_cast<float>(inserts[i]));
  }

  // Calculate the density and determine the peak
  density_max_x = 0;
  density_max_y = 0;

  min_x = kd->GetMin(0);
  max_x = kd->GetMax(0);
  // Bail if there turned out to be no variation in the stats
  if (fabs(min_x - max_x) < TOL) {
    std::cerr << "No variation in RNA stats." << std::endl;
    delete kd;
    return;
  }

  x_increment = (max_x - min_x) / 1000.0;
  for (double x = min_x; x <= max_x; x += x_increment) {
    double y = kd->Pdf(x);
    if (y > density_max_y) {
      density_max_y = y;
      density_max_x = x;
    }
  }

  uint32_t peak_insert = static_cast<uint32_t>(density_max_x + .499);
  //std::cerr << "RNA Kernel Density Peak: " << peak_insert << "\n";

  // Modify the list of inserts, discarding samples > P, and for every insert sample X
  // that is < P, add an artificial sample X' = P + (P - X)
  Inserts_c modifiedInserts;
  modifiedInserts.reserve(inserts.size());

  int64_t idx = 0;
  for (; idx < static_cast<int64_t>(inserts.size()); ++idx) {
    if (inserts[idx] > peak_insert) {
      break;
    }
    modifiedInserts.push_back(inserts[idx]);
  }

  // Proceed to add the reflected portion
  for (; idx >= 0; --idx) {
    if (inserts[idx] >= peak_insert) {
      continue;
    }
    modifiedInserts.push_back(2 * peak_insert - inserts[idx]);
  }

  // Update the inserts list
  inserts.swap(modifiedInserts);

  //  std::cerr << "calculateStatsRna output:\n";
  //  print_intlist(inserts);
  //  std::sort(inserts.begin(), inserts.end());
}

//--------------------------------------------------------------------------------adamb
// In response to SIGUSR1, dump information on the state of this object to the
// hang-dump log file
//
void StatsInterval::DumpHangInfo(std::ostream& os)
{
  os << "      RG " << m_rgId << ", interval " << (int)m_id << "\n";
  os << "        Valid? " << (m_valid ? "yes" : "no") << "\n";
  os << "        Num sent: " << m_numSent << "\n";
  os << "        Done sending? " << (m_doneSending ? "yes" : "no") << "\n";
  os << "        Num received: " << m_numReceived << std::endl;
}
