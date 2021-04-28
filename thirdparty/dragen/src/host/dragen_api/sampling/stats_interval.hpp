// Copyright 2016-2017 Edico Genome Corporation. All rights reserved.
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

#ifndef __STATS_INTERVAL_HPP__
#define __STATS_INTERVAL_HPP__

#include <inttypes.h>
#include <algorithm>
#include <cmath>
#include <ostream>
#include <vector>

class DbamHeader;
class InputDbamRecord;

//--------------------------------------------------------------------------------adamb
// The paired-end stats system has a phased system for applying stats.  When we
// send a given read, we map it using insert stats calculated based on (about)
// actual map/align results from 100,000 reads ago.  Then we that same read comes
// out of the mapper/aligner, we apply its insert-size to reads that will come
// (about) 100,000 reads later.
//
// To make this system deterministic, we divide the input reads into intervals.
// The read sender generates map/align inputs using stats from a fixed number of
// intervals in the past -- we call that the "interval delay". The reads are tagged
// with a number in the input record's FBLOB to direct it to the correct interval
// for tallying the output.  When we detect that all of the reads destined for an
// output interval have made it through the mapper/aligner, we gather up a group
// of insert sizes and calculate the updated stats.
//
// This class takes care of aggregating the insert-sizes and calculating the stats
// for a single interval.  The work of coordinating and timing which reads go to
// which intervals is handled in the ReadGroupInsertStats class.
//
class StatsInterval {
  // Let unit tests inspect the internals of this class:
  friend class StatsIntervalTest;

public:
  typedef std::vector<uint32_t>     Inserts_c;
  typedef Inserts_c::const_iterator Inserts_it;

private:
  static const int      DEFAULT_STDDEV;             // 10000
  static const float    OUTLIER_BOUND;              // 2.0
  static const double   MAPPING_BOUND;              // 3.0
  static const double   MAX_STDDEV;                 // 4.0
  static const double   DEFAULT_RESCUE_SIGMAS;      // 2.5
  static const uint32_t DEFAULT_RESCUE_MIN_INSERT;  // 1
  static const uint32_t DEFAULT_RESCUE_MAX_INSERT;  // 1000

private:
  enum WarningMode { WARN_OK, WARN_NONE, WARN_THREE, WARN_TWENTYEIGHT };

public:
  StatsInterval(
      const uint16_t rgId,     // read group id
      const uint8_t  id,       // index of interval within that read group
      const bool     rnaMode,  // true for RNA, false for DNA
      const uint32_t peOrientation,
      const double   rescueSigmas,      // num std deviations spread to allow rescue alignments
      const double   rescueCeilFactor,  // multiple of read len to compare to rescue sigmas
      const uint32_t rescueMinInsert,   // override of rescue min insert length
      const uint32_t rescueMaxInsert)   // override of rescue max insert length
      ;

  StatsInterval(                        // constructor to use when forcing stats based on config file
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
      const float    meanReadLen,       // average read length
      const bool     rnaMode)               // true for RNA, false for DNA
      ;

  StatsInterval(const StatsInterval& other);

  void copyStats(const StatsInterval& other);

  void reset();

  bool receivedAll() const { return (m_doneSending && (m_numReceived == m_numSent)); }

  bool isValid() const { return m_valid; }

  void addInsert(const uint32_t insertSize) {
    m_inserts.push_back(insertSize);
  }

  void clearInserts()
  {
    // Completely clear out the inserts list, and any of its capacity:
    Inserts_c emptyInserts;
    m_inserts.swap(emptyInserts);
  }

  void calculate(                // set stats for this interval
      Inserts_c&   inserts,      // the insert sizes to use for calculating the distribution
      const double meanReadLen)  // average length of the reads in this sample
      ;

  void fillPeInsertStats(InputDbamRecord& dbr);

  void setDoneSending() { m_doneSending = true; }

  uint64_t getNumSent() const { return m_numSent; }

  void incNumSent() { ++m_numSent; }

  void incNumReceived() { ++m_numReceived; }

  uint64_t getNumReceived() const { return m_numReceived; }

  void log(std::ostream& os) const;

  void updateRunStats();

  void printDescriptiveLog(std::ostream& os) const;

  void logError() const;

  void warnIfTooLittleData(std::ostream& os) const;

  // Return iterators over the pair-insert-lengths:
  std::pair<Inserts_it, Inserts_it> getInserts() const
  {
    return std::make_pair(m_inserts.begin(), m_inserts.end());
  }

  void enableRescueAlignments(const bool enableit) { m_enableRescueAlignments = enableit; }

  void DumpHangInfo(std::ostream& os);  // Dump state in response to SIGUSR1

private:
  void setDefaultStats();

  void updateOutlierFilter();

  void updateRescueInsertLimits(const float meanReadLen);

  // Calculate the minimum insert size for proper pairing
  uint32_t getMinInsert() const
  {
    return std::round(
        std::max(1.0, std::min(m_q25 - (m_s50 * MAPPING_BOUND), m_mean - (MAX_STDDEV * m_stddev))));
  }

  // Calculate the maximum insert size for proper pairing
  uint32_t getMaxInsert() const
  {
    return std::min(
        65535.0, round(std::max(m_q75 + (m_s50 * MAPPING_BOUND), m_mean + (MAX_STDDEV * m_stddev))));
  }

  uint16_t getSigmaFactor() const
  {
    return std::min(static_cast<uint16_t>(0xFFFF), static_cast<uint16_t>(round(0x2F200 / m_stddev)));
  }

  void calculateStatsRna(Inserts_c& inserts)  // the insert sizes to use for calculating the distribution
      ;

  uint32_t getHistogramPeak(const Inserts_c& inserts);

private:
  // General configuration that stays the same the whole run:
  const uint16_t m_rgId;                    // read group id.
  const uint8_t  m_id;                      // index of this interval within its read group
  bool           m_enableRescueAlignments;  // if false, force off all rescue alignments
  const bool     m_rnaMode;                 // true for RNA, false for DNA
  const uint32_t m_peOrientation;           // FF, FR, FF, RR from config file

  // Management of whether this interval is still filling, or is ready for just
  // handing out results
  bool     m_valid;        // whether the stats have been calculated and are valid
  uint64_t m_numSent;      // number of reads sent to the mapper/aligner on input side
  bool     m_doneSending;  // whether we have sent all of the reads into the mapper/aligner
  uint64_t m_numReceived;  // number of reads that we've received from map/align on the output side.

  Inserts_c m_inserts;

  // Statistics from sampling:
  int         m_q25;               // first insert size quartile
  int         m_q50;               // second insert size quartile
  int         m_q75;               // third insert size quartile
  int         m_s50;               // The center 50 span
  int         m_low;               // Low boundary for computing mean and standard deviation
  int         m_high;              // High boundary for computing mean and standard deviation
  int         m_numInsertsInMean;  // How many inserts were included in calculating the mean
  double      m_mean;              // The mean insert size
  double      m_stddev;            // One standard deviation of the insert
  uint32_t    m_minInsert;         // Minimum insert size for proper pairing
  uint32_t    m_maxInsert;         // Maximum insert size for proper pairing
  uint16_t    m_insertSigmaFactor;
  bool        m_rescueSigmasOverridden;  // whether rescueSigmas was set by the user
  double      m_rescueSigmas;
  double      m_rescueCeilFactor;
  double      m_rescueRadius;
  bool        m_rescueMinOverridden;  // whether rescueMinInsert was manually overridden
  uint32_t    m_rescueMinInsert;
  bool        m_rescueMaxOverridden;  // whether rescueMinInsert was manually overridden
  uint32_t    m_rescueMaxInsert;
  WarningMode m_warningMode;
};

#endif
