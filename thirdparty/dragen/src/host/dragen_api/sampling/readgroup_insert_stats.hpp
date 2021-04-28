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

#ifndef __READGROUP_INSERT_STATS_HPP__
#define __READGROUP_INSERT_STATS_HPP__

#include <boost/thread.hpp>
#include <cmath>
#include <deque>
#include <iostream>
#include <vector>

#include "stats_interval.hpp"
#include "input_dbam_record.hpp"
#include "input_dbam_remapper_interface.hpp"
#include "output_dbam_header.hpp"

//--------------------------------------------------------------------------------adamb
// This module contains logic to track and calculate insert-size statistics for a
// single read group.
//
class ReadGroupInsertStats {
  friend class AlignedReadsSamplerTest;

  enum InitializationState { SENDING = 0, WAITING = 1, DONE = 2 };
  enum SenderState { NORMAL = 0, STALLED = 1 };
  enum { INITIAL_HISTOGRAM_SIZE = 1000 };

public:
  ReadGroupInsertStats(
      const uint16_t rgId,              // read group index
      const bool     rnaMode,           // true for RNA, false for DNA
      const bool     continuousUpdate,  // whether to update stats continuously
      const bool     updateLogOnly,  // whether updates should just be logged, or also applied to input reads
      const uint32_t readsPerInterval,   // how many reads to accumulate per interval
      const uint32_t sampleSize,         // how many inserts to include in each calculation
      const uint8_t  maxIntervalMemory,  // max # of intervals to include in each calculation
      const uint32_t intervalDelay,      // lag between calculating and using stats
      const uint32_t mapqMin,            // minimum MAPQ required for inclusion in stats
      const uint32_t peOrientation,      // expected orientation of the pairs
      const double   rescueSigmas,       // num std deviations spread to allow rescue alignments
      const double   rescueCeilFactor,   // multiple of read len to compare to rescue sigmas
      const uint32_t rescueMinInsert,    // override of rescue min insert length
      const uint32_t rescueMaxInsert,    // override of rescue max insert length
      std::ostream&  logStream)           // where to log info on stats
      ;

  ~ReadGroupInsertStats();

public:
  uint16_t getReadGroupId() const { return m_readGroupId; }

  void sample(const DbamHeader* dbh);

  void printConfig(std::ostream& os, const uint8_t interval) const;

  void finalLog();

  void logError() const;

  void DumpHangInfo(std::ostream& os);

  void fillPeInsertStats(InputDbamRecord& dbr);

  void waitForValidInterval();

  void enableRescueAlignments(const bool enableit);

  void setReadSender(InputDbamRemapper* readSender) { m_readSender = readSender; }

  void setDummyRead(char* dummy) { memcpy((void*)m_dummyRead, (void*)dummy, DUMMY_READ_SIZE); }

public:
  bool justSentAllInitRecords()
  {
    const bool sentall = ((m_initState == WAITING) && (m_intervals[0]->getNumSent() == INIT_INTERVAL_SIZE));

    if (!sentall) return false;

    // We will only ever return "true" once, to make sure that the caller is only seeing
    // the state transition.
    if (m_reportedWaiting) return false;

    m_reportedWaiting = true;
    return true;
  }

  void saveForRemapping(InputDbamRecord& dbr);

  //--------------------------------------------------------------------------------adamb
  // Have we sent all of the init records for mapping against default stats?  If
  // so, return the count of reads that need to be remapped with proper stats
  size_t getNumRemapRecords() const
  {
    if (m_initState != DONE) return 0;

    return m_remapRecords.size();
  }

  //--------------------------------------------------------------------------------adamb
  // Return pointers to all of the records that need remapping with non-default stats
  template <class OutputIter>
  OutputIter getRemapRecords(OutputIter o)
  {
    if (m_remapRecords.empty()) return o;

    assert(m_sendingInterval == 0);
    for (auto it = m_remapRecords.begin(); it != m_remapRecords.end(); ++it) {
      *o++ = (*it);
    }

    // Make sure that we only return the remapping reads once:
    m_remapRecords.clear();

    return o;
  }

  void setInitDoneSending();

private:
  void completeInitialization();

  void checkIntervalForCompletion(const uint8_t interval);

public:
  void checkForInitComplete();

  bool isInitDone()
  {
    boost::unique_lock<boost::mutex> lock(m_mutex);
    return (m_initState == DONE);
  }

private:
  void waitUntilInitDone();

  //--------------------------------------------------------------------------------adamb

private:
  void flushUntilIntervalValid(uint32_t useInterval);
  bool shouldUseRead(  // determine whether a record should be included in stats calculation
      const DbamHeader* dbam);

  uint32_t getReadsPerInterval()
  {
    if (m_initState == SENDING) return static_cast<uint32_t>(INIT_INTERVAL_SIZE);

    return m_readsPerInterval;
  }

  uint32_t getInsertSize(     // Calculate the insert size for the given read pair.
      const DbamHeader* dbh)  // the forward-oriented member of the pair
      ;

  void calculate(  // after sampling is done, calculate the insert-size stats
      const uint8_t interval);

  void logInterval(const uint8_t interval) const;

  double getMeanReadLen() const
  {
    if (!m_numReadsBasesAdded) return 1.0;
    return (static_cast<double>(m_basesAdded) / m_numReadsBasesAdded);
  }

  void fetchInserts(const uint8_t interval, StatsInterval::Inserts_c& inserts);

  uint32_t roundUpToMultiple(uint32_t numToRound, uint32_t multiple);

private:
  typedef std::vector<StatsInterval*> Interval_c;
  typedef std::deque<InputDbamRecord> InputRecord_c;

private:
  enum { DUMMY_READ_SIZE = 128 };
  InputDbamRemapper*        m_readSender;
  char                      m_dummyRead[DUMMY_READ_SIZE];  // raw memory for the "init-flush" read
  boost::mutex              m_mutex;                       // protect the condition variable below
  boost::condition_variable m_calculatedAnInterval;        // signal when an interval becomes valid
  boost::mutex              m_init_mutex;                  // protect the condition variable below
  boost::condition_variable m_finishedInitialization;      // signal when we got through init run

  uint64_t m_numStalls;  // How many times have we stalled waiting for interval
                         // to become valid.

  //--------------------------------------------------------------------------------adamb
  // The first 100,000 reads we send to interval 0, but save a copy for remapping
  // for the moment when we get those results.  This allows us to avoid bad mappings
  // due to bad stats for the first reads in each input data set
  InitializationState m_initState;
  SenderState         m_senderState;
  bool                m_reportedWaiting;  // whether we told anyone we're already waiting
  InputRecord_c       m_remapRecords;

  const uint16_t m_readGroupId;       // index of the read group
  const bool     m_rnaMode;           // true if we're running RNA, false for DNA
  const bool     m_continuousUpdate;  // whether to contuously update stats after the initial interval
  const bool     m_updateLogOnly;     // whether to apply continuous updates (false) or just log them (true)
  const uint32_t m_mapqMin;           // minimum MAPQ required to include a read in the stats
  const uint32_t m_peOrientation;     // the paired-end orientation flag, plumbed from the config file
  const uint32_t m_readsPerInterval;  // recalculate stats each time we've received this many reads
  const uint32_t m_sampleSize;        // each time we recalculate stats, use this many past pairs

  Interval_c m_intervals;  // stats intervals
  const uint8_t
                 m_intervalDelay;  // For input side: delay interval we're sending and the one we grab stats from
  const size_t   INIT_INTERVAL_SIZE;
  const uint8_t  m_maxIntervalMemory;   // max # of intervals to include in each calculation
  uint8_t        m_sendingInterval;     // which interval are we currently collecting stats for
  uint64_t       m_basesAdded;          // total number of bases we've added to stats since the start
  uint64_t       m_numReadsBasesAdded;  // number of reads whose bases we've added to stats since the start
  std::ostream&  m_logStream;           // where to log info on stats calculated
  bool           m_printedFirstDetectedStats;
  StatsInterval* m_fixedStats;       // use these unchanging stats instead of automatic detection
  uint64_t       m_dummySendCount;   // total number of dummy records sent
  uint64_t       m_recordsReceived;  // total number of non-flush records processed by insert stats
  uint64_t       m_flushesReceived;  // total number of flush records processed by insert stats
  std::vector<uint64_t> m_insertSizeHistogram;
};

#endif
