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

#include "readgroup_insert_stats.hpp"

#include <algorithm>
#include <ctime>
//#include "api/BamConstants.h"
#include "dragen_run_log.hpp"
#include "input_dbam_record.hpp"
#include "output_dbam_header.hpp"
#include "read_group_list.hpp"
#include "run_stats.hpp"

//--------------------------------------------------------------------------------adamb
// constructor
//
ReadGroupInsertStats::ReadGroupInsertStats(
    const uint16_t rgid,              // index of the read group
    const bool     rnaMode,           // true for RNA, false for DNA
    const bool     continuousUpdate,  // whether to update stats continuously
    const bool     updateLogOnly,     // whether updates should just be logged, or also applied to input reads
    const uint32_t readsPerInterval,  // how many reads to accumulate per interval
    const uint32_t sampleSize,        // how many inserts to include in each calculation
    const uint8_t  maxIntervalMemory,  // max # of intervals to include in each calculation
    const uint32_t intervalDelay,      // lag between calculating and using stats
    const uint32_t mapqMin,            // minimum MAPQ required for inclusion in stats
    const uint32_t peOrientation,      // expected orientation of the pairs
    const double   rescueSigmas,       // num std deviations spread to allow rescue alignments
    const double   rescueCeilFactor,   // multiple of read len to compare to rescue sigmas
    const uint32_t rescueMinInsert,    // override of rescue min insert length
    const uint32_t rescueMaxInsert,    // override of rescue max insert length
    std::ostream&  logStream)           // where to log info on stats
  : m_readSender(0),
    m_numStalls(0),
    m_initState(SENDING),
    m_senderState(NORMAL),
    m_reportedWaiting(false),
    m_readGroupId(rgid),
    m_rnaMode(rnaMode),
    m_continuousUpdate(continuousUpdate),
    m_updateLogOnly(updateLogOnly && !m_rnaMode),
    m_mapqMin(mapqMin),
    m_peOrientation(peOrientation),
    m_readsPerInterval(readsPerInterval),
    m_sampleSize(sampleSize),
    m_intervals(std::numeric_limits<uint8_t>::max()),  // max id is 254
    m_intervalDelay(intervalDelay),
    INIT_INTERVAL_SIZE(m_readsPerInterval * (m_intervalDelay - 1)),
    m_maxIntervalMemory(maxIntervalMemory),
    m_sendingInterval(0),
    m_basesAdded(0),
    m_numReadsBasesAdded(0),
    m_logStream(logStream),
    m_printedFirstDetectedStats(false),
    m_fixedStats(0),
    m_dummySendCount(0),
    m_recordsReceived(0),
    m_flushesReceived(0),
    m_insertSizeHistogram(INITIAL_HISTOGRAM_SIZE, 0)
{
  for (size_t i = 0; i < m_intervals.size(); ++i) {
    m_intervals[i] = new StatsInterval(
        rgid, i, rnaMode, m_peOrientation, rescueSigmas, rescueCeilFactor, rescueMinInsert, rescueMaxInsert);
  }

  // Mark the first couple of intervals as valid, with default stats, to get us
  // started.  After the initialization phase has completed, these default stats
  // will be overwritten with the initial sampled stats.
  const size_t to   = std::numeric_limits<uint8_t>::max();
  const size_t from = to - m_intervalDelay - 1;
  for (size_t i = from; i < to; ++i) {
    calculate(i);
  }
}

//--------------------------------------------------------------------------------adamb
// Destructor
//
ReadGroupInsertStats::~ReadGroupInsertStats()
{
  for (auto it = m_intervals.begin(); it != m_intervals.end(); ++it) {
    delete (*it);
  }

  if (m_fixedStats) delete m_fixedStats;
}

//--------------------------------------------------------------------------------adamb
// Either completely disable rescue alignments, or leave them enabled such that they
// will be applied when the stats justify their use
//
void ReadGroupInsertStats::enableRescueAlignments(const bool enableit)
{
  for (auto it = m_intervals.begin(); it != m_intervals.end(); ++it) {
    (*it)->enableRescueAlignments(enableit);
  }
}

//--------------------------------------------------------------------------------cobus
// Send a bunch of non-mapping, non-outputtable records to flush out the
// insert-stats system. Sending 64KB chunks at a time and waiting in between.
// Poll/timeout/assert if we send too many and do not get valid interval
//
// There will be records in Bufferizer to flush out to S2C, and once those records
// have been mapped we also flush until partially filled C2S is output
//
void ReadGroupInsertStats::flushUntilIntervalValid(uint32_t useInterval)
{
  uint32_t wait_timeout = 0;

  while (!m_intervals[useInterval]->isValid()) {
    m_senderState = STALLED;

    InputDbamRecord dummy = InputDbamRecord(reinterpret_cast<uint8_t*>(m_dummyRead));
    for (size_t i = 0; i < 65536; i += DUMMY_READ_SIZE) {
      // send 64KB of dummy reads at a time
      m_readSender->sendRead(dummy, false);
      m_dummySendCount++;
    }

    std::clock_t timestamp = std::clock();
    // wait for 100 msec. arbitrary time to allow some HW processing to finish
    while ((std::clock() >= timestamp) && (std::clock() < (timestamp + CLOCKS_PER_SEC / 10)))
      ;
    // limit to 100MB of dummy records and ~2.7 minutes
    ++wait_timeout;
    assert(wait_timeout < 1600);
  }
  m_numStalls += m_senderState;
  m_senderState = NORMAL;
}

//--------------------------------------------------------------------------------adamb
// Copy the correct insert stats into a mapper input record
// Caller must ensure that the stats interval is valid, via prior call to
// waitForValidInterval()
//
void ReadGroupInsertStats::fillPeInsertStats(InputDbamRecord& dbr)
{
  // Check if we've overridden the automatic stats detection for this read group:
  if (m_fixedStats && !m_updateLogOnly) {
    dbr.getBlobs().getFixedBlob()->peStatsInterval = 0;
    m_fixedStats->fillPeInsertStats(dbr);
    return;
  }

  // m_sendingInterval indicates which interval is going to get to tally
  // up the stats for the current record, dbr.
  //
  // useInterval indicates which interval dbr should use to get stats
  // for its own mapping.
  //
  // First check if we can continue sending to the current interval, or
  // if it's time to move on to the next one:

  // All initialization records should have been sent, and
  // should have arrived on the receive side, and init stats should have been
  // valid before calling this API
  assert(m_initState != WAITING);

  // Keep a running tally of the number of bases added, so we can find out the mean
  // read length:
  m_basesAdded += dbr.getHeader()->getSequenceLen();
  ++m_numReadsBasesAdded;

  // If we've sent enough for one interval, advance to the next:
  if (m_intervals[m_sendingInterval]->getNumSent() == getReadsPerInterval()) {
    m_intervals[m_sendingInterval]->setDoneSending();
    m_sendingInterval = (m_sendingInterval + 1) % m_intervals.size();
    m_intervals[m_sendingInterval]->reset();
  }

  dbr.getBlobs().getFixedBlob()->peStatsInterval = m_sendingInterval;
  m_intervals[m_sendingInterval]->incNumSent();

  // We use stats from "m_intervalDelay" intervals ago:
  uint32_t useInterval = (m_sendingInterval + m_intervals.size() - m_intervalDelay) % m_intervals.size();

  // Copy the actual stats into the input record:
  StatsInterval* pint = (m_fixedStats /* && m_updateLogOnly*/) ? m_fixedStats : m_intervals[useInterval];

  // If we're still in "initialization mode" then save up the record to remap
  // after we've received proper stats
  if (m_initState != SENDING) {
    // Not initializing.  The interval must be valid before calling this API
    // valid interval is ensured by prior call to waitForValidInterval()
    assert(pint->isValid());
  }
  pint->fillPeInsertStats(dbr);
}

//--------------------------------------------------------------------------------cobus
// Wait for valid stats from prior interval to become valid
// This function will wait and flush with dummy records until the interval
// is valid
//
void ReadGroupInsertStats::waitForValidInterval()
{
  if (!m_fixedStats) {
    // m_sendingInterval indicates which interval is going to get to tally
    // up the stats for the current record, dbr.
    //
    // useInterval indicates which interval dbr should use to get stats
    // for its own mapping.
    //
    if (m_initState == WAITING) {
      // All initialization records have been sent, but haven't yet
      // arrived on the receive side.  Wait until init stats are valid.
      waitUntilInitDone();
    }

    // We use stats from "m_intervalDelay" intervals ago:
    uint32_t useInterval =
        (m_sendingInterval + 1 + m_intervals.size() - m_intervalDelay) % m_intervals.size();

    if (m_initState != SENDING) {
      // Not initializing.
      // Wait until the interval is valid, before we use it in fillPeInsertStats
      // If we have a readSender, flush with dummy records and poll
      // If not, we cannot flush, just wait.
#if defined(DMA_REAL_DRIVER)
      if (m_readSender) {
        flushUntilIntervalValid(useInterval);
      } else {
#endif
        boost::unique_lock<boost::mutex> lock(m_mutex);
        while (!m_intervals[useInterval]->isValid()) {
          ++m_numStalls;
          m_calculatedAnInterval.wait(lock);
        }
#if defined(DMA_REAL_DRIVER)
      }
#endif
    }
  }
}

//--------------------------------------------------------------------------------adamb
// Each read group's stats are initialized with the first 100,000 input reads.
// After the stats are calculated, we remap those records.  This is where
// we save a copy of the records for remapping later.
//
void ReadGroupInsertStats::saveForRemapping(InputDbamRecord& dbr)
{
  if (m_initState != SENDING || !dbr.hasMate()) return;

  // Save a copy of the record for remapping later.  Whoever remaps it will
  // be responsible for freeing it.
  char* recCopy = reinterpret_cast<char*>(malloc(dbr.getRecordLen()));
  memcpy(recCopy, dbr.getRawData(), dbr.getRecordLen());
  m_remapRecords.push_back(InputDbamRecord(reinterpret_cast<uint8_t*>(recCopy)));

  // Because the record will be remapped later when we have good stats, make sure
  // that it doesn't get sent to output after its trip through the mapper/aligner.
  dbr.getHeader()->flag |= DbamHeader::SUPPRESS_OUTPUT;

  // If we just completed the initialization, switch into a state
  // where future reads will wait for the stats to flush.
  if (m_intervals[0]->getNumSent() == getReadsPerInterval()) {
    setInitDoneSending();
  }
}

//--------------------------------------------------------------------------------adamb
// If the input file was really small, we need to get the init phase to conclude,
// wait for initial stats to be computed, and prepare to remap any saved up
// records.
void ReadGroupInsertStats::setInitDoneSending()
{
  boost::unique_lock<boost::mutex> lock(m_init_mutex);

  if (m_initState != SENDING) return;

  // We have not yet sent 100,000 'init' records with default stats.  Regardless,
  // that's all we're going to get, so tell the receive side not to expect
  // any more stats to arrive after the current set flushes.
  m_intervals[m_sendingInterval]->setDoneSending();
  if (0 == m_intervals[m_sendingInterval]->getNumSent()) {
    m_initState = DONE;
    m_finishedInitialization.notify_all();
    if (m_rnaMode || !m_continuousUpdate || m_updateLogOnly) {
      ASSERT(!m_fixedStats);
      m_fixedStats = new StatsInterval(*m_intervals[m_sendingInterval]);
    }
  } else {
    m_initState = WAITING;
  }
}

//-------------------------------------------------------------------------------swhitmore/adamb
// Returns true if the read #dbam# should be used in calculating sample stats
//
bool ReadGroupInsertStats::shouldUseRead(const DbamHeader* dbam)
{
  // Filter out any reads that don't match the correct flags
  static const uint16_t DONT_USE_ME_MASK =
      BamTools::Constants::BAM_ALIGNMENT_UNMAPPED | BamTools::Constants::BAM_ALIGNMENT_SECONDARY |
      static_cast<uint16_t>(DbamHeader::ALIGNMENT_FLAG_SUPPLEMENTARY) |
      BamTools::Constants::BAM_ALIGNMENT_QC_FAILED | BamTools::Constants::BAM_ALIGNMENT_DUPLICATE |
      BamTools::Constants::BAM_ALIGNMENT_READ_2;  // second-in-pair doesn't have a target interval

  // Bail if any of the "don't use" flags are turned on
  if (dbam->getFlag() & DONT_USE_ME_MASK) {
    return false;
  }

  // Bail if it's single ended, or mate is unmapped, or if it's part of an improper pair.
  if (!dbam->hasMate() || dbam->isMateUnmapped() || !dbam->isProperlyPaired()) {
    return false;
  }

  // Bail if it's a low-quality mapping
  if (dbam->getMapQuality() < m_mapqMin) {
    return false;
  }

  return true;
}

//--------------------------------------------------------------------------------adamb
// Check if the header represents a good pair that we want to include in the stats.
// If so, add its info to our aggregate record.
//
void ReadGroupInsertStats::sample(const DbamHeader* dbh)
{
  // If we've overridden automatic stats detection, bail.
  if (m_fixedStats && !m_updateLogOnly) return;

  // Keep track of the number of primary alignments returned by the aligner, which
  // should match the number of reads sent into the aligner.
  if (dbh->isPrimary() && dbh->isFirstInPair()) {
    // record counts for debug
    if (dbh->suppressOutput()) {
      m_flushesReceived++;
    } else {
      m_recordsReceived++;
    }

    // Figure out which stats interval this record should be tallied in, and add it:
    const uint8_t interval = dbh->getPeStatsInterval();
    if (shouldUseRead(dbh)) {
      uint32_t insertSize = getInsertSize(dbh);
      while (insertSize >= m_insertSizeHistogram.size()) {
        m_insertSizeHistogram.resize(m_insertSizeHistogram.size() * 2, 0);
      }
      m_insertSizeHistogram[insertSize]++;
      m_intervals[interval]->addInsert(insertSize);
    }

    // First end of a pair had insert stats on the "send" side, so incremented the
    // stats interval's m_numSend:
    m_intervals[interval]->incNumReceived();

    checkIntervalForCompletion(interval);
  }
}

//--------------------------------------------------------------------------------adamb
// Check if a particular interval has received all of its data.  If so, calculate
// its stats and do whatever other followup is required, e.g. printing reports
// and noting that we have initial stats so remapping of sample records can proceed.
//
void ReadGroupInsertStats::checkIntervalForCompletion(const uint8_t interval)
{
  // Once we've received the primary alignment for all of the reads that were sent
  // for an interval, we can calculate the resulting stats.
  if (m_intervals[interval]->receivedAll()) {
    calculate(interval);

    if (!m_printedFirstDetectedStats) {
      printConfig(std::cout, interval);
      m_printedFirstDetectedStats = true;
    }

    // The first 100000 reads are sent with default stats into the
    // m_intervalDelay'th interval.   When we've received the last of them,
    // we calculate the stats for that sample, and copy the results into
    // the first n intervals.
    if (m_initState == WAITING) {
      completeInitialization();
    }
  }
}

//--------------------------------------------------------------------------------adamb
// Bug 1950 uncovered an interesting situation, where a small read group near the
// end of multi-RG input file would never complete initialization. Turned out this
// was because that read group would receive all of its output records BEFORE
// we turn on StatsInterval::m_doneSending.  The fix is to put an extra check
// in AlignedReadsSampler::sample(), to make sure we don't get stuck in that state.
//
void ReadGroupInsertStats::checkForInitComplete()
{
  if (m_initState != WAITING) {
    return;
  }

  checkIntervalForCompletion(0);
}

//--------------------------------------------------------------------------------adamb
// After the initialization reads have all been sent we wait for them to make
// it through the mapper.
//
void ReadGroupInsertStats::completeInitialization()
{
  if (m_initState != WAITING) return;

  const size_t to   = std::numeric_limits<uint8_t>::max();
  const size_t from = to - m_intervalDelay - 1;

  for (uint8_t i = from; i < to; ++i) {
    m_intervals[i]->copyStats(*m_intervals[0]);
  }

  // For RNA mode, we need to keep the same stats throughout multiple
  // runs.  So we have to switch the stats-collector into fixed stats mode.
  if (m_rnaMode || !m_continuousUpdate || m_updateLogOnly) {
    assert(!m_fixedStats);
    m_fixedStats = new StatsInterval(*m_intervals[from]);
  }

  // Now set things up so that we will start re-sending any initialization
  // reads, but this time using the stats from the intervals "from" to "to"
  m_sendingInterval = 0;
  m_intervals[m_sendingInterval]->reset();
  m_intervals[m_sendingInterval]->clearInserts();
  boost::unique_lock<boost::mutex> lock(m_init_mutex);
  m_initState = DONE;
  m_finishedInitialization.notify_all();
}

//--------------------------------------------------------------------------------adamb
// When we know we've sent all of the records for a read group, someone may need
// to wait until the results have arrived on the other side.  This does that.
//
void ReadGroupInsertStats::waitUntilInitDone()
{
  if (m_initState == DONE) return;

  // We're done sending the init input records, but they haven't come back from
  // the board yet. Wait until they do.
  boost::unique_lock<boost::mutex> lock(m_init_mutex);
  while (m_initState == WAITING) {
    m_finishedInitialization.wait(lock);
  }
}

//--------------------------------------------------------------------------------adamb
// Return a number that is #numToRound#, rounded up to the nearest
// multiple of #multiple#
uint32_t ReadGroupInsertStats::roundUpToMultiple(uint32_t numToRound, uint32_t multiple)
{
  if (multiple == 0) return numToRound;

  uint32_t remainder = numToRound % multiple;
  if (remainder == 0) return numToRound;

  return numToRound + multiple - remainder;
}

//--------------------------------------------------------------------------------adamb
// Calculate the stats for the specified interval
//
void ReadGroupInsertStats::calculate(const uint8_t interval)
{
  if (m_intervals[interval]->isValid()) {
    return;
  }

  // Get the ~100000 most recent inserts, which may be spread across some number
  // of intervals, and calculate stats based on that set.
  StatsInterval::Inserts_c inserts;
  inserts.reserve(roundUpToMultiple(m_sampleSize + 1, m_readsPerInterval));
  fetchInserts(interval, inserts);
  boost::unique_lock<boost::mutex> lock(m_mutex);
  m_intervals[interval]->calculate(inserts, getMeanReadLen());
  logInterval(interval);

  // There is a limit to how far back we'll lock for stats -- we want it to remain
  // responsive.  So roll off the oldest one we're still looking at:
  const uint16_t rolloffInterval = (interval + m_intervals.size() - m_maxIntervalMemory) % m_intervals.size();
  m_intervals[rolloffInterval]->clearInserts();

  m_calculatedAnInterval.notify_all();
}

//--------------------------------------------------------------------------------adamb
// We are ready to calculate the stats after having received all of the expected
// data for a particular #interval#.  Fetch all of the inserts into the #inserts#
// container
void ReadGroupInsertStats::fetchInserts(const uint8_t interval, StatsInterval::Inserts_c& inserts)
{
  // Fetch inserts from full intervals until we've either gotten a full target
  // number, or we've looked as far into the past as we're willing to look.
  // Note that we always have to fetch full intervals, never just partial
  // intervals, because the order of records changes from run to run and
  // we want to have deterministic results.
  for (uint8_t i = 0; i < m_maxIntervalMemory; ++i) {
    int32_t idx = (interval + m_intervals.size() - i) % m_intervals.size();
    if ((i == 0) || m_intervals[idx]->isValid()) {
      auto add_ins = m_intervals[idx]->getInserts();
      inserts.insert(inserts.end(), add_ins.first, add_ins.second);
    }

    if (inserts.size() >= m_sampleSize) break;
  }
}

//--------------------------------------------------------------------------------adamb
// Print out information on a particular interval
//
void ReadGroupInsertStats::logInterval(const uint8_t interval) const
{
  m_intervals[interval]->log(m_logStream);
  m_logStream << m_numStalls << "\t";

  const uint32_t delay = (m_sendingInterval - interval + m_intervals.size()) % m_intervals.size();
  m_logStream << delay << "\t\n";

  m_logStream.flush();
}

//-------------------------------------------------------------------------------adamb
// GetInsertSize - Calculate the insert size for the given read pair, based on
// information contained in the forward-oriented member of the pair
//
uint32_t ReadGroupInsertStats::getInsertSize(const DbamHeader* dbh)
{
  if ((m_peOrientation > 0) && !m_rnaMode) {
    return dbh->getTemplateLen();
  }

  const int64_t myCoord = dbh->getUnclippedAlignmentCoordinate();
  if (myCoord >= static_cast<int64_t>(dbh->getMateCoordinate()))
    return myCoord - dbh->getMateCoordinate() + 1;
  else
    return dbh->getMateCoordinate() - myCoord + 1;
}

//-------------------------------------------------------------------------------swhitmore
// compare - comparison for qsort
//
inline int compare(const void* a, const void* b)
{
  return (*(uint32_t*)a - *(uint32_t*)b);
}

////--------------------------------------------------------------------------------swhitmore/adamb
//// print the stats for this read group to the specified output stream
////
void ReadGroupInsertStats::printConfig(std::ostream& os, const uint8_t interval) const
{
  os << (!m_printedFirstDetectedStats ? "Initial" : "Final")
     << " paired-end statistics detected for read group "
     << GetReadGroupList().getReadGroupName(m_readGroupId);
  m_intervals[interval]->printDescriptiveLog(os);
}

//--------------------------------------------------------------------------------adamb
// When a run is complete, write out the last set of stats we had detected.
//
void ReadGroupInsertStats::finalLog()
{
  if (m_fixedStats) {
    m_fixedStats->updateRunStats();
  }

  if (m_rnaMode) {
    return;
  }

  // Add insert size histogram in RunStats
//  RunStats::Instance()->addInsertSizeHistogram(m_insertSizeHistogram);

  if (m_fixedStats) {
    std::cout << "Paired-end statistics were based on initial reads for read group "
              << GetReadGroupList().getReadGroupName(m_readGroupId);
    m_fixedStats->printDescriptiveLog(std::cout);
    return;
  }

  // Conclude the stats calculation for the final interval.  The only situation in
  // which it would already be valid at this point would be if the last read happened
  // to be the last one needed to complete the stats for the region.
  if (!m_intervals[m_sendingInterval]->isValid()) {
    calculate(m_sendingInterval);
  }

  // Even if we never printed the initial stats (e.g. if we have a small data set),
  // when we print out this info call it "Final" rather than "Initial"
  m_printedFirstDetectedStats = true;
//  printConfig(std::cout, m_sendingInterval);
  m_intervals[m_sendingInterval]->warnIfTooLittleData(std::cout);

  // update runStats to be printed in metrics log --theoh
  m_intervals[m_sendingInterval]->updateRunStats();
}

//--------------------------------------------------------------------------------adamb
// Dump some info to cerr to give a clue after we've hit some kind of serious error
//
void ReadGroupInsertStats::logError() const
{
  std::cerr << "Read group " << m_readGroupId << std::endl;
  std::cerr << "  Sending interval: " << static_cast<uint32_t>(m_sendingInterval) << std::endl;
  m_intervals[m_sendingInterval]->logError();
}

//--------------------------------------------------------------------------------adamb
// In response to SIGUSR1, we dump a bunch of information on the state of the system
// to a hang_dump log file.
void ReadGroupInsertStats::DumpHangInfo(std::ostream& os)
{
  os << "    Read group id: " << m_readGroupId << "\n";
  os << "    RNA mode: " << (m_rnaMode ? "yes" : "no") << "\n";
  os << "    Num prior stalls: " << m_numStalls << "\n";
  os << "    Num flushes sent: " << m_dummySendCount << "\n";
  os << "    Stalled? " << ((m_senderState == STALLED) ? "yes" : "no") << "\n";
  os << "    Num records received: " << m_recordsReceived << "\n";
  os << "    Num flushes received: " << m_flushesReceived << "\n";
  os << "    Init state: ";
  switch (m_initState) {
  case (SENDING):
    os << "SENDING\n";
    break;
  case (WAITING):
    os << "WAITING\n";
    break;
  case (DONE):
    os << "DONE\n";
    break;
  }
  os << "    Reported waiting? " << (m_reportedWaiting ? "yes" : "no") << "\n";
  os << "    Num remap records: " << m_remapRecords.size() << "\n";
  if (m_fixedStats) {
    os << "    USING FIXED STATS\n";
  } else {
    os << "    Sending interval: " << (int)m_sendingInterval << "\n";
  }

  os << "    INTERVALS:\n";
  for (auto it = m_intervals.begin(); it != m_intervals.end(); ++it) {
    (*it)->DumpHangInfo(os);
  }
  os << std::endl;
}
