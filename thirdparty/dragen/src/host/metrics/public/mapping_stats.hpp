// Copyright 2018 Edico Genome Corporation. All rights reserved.
// Copyright 2018-2020 Illumina, Inc. All rights reserved.

#ifndef __MA_STATS__
#define __MA_STATS__

#include <string>

#define _TARGET_X86_
#if defined(_TARGET_X86_)
#include <mmintrin.h>
#include <xmmintrin.h>
#elif defined(_TARGET_PPC_)
/*#include "mmintrin_types.h"
#include "vec128int.h"
#include "vec_defines.h"
#include "vecs_256_and_64_bit.h"
#undef bool
#undef vector*/
#elif defined(_TARGET_ARM_)
#include "SSE2NEON.h"
#else
#error Target not recognized, porting needed
#endif
#include "infra_compiler.h"

#include "print_metrics.hpp"
#include "output_dbam_header.hpp"
const uint8_t MAPQ_HIST_NR_BINS = 5;

//--------------------------------------------------------------------------------
// PairedEndInsertSizeStats - a small struck to keep track of the estimated insert sizes per per read group
struct RGinsertSizeStats {
  RGinsertSizeStats()
    : m_errMsg(""),
      m_numInsertsInMean(0),
      m_q25(0),
      m_q50(0),
      m_q75(0),
      m_mean(0),
      m_stddev(0),
      m_low(0),
      m_high(0),
      m_minInsert(0),
      m_maxInsert(0)
  {
  }
  std::string m_errMsg;
  int         m_numInsertsInMean;
  int         m_q25;
  int         m_q50;
  int         m_q75;
  double      m_mean;
  double      m_stddev;
  int         m_low;
  int         m_high;
  uint32_t    m_minInsert;
  uint32_t    m_maxInsert;
};

// ------------------------------------------------------------------
// Types of record sets; Used to update all or a subset of the record counts.
// TODO: what we really want longer-term is multiple separate
// partial-aggregate classes, each with its own update() method.
// Call them AlignmentStats, and DupmarkStats, say. Then the RunStats
// class knows how to aggregate each of these specialized stats-gatherers.
enum RecordCountType {
  COUNT_ALL_BUT_DUPS = 1,                                      // count all records except duplicates
  COUNT_DUPLICATES   = 2,                                      // count duplicates
  COUNT_ALL          = COUNT_ALL_BUT_DUPS | COUNT_DUPLICATES,  // count all records
};

//--------------------------------------------------------------------------------
// ReadGroupAlignmentCounts - a class that tracks map/align statistics for an
// individual read group.
//
struct ReadGroupAlignmentCounts {
public:
  ReadGroupAlignmentCounts(std::ostream&  logStream)
    : m_numRecords(0),
      m_numRecordsR1(0),
      m_numRecordsR2(0),
      m_total_mapq_gt_0(0),
      m_numDuplicatesMarked(0),
      m_numDuplicatesRemoved(0),
      m_dup_mapq_gt_0(0),
      m_numSecondary(0),
      m_numSupplementary(0),
      m_unmapped(0),
      m_filteredContig(0),
      m_nonrefDecoy(0),
      m_unmappedR1(0),
      m_unmappedR2(0),
      m_pairedReadDiffContig(0),
      m_pairedReadDiffContigMapq10(0),
      m_pairedDiscordant(0),
      m_concordantMultiple(0),
      m_concordantOnce(0),
      m_unpairedMultiple(0),
      m_unpairedOnce(0),
      m_QCfailed(0),
      m_singleton(0),
      m_numValidRecords(0),
      m_sumSeqLength(0),
      m_sumSeqLengthR1(0),
      m_sumSeqLengthR2(0),
      m_sumUnmappedSeqLengthR1(0),
      m_sumUnmappedSeqLengthR2(0),
      m_numSoftClippedR1(0),
      m_numSoftClippedR2(0),
      m_numIndelReadsR1(0),
      m_numIndelReadsR2(0),
      m_splicedReads(0),
      m_numIndelBasesR1(0),
      m_numIndelBasesR2(0),
      m_numMismatchesR1(0),
      m_numMismatchesR2(0),
      m_numAllQ30BasesR1(0),
      m_numAllQ30BasesR2(0),
      m_numNonDupNonClippedQ30Bases(0),
      m_hasMate(0),
      m_suppressed(0),
      m_isTumor(false),
      m_logStream(logStream)
  {
  }

  //--------------------------------------------------------------------------------
  // update alignment counts for someone who's converting dbam to bam
  // TODO: what we really want longer-term is multiple separate partial-aggregate
  // classes, each with its own update() method.  Call them AlignmentStats,
  // and DupmarkStats, say.  Then the RunStats class knows how to aggregate each
  // of these specialized stats-gatherers.
  //
  void update(const DbamHeader* dbh, RecordCountType type, const bool extendedBamMetrics);
  void updateDups(const DbamHeader* dbh, const bool extendedBamMetrics);
  void updateNonDups(const DbamHeader* dbh, const bool extendedBamMetrics);

  // update alignment stats based on a particular record, but decrement counters
  // void decrement(const DbamHeader* dbh)  // the record being processed
  // ;

  //bridge function for dragmap
  void addRecord(const dragenos::align::SerializedAlignment& alignment, const dragenos::sequences::SerializedRead& read)
  {
    const DbamHeader dbh(alignment,read);
    this->update(&dbh,COUNT_ALL,true);
  }

  void add(const ReadGroupAlignmentCounts& other)
  {
    m_numRecords += other.m_numRecords;  // total number of reads - suppressed reads
    m_numRecordsR1 += other.m_numRecordsR1;
    m_numRecordsR2 += other.m_numRecordsR2;
    m_total_mapq_gt_0 += other.m_total_mapq_gt_0;
    m_numDuplicatesMarked += other.m_numDuplicatesMarked;
    m_numDuplicatesRemoved += other.m_numDuplicatesRemoved;
    m_dup_mapq_gt_0 += other.m_dup_mapq_gt_0;
    m_numSecondary += other.m_numSecondary;
    m_numSupplementary += other.m_numSupplementary;
    m_unmapped += other.m_unmapped;
    m_filteredContig += other.m_filteredContig;
    m_nonrefDecoy += other.m_nonrefDecoy;
    m_unmappedR1 += other.m_unmappedR1;
    m_unmappedR2 += other.m_unmappedR2;
    m_pairedReadDiffContig += other.m_pairedReadDiffContig;
    m_pairedReadDiffContigMapq10 += other.m_pairedReadDiffContigMapq10;
    m_pairedDiscordant += other.m_pairedDiscordant;
    m_concordantMultiple += other.m_concordantMultiple;
    m_concordantOnce += other.m_concordantOnce;
    m_unpairedMultiple += other.m_unpairedMultiple;
    m_unpairedOnce += other.m_unpairedOnce;
    m_QCfailed += other.m_QCfailed;
    for (int i = 0; i < MAPQ_HIST_NR_BINS; i++) {
      m_mapq_hist[i] += other.m_mapq_hist[i];
    }
    m_singleton += other.m_singleton;
    m_numValidRecords += other.m_numValidRecords;
    m_sumSeqLength += other.m_sumSeqLength;
    m_sumSeqLengthR1 += other.m_sumSeqLengthR1;
    m_sumSeqLengthR2 += other.m_sumSeqLengthR2;
    m_sumUnmappedSeqLengthR1 += other.m_sumUnmappedSeqLengthR1;
    m_sumUnmappedSeqLengthR2 += other.m_sumUnmappedSeqLengthR2;
    m_numSoftClippedR1 += other.m_numSoftClippedR1;
    m_numSoftClippedR2 += other.m_numSoftClippedR2;
    m_numIndelReadsR1 += other.m_numIndelReadsR1;
    m_numIndelReadsR2 += other.m_numIndelReadsR2;
    m_splicedReads += other.m_splicedReads;
    m_numIndelBasesR1 += other.m_numIndelBasesR1;
    m_numIndelBasesR2 += other.m_numIndelBasesR2;
    m_numMismatchesR1 += other.m_numMismatchesR1;
    m_numMismatchesR2 += other.m_numMismatchesR2;
    m_numAllQ30BasesR1 += other.m_numAllQ30BasesR1;
    m_numAllQ30BasesR2 += other.m_numAllQ30BasesR2;
    m_numNonDupNonClippedQ30Bases += other.m_numNonDupNonClippedQ30Bases;
    m_hasMate += other.m_hasMate;
    m_suppressed += other.m_suppressed;
    m_isTumor = other.m_isTumor;
  }

  void reset()
  {
    m_numRecords                 = 0;
    m_numRecordsR1               = 0;
    m_numRecordsR2               = 0;
    m_total_mapq_gt_0            = 0;
    m_numDuplicatesMarked        = 0;
    m_numDuplicatesRemoved       = 0;
    m_dup_mapq_gt_0              = 0;
    m_numSecondary               = 0;
    m_numSupplementary           = 0;
    m_unmapped                   = 0;
    m_filteredContig             = 0;
    m_nonrefDecoy                = 0;
    m_unmappedR1                 = 0;
    m_unmappedR2                 = 0;
    m_pairedReadDiffContig       = 0;
    m_pairedReadDiffContigMapq10 = 0;
    m_pairedDiscordant           = 0;
    m_concordantMultiple         = 0;
    m_concordantOnce             = 0;
    m_unpairedMultiple           = 0;
    m_unpairedOnce               = 0;
    m_QCfailed                   = 0;
    for (int i = 1; i < MAPQ_HIST_NR_BINS; i++) {
      m_mapq_hist[i] = 0;
    }
    m_singleton                   = 0;
    m_numValidRecords             = 0;
    m_sumSeqLength                = 0;
    m_sumSeqLengthR1              = 0;
    m_sumSeqLengthR2              = 0;
    m_sumUnmappedSeqLengthR1      = 0;
    m_sumUnmappedSeqLengthR2      = 0;
    m_numSoftClippedR1            = 0;
    m_numSoftClippedR2            = 0;
    m_numIndelReadsR1             = 0;
    m_numIndelReadsR2             = 0;
    m_splicedReads                = 0;
    m_numIndelBasesR1             = 0;
    m_numIndelBasesR2             = 0;
    m_numMismatchesR1             = 0;
    m_numMismatchesR2             = 0;
    m_numAllQ30BasesR1            = 0;
    m_numAllQ30BasesR2            = 0;
    m_numNonDupNonClippedQ30Bases = 0;

    m_hasMate    = 0;
    m_suppressed = 0;
  }

  uint64_t totalNumPaired() const
  {
    return (m_pairedDiscordant + m_concordantMultiple + m_concordantOnce) / 2;
  }


  void setIsTumor(bool isTumor) { m_isTumor = isTumor; }

  bool isTumor() const { return m_isTumor; }

  // void setExtendedMetrics(const bool enable ) { m_extendedMetrics = enable; };

  // Record counts, for summary output
  uint64_t m_numRecords;
  uint64_t m_numRecordsR1;
  uint64_t m_numRecordsR2;
  uint64_t m_total_mapq_gt_0;
  uint64_t m_numDuplicatesMarked;
  uint64_t m_numDuplicatesRemoved;
  uint64_t m_dup_mapq_gt_0;
  uint64_t m_numSecondary;
  uint64_t m_numSupplementary;
  uint64_t m_unmapped;
  uint64_t m_filteredContig;
  uint64_t m_nonrefDecoy;
  uint64_t m_unmappedR1;
  uint64_t m_unmappedR2;
  uint64_t m_pairedReadDiffContig;
  uint64_t m_pairedReadDiffContigMapq10;
  uint64_t m_pairedDiscordant;
  uint64_t m_concordantMultiple;
  uint64_t m_concordantOnce;
  uint64_t m_unpairedMultiple;
  uint64_t m_unpairedOnce;
  uint64_t m_QCfailed;
  uint64_t m_mapq_hist[MAPQ_HIST_NR_BINS] = {0};
  uint64_t m_singleton;
  uint64_t m_numValidRecords;
  uint64_t m_sumSeqLength;
  uint64_t m_sumSeqLengthR1;
  uint64_t m_sumSeqLengthR2;
  uint64_t m_sumUnmappedSeqLengthR1;
  uint64_t m_sumUnmappedSeqLengthR2;
  uint64_t m_numSoftClippedR1;  // Number of soft-clipped bases for R1
  uint64_t m_numSoftClippedR2;  // Number of soft-clipped bases for R2
  uint64_t m_numIndelReadsR1;   // Number of R1 reads containing indels
  uint64_t m_numIndelReadsR2;   // Number of R2 reads containing indels
  uint64_t m_splicedReads;      // For RNA data, reads aligning across a splice junction
  uint64_t m_numIndelBasesR1;   // Number of R1 reads containing indels
  uint64_t m_numIndelBasesR2;   // Number of R2 reads containing indels
  uint64_t m_numMismatchesR1;   // NM edit distance
  uint64_t m_numMismatchesR2;
  uint64_t m_numAllQ30BasesR1;
  uint64_t m_numAllQ30BasesR2;
  uint64_t m_numNonDupNonClippedQ30Bases;

  uint64_t         m_hasMate;
  uint64_t         m_suppressed;
  bool             m_isTumor;
  bool             m_extended_metrics;


void printStats(std::chrono::duration<float>  mapDuration);

private:

  void updateUsingCigar(const DbamHeader* dbh, const bool hasMate, const bool isFirstInPair)
  {
    uint16_t n_softclipped  = 0;      // number of soft-clipped bases
    uint16_t n_indelBases   = 0;      // number of indel bases
    bool     containsSplice = false;  // Read CIGAR spans a refskip (RNA splice junction)
    bool     containsIndel  = false;  // read contains indel?

    // Walk the cigar
    

    const auto& cigar = dbh->getCigar();
    const int16_t   n_cigar_recs = cigar.end() - cigar.begin();
    auto it = cigar.begin();
    
    for (int32_t j = 0;  j < n_cigar_recs; ++j,  ++it ){

  
      auto & operation = *it;
      const unsigned count = operation.second;

      switch (operation.first) {
      case dragenos::align::Cigar::DELETE:
      case dragenos::align::Cigar::INSERT:
        containsIndel = true;
        n_indelBases += count;
        break;
      case dragenos::align::Cigar::SOFT_CLIP:
        n_softclipped += count;
        break;
      case dragenos::align::Cigar::SKIP:
        // Cigar operations of type SKIP consume two records.
        ++it;++j;
        containsSplice = true;
        break;
      default:
        break;
      }  // switch
    }    // for

    if (!hasMate || isFirstInPair) {
      m_numSoftClippedR1 += n_softclipped;
      if (containsIndel) ++m_numIndelReadsR1;
      m_numIndelBasesR1 += n_indelBases;
    } else {
      m_numSoftClippedR2 += n_softclipped;
      if (containsIndel) ++m_numIndelReadsR2;
      m_numIndelBasesR2 += n_indelBases;
    }
    if (containsSplice) ++m_splicedReads;
  }

public:  // public to test countBQ30Bases
  static uint32_t countBQ30Bases(const uint8_t* seq, const uint32_t len)
  {
    uint32_t n_q30      = 0;
    uint32_t bytes_left = len;


#ifndef _TARGET_PPC_
    // Process as much of the sequence as we can 16 bytes at a time.
    // We'll be checking for any bases (in DBAM format, so 2 low-order bits of basecall and 6
    // high-order bits of quality) that are >= qual 30.
    static const __m128i QUAL_30 = _mm_set1_epi8(30); //GR for dragmap we get actual qualities

    // Set up loop to grab 16 bytes at a time:
    const __m128i* seqptr = reinterpret_cast<const __m128i*>(seq);

    for (; bytes_left >= 16; bytes_left -= 16) {
      // Stash a copy of the current 16 bytes
      const __m128i S = _mm_loadu_si128(seqptr++);

      // For each byte, check if its qual is >= 30
      // Note that (x >= y) is equivalent to (x == max(x,y)).  Cool, right?
      const __m128i highbytes = _mm_cmpeq_epi8(_mm_max_epu8(S, QUAL_30), S);

      // The previous operation checks each byte, and for each byte that evaluated
      // true, it sets the resulting output byte to 0xFF.  0 otherwise.  So we
      // can now count how many bits are set.
      const uint64_t* lo = reinterpret_cast<const uint64_t*>(&highbytes);
      const uint64_t* hi = lo + 1;
      n_q30 += __POPCOUNTLL(*lo) + __POPCOUNTLL(*hi);
    }

    // The previous loop counted 8 bits for each byte that evaluated as high quality.
    n_q30 /= 8;
#endif



    // Count over the final (up to) 8 bases one by one:
    for (uint32_t i = len - bytes_left; i < len; ++i) {
      if (seq[i] >= 30) n_q30++;
    }

    return n_q30;
  }

private:
  std::ostream  & m_logStream;           // where to log info on stats calculated


void printMapAlignCommonStats(
    PrintMetrics&             ma_printer,
    ReadGroupAlignmentCounts& stats,
    const std::string&        runName,
    const std::string&        sectionName,
    const std::string&        sampleName,
    const std::string&        inputReadsDescription,
    const bool                setPrintedSampleName,
    float                     mapDurationSeconds);


  void updateNonDupsBQ30counts(const DbamHeader* dbh)
  {
   // const int16_t   n_cigar_recs = static_cast<int16_t>(dbh->getCigarLen() / 2);
   // const uint16_t* pCigar       = dbh->getConstCigarRecordPtr();
    const uint8_t*  seq          = dbh->getConstQualities();
    uint32_t        q30_bases    = 0;
    uint32_t        inSeqIdx     = 0;


    const auto& cigar = dbh->getCigar();
    const int16_t   n_cigar_recs = cigar.end() - cigar.begin();
    auto it = cigar.begin();
        
    for (int32_t j = 0;  j < n_cigar_recs; ++j,  ++it ){
       auto & operation = *it;
       const unsigned oplen = operation.second;

      switch (operation.first) {
      case dragenos::align::Cigar::ALIGNMENT_MATCH:
      case dragenos::align::Cigar::SEQUENCE_MISMATCH:
      case dragenos::align::Cigar::SEQUENCE_MATCH:
      case dragenos::align::Cigar::INSERT:
        q30_bases += countBQ30Bases(seq + inSeqIdx, oplen);
        inSeqIdx += oplen;
        break;
      case dragenos::align::Cigar::SOFT_CLIP:
        inSeqIdx += oplen;
        break;
      case dragenos::align::Cigar::SKIP:
        // Cigar operations of type SKIP consume two records.
        ++it;j++;
        break;
      case dragenos::align::Cigar::DELETE:
      case dragenos::align::Cigar::HARD_CLIP:
      case dragenos::align::Cigar::PAD:
      default:
        break;
      }
    }

    m_numNonDupNonClippedQ30Bases += q30_bases;
  }

};

#endif

