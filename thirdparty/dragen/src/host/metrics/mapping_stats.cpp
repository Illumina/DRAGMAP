// Copyright 2018 Edico Genome Corporation. All rights reserved.
// Copyright 2018-2020 Illumina, Inc. All rights reserved.

#include "mapping_stats.hpp"
#include "print_metrics.hpp"

void ReadGroupAlignmentCounts::updateDups(const DbamHeader* dbh, const bool extendedBamMetrics)
{
  if (dbh->isSecondary() || dbh->isDisqualified()) return;

  // count non dup Q30 bases ( include supplementary bases )
  if (!dbh->isUnmapped()) {
    // if mapped, then check for and skip duplicates
    if (!dbh->isDuplicate() && extendedBamMetrics) updateNonDupsBQ30counts(dbh);
  } else {
    // unmapped reads cannot be duplicates, supplementary or have soft clipped reads
    // we still want to count bases
  
    const uint8_t* seq       = dbh->getConstQualities();
    const uint32_t seqLength = dbh->getSequenceLen();

    if (extendedBamMetrics) {
      m_numNonDupNonClippedQ30Bases += countBQ30Bases(seq, seqLength);
    }

  }

  if (dbh->isSupplementary()) return;

  // process duplicates
  if (dbh->isDuplicate()) {
    ++m_numDuplicatesRemoved;
    // when DRAGEN mark and removes dups, some unmapped mates end up with the dup flag set to true
    if (!dbh->isUnmapped()) {
      ++m_numDuplicatesMarked;
    }
    if ((dbh->getMapQuality() > 0) && (!(dbh->getFlag() & BamTools::Constants::BAM_ALIGNMENT_UNMAPPED))) {
      ++m_dup_mapq_gt_0;
    }
  }
};

void ReadGroupAlignmentCounts::updateNonDups(const DbamHeader* dbh, const bool extendedBamMetrics)
{
  if (dbh->isDisqualified()) {
    ++m_QCfailed;
    return;
  }

  if (dbh->isSecondary()) {
    ++m_numSecondary;
    return;
  }

  if (dbh->isSupplementary()) {
    ++m_numSupplementary;
    return;
  }

  const uint8_t* seq           = dbh->getConstQualities();
  const uint8_t  mapq          = dbh->getMapQuality();
  const bool     hasMate       = dbh->hasMate();
  const uint32_t seqLength     = dbh->getSequenceLen();
  const bool     isFirstInPair = dbh->isFirstInPair();
  ++m_numRecords;

  if (!hasMate || isFirstInPair) {
    ++m_numRecordsR1;
    m_numMismatchesR1 += dbh->getEditDistance();
  } else {
    ++m_numRecordsR2;
    m_numMismatchesR2 += dbh->getEditDistance();
  }

  // total number of reads mapped, with mapq>0, which are not
  // secondary, supplementary, or disqualified
  if ((mapq > 0) && (!(dbh->getFlag() & BamTools::Constants::BAM_ALIGNMENT_UNMAPPED))) {
    ++m_total_mapq_gt_0;
  }


  // estimate read length
  m_sumSeqLength += seqLength;
  if (!hasMate || isFirstInPair)
    m_sumSeqLengthR1 += seqLength;
  else
    m_sumSeqLengthR2 += seqLength;

  // BQ30 bases for all reads
  uint16_t lcnt = 0;
  if (extendedBamMetrics) {

    lcnt = countBQ30Bases(seq, seqLength);

    if (!hasMate || isFirstInPair) {
      m_numAllQ30BasesR1 += lcnt;
    } else {
      m_numAllQ30BasesR2 += lcnt;
    }
    // end BQ30 bases for all reads

}

  const bool isUnmapped     = dbh->isUnmapped();
  const bool isMateUnmapped = dbh->isMateUnmapped();
  if (hasMate) {
    ++m_hasMate;
    if (!isUnmapped && isMateUnmapped) ++m_singleton;
  }

  if (isUnmapped) {
   //GR I dont think this applies to dragmap
   // if (dbh->getZSTag() == DbamHeader::ZS_FILTERED_CONTIG) ++m_filteredContig;
  //  if (dbh->getZSTag() == DbamHeader::ZS_NONREF_DECOY) ++m_nonrefDecoy;
    ++m_unmapped;
    if (!hasMate || isFirstInPair) {
      ++m_unmappedR1;
      m_sumUnmappedSeqLengthR1 += seqLength;
    } else {
      ++m_unmappedR2;
      m_sumUnmappedSeqLengthR2 += seqLength;
    }
    return;
  }

  // MAPQ histogram
  uint16_t idx = mapq / 10;                                             // floor, e.g. first bin = [0 : 10)
  idx = (idx > MAPQ_HIST_NR_BINS - 1) ? (MAPQ_HIST_NR_BINS - 1) : idx;  // last bin = [mapq lower: inf)
  m_mapq_hist[idx]++;

  // Go over CIGAR
  if (extendedBamMetrics) {
    updateUsingCigar(dbh, hasMate, isFirstInPair);
  }

  // To count multiples, we first check if the NH tag is non-zero. If it is non-zero, then we can use the
  // NH tag to determine if it is a unique or multimapping alignment. If the NH tag is zero, then this
  // field was not populated and we fall back to the old method of using MAPQ as a best guess..
  //const bool multiple = dbh->getNumHits() ? (not dbh->isUnique()) : (mapq == 0);
  //GR falling back to mapq==0 method for dragmap 
  const bool multiple = mapq == 0 ;

  const bool isPaired = dbh->isPairMapped();
  // mate not on same contig
  if (!isUnmapped && isPaired && !isMateUnmapped && !dbh->isMateOnSameContig()) {
    m_pairedReadDiffContig++;
    if (mapq >= 10) {
      m_pairedReadDiffContigMapq10++;
    }
  }

  // const bool concordant = (dbh->isPairMapped() && dbh->isProperlyPaired());
  if (isPaired) {
    if (dbh->isProperlyPaired()) {
      if (multiple)
        ++m_concordantMultiple;
      else
        ++m_concordantOnce;
    } else
      ++m_pairedDiscordant;
  } else {
    if (multiple)
      ++m_unpairedMultiple;
    else
      ++m_unpairedOnce;
  }
};

void ReadGroupAlignmentCounts::update(
    const DbamHeader* dbh, RecordCountType type, const bool extendedBamMetrics)
{
   //GR hmm no idea what it is for, removing
   //  (const_cast<DbamHeader*>(dbh))->clearCountedOnTarget();

  if (type & COUNT_ALL_BUT_DUPS) {
    updateNonDups(dbh, extendedBamMetrics);
  } else if (type & COUNT_DUPLICATES) {
    updateDups(dbh, extendedBamMetrics);
  } else if (type & COUNT_ALL) {
    updateDups(dbh, extendedBamMetrics);
    updateNonDups(dbh, extendedBamMetrics);
  }
}



void ReadGroupAlignmentCounts::printStats(std::chrono::duration<float> mapDuration)
{

   float mapDurationSeconds = static_cast<float>(std::chrono::duration_cast<std::chrono::milliseconds> (mapDuration).count()) / 1000;

   std::ostream nullStream(nullptr);
   auto prt =  std::make_shared<PrintMetrics>(nullStream, m_logStream, 80);
   printMapAlignCommonStats(*prt, *this, "SUMMARY","MAPPING/ALIGNING ", "", "Total input reads",false,mapDurationSeconds);

   m_logStream.flush();
}

//-------------------------------------------------------------------------------
// print out the map/align stats (common to summary and per group stats printing)
//
void ReadGroupAlignmentCounts::printMapAlignCommonStats(
    PrintMetrics&             ma_printer,
    ReadGroupAlignmentCounts& stats,
    const std::string&        runName,
    const std::string&        sectionName,
    const std::string&        sampleName,
    const std::string&        inputReadsDescription,
    const bool                setPrintedSampleName,
    float mapDurationSeconds)
{
  uint64_t readsWithoutMateSequenced = stats.m_numRecords - stats.m_hasMate;
  uint64_t pairedReads     = (stats.m_concordantOnce + stats.m_concordantMultiple + stats.m_pairedDiscordant);
  uint64_t properlyPaired  = stats.m_concordantOnce + stats.m_concordantMultiple;
  uint64_t mapped          = stats.m_numRecords - stats.m_unmapped;
  uint64_t totalAlignments = mapped + stats.m_numSecondary + stats.m_numSupplementary;
  float    estReadLen      = static_cast<float>(stats.m_sumSeqLength) / stats.m_numRecords;

  ma_printer.setSection(sectionName + runName);
  const bool skipJson = false;
  ma_printer.setSampleName(sampleName, skipJson, setPrintedSampleName);
  ma_printer.setDenominator(stats.m_numRecords);
  ma_printer.printCountAndPercentage(inputReadsDescription, stats.m_numRecords);

  ma_printer.printCountAndPercentage("Number of duplicate marked reads", stats.m_numDuplicatesMarked);
 // if (m_removeDuplicates) {
    ma_printer.printCountAndPercentage(
        "Number of duplicate marked and mate reads removed", stats.m_numDuplicatesRemoved);
  //} else {
 //   ma_printer.printNA("Number of duplicate marked and mate reads removed");
//  }
  ma_printer.printCountAndPercentage(
      "Number of unique reads (excl. duplicate marked reads)",
      stats.m_numRecords - stats.m_numDuplicatesMarked);

  ma_printer.printCountAndPercentage("Reads with mate sequenced", stats.m_hasMate);
  ma_printer.printCountAndPercentage("Reads without mate sequenced", readsWithoutMateSequenced);
  ma_printer.printCountAndPercentage("QC-failed reads", stats.m_QCfailed);

  ma_printer.printCountAndPercentage("Mapped reads", mapped);
  ma_printer.printCountAndPercentage(
      "Mapped reads adjusted for filtered mapping", mapped + stats.m_filteredContig + stats.m_nonrefDecoy);

  uint64_t mappedR1 = stats.m_numRecordsR1 - stats.m_unmappedR1;
  uint64_t mappedR2 = stats.m_numRecordsR2 - stats.m_unmappedR2;
  ma_printer.printCountAndPercentage("Mapped reads R1", mappedR1, stats.m_numRecordsR1);
  ma_printer.printCountAndPercentage("Mapped reads R2", mappedR2, stats.m_numRecordsR2);
  ma_printer.printCountAndPercentage(
      "Number of unique & mapped reads (excl. duplicate marked reads)", mapped - stats.m_numDuplicatesMarked);

  ma_printer.printCountAndPercentage("Unmapped reads", stats.m_unmapped);
  ma_printer.printCountAndPercentage(
      "Unmapped reads adjusted for filtered mapping",
      stats.m_unmapped - stats.m_filteredContig - stats.m_nonrefDecoy);
  ma_printer.printCountAndPercentage(
      "Adjustment of reads matching non-reference decoys", stats.m_nonrefDecoy);

  ma_printer.printCountAndPercentage("Singleton reads (itself mapped; mate unmapped)", stats.m_singleton);
  ma_printer.printCountAndPercentage("Paired reads (itself & mate mapped)", pairedReads);
  ma_printer.printCountAndPercentage("Properly paired reads", properlyPaired);
  ma_printer.printCountAndPercentage("Not properly paired reads (discordant)", stats.m_pairedDiscordant);
  ma_printer.printCountAndPercentage(
      "Paired reads mapped to different chromosomes", stats.m_pairedReadDiffContig, pairedReads);
  ma_printer.printCountAndPercentage(
      "Paired reads mapped to different chromosomes (MAPQ>=10)",
      stats.m_pairedReadDiffContigMapq10,
      pairedReads);
  ma_printer.printCountAndPercentage("Reads with MAPQ [40:inf)", stats.m_mapq_hist[4]);
  ma_printer.printCountAndPercentage("Reads with MAPQ [30:40)", stats.m_mapq_hist[3]);
  ma_printer.printCountAndPercentage("Reads with MAPQ [20:30)", stats.m_mapq_hist[2]);
  ma_printer.printCountAndPercentage("Reads with MAPQ [10:20)", stats.m_mapq_hist[1]);
  ma_printer.printCountAndPercentage("Reads with MAPQ [ 0:10)", stats.m_mapq_hist[0]);
  ma_printer.printCountAndPercentage("Reads with MAPQ NA (Unmapped reads)", stats.m_unmapped);

  ma_printer.printCountAndPercentage("Reads with indel R1", stats.m_numIndelReadsR1, mappedR1);
  if (mappedR2 > 0)
    ma_printer.printCountAndPercentage("Reads with indel R2", stats.m_numIndelReadsR2, mappedR2);

  // is this a paired run and should we print R2 metrics?
  bool paired_run = (stats.m_sumSeqLengthR2 > 0);

  ma_printer.printCount("Total bases", stats.m_sumSeqLengthR1 + stats.m_sumSeqLengthR2);
  ma_printer.printCount("Total bases R1", stats.m_sumSeqLengthR1);
  if (paired_run) ma_printer.printCount("Total bases R2", stats.m_sumSeqLengthR2);

  uint64_t sumMappedSeqLengthR1 = stats.m_sumSeqLengthR1 - stats.m_sumUnmappedSeqLengthR1;
  uint64_t sumMappedSeqLengthR2 = stats.m_sumSeqLengthR2 - stats.m_sumUnmappedSeqLengthR2;

  ma_printer.printCount("Mapped bases R1", sumMappedSeqLengthR1);
  if (paired_run) ma_printer.printCount("Mapped bases R2", sumMappedSeqLengthR2);

  ma_printer.printCountAndPercentage("Soft-clipped bases R1", stats.m_numSoftClippedR1, sumMappedSeqLengthR1);
  if (paired_run) {
    ma_printer.printCountAndPercentage(
        "Soft-clipped bases R2", stats.m_numSoftClippedR2, sumMappedSeqLengthR2);
  } else {
    ma_printer.printNA("Soft-clipped bases R2");
  }

  ma_printer.printCountAndPercentage("Mismatched bases R1", stats.m_numMismatchesR1, sumMappedSeqLengthR1);
  if (paired_run) {
    ma_printer.printCountAndPercentage("Mismatched bases R2", stats.m_numMismatchesR2, sumMappedSeqLengthR2);
  } else {
    ma_printer.printNA("Mismatched bases R2");
  }

  if (stats.m_numMismatchesR1 > stats.m_numIndelBasesR1 || stats.m_numIndelBasesR1 == 0) {
    uint64_t mismatchedBasesR1exclIndels = stats.m_numMismatchesR1 - stats.m_numIndelBasesR1;
    ma_printer.printCountAndPercentage(
        "Mismatched bases R1 (excl. indels)", mismatchedBasesR1exclIndels, sumMappedSeqLengthR1);
  } else {
    ma_printer.printNA("Mismatched bases R1 (excl. indels)");
  }

  // if no MN flags in BAM then numMisMatches will be 0
  if (paired_run && (stats.m_numMismatchesR2 > stats.m_numIndelBasesR2 || stats.m_numIndelBasesR2 == 0)) {
    uint64_t mismatchedBasesR2exclIndels = stats.m_numMismatchesR2 - stats.m_numIndelBasesR2;
    ma_printer.printCountAndPercentage(
        "Mismatched bases R2 (excl. indels)", mismatchedBasesR2exclIndels, sumMappedSeqLengthR2);
  } else {
    ma_printer.printNA("Mismatched bases R2 (excl. indels)");
  }

  ma_printer.printCountAndPercentage(
      "Q30 bases",
      stats.m_numAllQ30BasesR1 + stats.m_numAllQ30BasesR2,
      stats.m_sumSeqLengthR1 + stats.m_sumSeqLengthR2);

  ma_printer.printCountAndPercentage("Q30 bases R1", stats.m_numAllQ30BasesR1, stats.m_sumSeqLengthR1);
  if (paired_run)
    ma_printer.printCountAndPercentage("Q30 bases R2", stats.m_numAllQ30BasesR2, stats.m_sumSeqLengthR2);

  //if (m_dupMetricsEnabled)
  //  ma_printer.printCount("Q30 bases (excl. dups & clipped bases)", stats.m_numNonDupNonClippedQ30Bases);

  ma_printer.printCount("Total alignments", totalAlignments);
  ma_printer.printCount("Secondary alignments", stats.m_numSecondary);
  ma_printer.printCount("Supplementary (chimeric) alignments", stats.m_numSupplementary);
  ma_printer.printFloat("Estimated read length", estReadLen);

  if (mapDurationSeconds > 0.005) {
    const float thousand_reads_per_second = mapped / mapDurationSeconds / 1000;
    ma_printer.printFloat("DRAGEN mapping rate [thousand reads/second]", thousand_reads_per_second);
  } else {
    ma_printer.printNA("DRAGEN mapping rate [thousand reads/second]");
  }
}
