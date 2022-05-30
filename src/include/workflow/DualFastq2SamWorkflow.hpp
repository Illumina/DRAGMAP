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
#include "fastq/FastqNRecordReader.hpp"
#include "options/DragenOsOptions.hpp"
#include "reference/Hashtable.hpp"
#include "reference/ReferenceDir.hpp"

namespace dragenos {
namespace workflow {

class DualFastq2SamWorkflow {
  const options::DragenOsOptions&     options_;
  const reference::ReferenceSequence& refSeq_;
  const reference::HashtableConfig&   htConfig_;
  const reference::Hashtable&         hashtable_;
  // IMPORTANT: this has to divide INIT_INTERVAL_SIZE without remainder. Else the whole insert
  // size stats detection will hang because it depends on processing alignment results exactly
  // after sending INIT_INTERVAL_SIZE into the aligner.
  static const int RECORDS_AT_A_TIME_ = 100000;

  // ensure FIFO
  int blockToStart_          = 0;
  int blockToGetInsertSizes_ = 0;
  int blockToAlign_          = 0;
  int blockToRead_           = 0;
  int blockToStore_          = 0;

  bool r1Eof_ = false;
  bool r2Eof_ = false;

public:
  DualFastq2SamWorkflow(
      const options::DragenOsOptions&     options,
      const reference::ReferenceSequence& refSeq,
      const reference::HashtableConfig&   htConfig,
      const reference::Hashtable&         hashtable)
    : options_(options), refSeq_(refSeq), htConfig_(htConfig), hashtable_(hashtable)
  {
  }

  void parseDualFastq(
      std::ostream& os, std::ostream& insertSizeDistributionLogStream, std::ostream& mappingMetricsLogStream);

private:
  align::InsertSizeParameters requestInsertSizeInfo(
      align::InsertSizeDistribution& insertSizeDistribution, std::istream& inputR1, std::istream& inputR2);

  void alignDualFastqBlock(
      common::ThreadPool::lock_type&         lock,
      std::size_t&                           cpuThreads,
      std::size_t&                           threadID,
      std::vector<ReadGroupAlignmentCounts>& mappingMetricsVector,
      align::InsertSizeDistribution&         insertSizeDistribution,
      fastq::FastqNRecordReader&             r1Reader,
      fastq::FastqNRecordReader&             r2Reader,
      std::ostream&                          os,
      const align::SinglePicker&             singlePicker,
      const align::SimilarityScores&         similarity,
      const sam::SamGenerator&               sam);

  template <typename StoreOp>
  void alignDualFastq(
      align::InsertSizeParameters& insertSizeParameters,
      std::istream&                inputR1,
      std::istream&                inputR2,
      align::Aligner&              aligner,
      const align::SinglePicker&   singlePicker,
      const align::PairBuilder&    pairBuilder,
      StoreOp                      store);

  //  void parseDualFastq(
  //    const align::InsertSizeDistribution& insertSizeDistribution,
  //    std::istream& inputR1, std::istream& inputR2, align::Aligner& aligner, std::ostream& output);
  void readBlockThread(
      common::ThreadPool::lock_type& lock,
      const int                      ourBlock,
      int&                           r1Records,
      std::vector<char>&             r1Block,
      fastq::FastqNRecordReader&     r1Reader,
      int&                           r2Records,
      std::vector<char>&             r2Block,
      fastq::FastqNRecordReader&     r2Reader);

  void parseDualFastq(
      std::istream& r1Stream,
      std::istream& r2Stream,
      std::ostream& os,
      std::ostream& insertSizeDistributionLogStream,
      std::ostream& mappingMetricsLogStream);
};

}  // namespace workflow
}  // namespace dragenos
