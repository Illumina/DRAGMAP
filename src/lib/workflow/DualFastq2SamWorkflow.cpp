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

#include <fstream>
#include <limits>

#include <boost/iostreams/device/back_inserter.hpp>
// #include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>

#include "common/Debug.hpp"
#include "common/Threads.hpp"
#include "mapping_stats.hpp"

#include "align/Aligner.hpp"
#include "align/SinglePicker.hpp"
#include "fastq/Tokenizer.hpp"
#include "io/Fastq2ReadTransformer.hpp"
#include "sam/SamGenerator.hpp"

#include "workflow/DualFastq2SamWorkflow.hpp"

#include "workflow/alignment/AlignmentUtils.hpp"

namespace dragenos {
namespace workflow {

align::InsertSizeParameters DualFastq2SamWorkflow::requestInsertSizeInfo(
    align::InsertSizeDistribution& insertSizeDistribution, std::istream& inputR1, std::istream& inputR2)
{
  fastq::Tokenizer r1Tokenizer(inputR1);
  fastq::Tokenizer r2Tokenizer(inputR2);

  align::InsertSizeParameters ret;
  bool                        retDone = false;
  while (true) {
    // Always call next() on both tokenizers, to make sure both files are read to the end.
    const bool r1Next = r1Tokenizer.next();
    const bool r2Next = r2Tokenizer.next();
    if (!r1Next || !r2Next) {
      break;
    }

    const auto& r1Token = r1Tokenizer.token();
    const auto& r2Token = r2Tokenizer.token();
    //    std::cout << "token: " << r1Token << "-" << r2Token << "\n";
    assert(r1Token.valid() && r2Token.valid());

    if (!retDone) {
      ret     = insertSizeDistribution.getInsertSizeParameters(r1Token.readLength());
      retDone = true;
    } else {
      // keep calling this to make sure the insert stats
      // counts match properly
      insertSizeDistribution.getInsertSizeParameters(r1Token.readLength());
    }
  }

  // make sure both files have been read to end
  assert(inputR1.eof());
  assert(inputR2.eof());
  // make sure there is no case of one file having a good read and the other one not
  assert(!r1Tokenizer.token().valid() && !r2Tokenizer.next());

  return ret;
}

template <typename StoreOp>
void DualFastq2SamWorkflow::alignDualFastq(
    align::InsertSizeParameters& insertSizeParameters,
    std::istream&                inputR1,
    std::istream&                inputR2,
    align::Aligner&              aligner,
    const align::SinglePicker&   singlePicker,
    const align::PairBuilder&    pairBuilder,
    StoreOp                      store)
{
  fastq::Tokenizer r1Tokenizer(inputR1);
  fastq::Tokenizer r2Tokenizer(inputR2);

  align::AlignmentPairs alignmentPairs;

  io::FastqToReadTransformer fastq2Read(options_.inputQnameSuffixDelim_, options_.fastqOffset_);
  align::Aligner::ReadPair   pair;

  int64_t fragmentId = 0;
  while (true) {
    // Always call next() on both tokenizers, to make sure both files are read to the end.
    const bool r1Next = r1Tokenizer.next();
    const bool r2Next = r2Tokenizer.next();
    if (!r1Next || !r2Next) {
      break;
    }

    const auto& r1Token = r1Tokenizer.token();
    const auto& r2Token = r2Tokenizer.token();
    //    std::cout << "token: " << r1Token << "-" << r2Token << "\n";
    assert(r1Token.valid() && r2Token.valid());

    //    fragment.at(0) = fastq2Read(r1Token, 0, fragment.size());
    //    fragment.at(1) = fastq2Read(r2Token, 1, fragment.size());
    fastq2Read(r1Token, 0, fragmentId, pair.at(0));
    //    std::cerr << "r1:" << pair[0] << std::endl;
    fastq2Read(r2Token, 1, fragmentId, pair.at(1));
    //    std::cerr << "r2:" << pair[1] << std::endl;
    //    std::cout << "fragment: " << fragment << "\n";
    //    std::cout << "tada: " << fastq2Read.tmpName_.capacity() << std::endl;

    alignmentPairs.clear();
    alignment::alignAndStorePair(
        insertSizeParameters, pair, aligner, singlePicker, pairBuilder, alignmentPairs, store);
    ++fragmentId;
  }

  // make sure both files have been read to end
  assert(inputR1.eof());
  assert(inputR2.eof());
  // make sure there is no case of one file having a good read and the other one not
  assert(!r1Tokenizer.token().valid() && !r2Tokenizer.next());
}

void DualFastq2SamWorkflow::readBlockThread(
    common::ThreadPool::lock_type& lock,
    const int                      ourBlock,
    int&                           r1Records,
    std::vector<char>&             r1Block,
    fastq::FastqNRecordReader&     r1Reader,
    int&                           r2Records,
    std::vector<char>&             r2Block,
    fastq::FastqNRecordReader&     r2Reader)
{
  // if a thread is late to the party, both read blocks have been already done.
  // use < instead of != for wait
  while (!r1Eof_ && !r2Eof_ && blockToRead_ < ourBlock) {
    common::CPU_THREADS().waitForChange(lock);
  }

  if (blockToRead_ == ourBlock) {
    // in case reading takes time and more than one thread enters, reading will happen in parallel
    if (!r1Eof_ && -1 == r1Records) {
      r1Records = 0;  // we will be reading this block.
      {
        common::unlock_guard<common::ThreadPool::lock_type> unlock(lock);
        r1Block.clear();
        r1Records = r1Reader.read(std::back_inserter(r1Block), RECORDS_AT_A_TIME_);
        r1Eof_    = r1Reader.eof();
      }
    }

    if (!r2Eof_ &&
        // if r2Recrods != -1 that means another thread is already reading it or done reading, just return.
        -1 == r2Records) {
      r2Records = 0;
      {
        common::unlock_guard<common::ThreadPool::lock_type> unlock(lock);
        r2Block.clear();
        r2Records = r2Reader.read(std::back_inserter(r2Block), RECORDS_AT_A_TIME_);
        r2Eof_    = r2Reader.eof();
      }
    }
  }
  // else, we're supposed to read a block that will never be read because we've run out of input data
}

void DualFastq2SamWorkflow::parseDualFastq(
    std::ostream& os, std::ostream& insertSizeDistributionLogStream, std::ostream& mappingMetricsLogStream)
{
  std::ifstream r1Stream(options_.inputFile1_, std::ios_base::in | std::ios_base::binary);
  //  r1Stream.exceptions(std::ios_base::badbit | std::ios_base::failbit);
  boost::iostreams::filtering_istream r1Decomp;
  if (options_.inputFile1_.length() > 3 &&
      ".gz" == options_.inputFile1_.substr(options_.inputFile1_.length() - 3)) {
    r1Decomp.push(boost::iostreams::gzip_decompressor());
  }

  boost::iostreams::filtering_istream r2Decomp;
  std::ifstream r2Stream(options_.inputFile2_, std::ios_base::in | std::ios_base::binary);
  //  r2Stream.exceptions(std::ios_base::badbit | std::ios_base::failbit);
  if (options_.inputFile2_.length() > 3 &&
      ".gz" == options_.inputFile2_.substr(options_.inputFile2_.length() - 3)) {
    r2Decomp.push(boost::iostreams::gzip_decompressor());
  }

  //    r1Decomp.push(boost::iostreams::file_source(options_.fastqFile1_, BOOST_IOS::in | BOOST_IOS::binary));
  r1Decomp.push(r1Stream);
  r1Decomp.exceptions(std::ios_base::badbit);
  //    r2Decomp.push(boost::iostreams::file_source(options_.fastqFile2_, BOOST_IOS::in | BOOST_IOS::binary));
  r2Decomp.push(r2Stream);
  r2Decomp.exceptions(std::ios_base::badbit);
  try {
    parseDualFastq(r1Decomp, r2Decomp, os, insertSizeDistributionLogStream, mappingMetricsLogStream);
  } catch (boost::iostreams::gzip_error& e) {
    BOOST_THROW_EXCEPTION(std::runtime_error(
        e.what() + std::string(" ") + std::to_string(e.error()) +
        (boost::iostreams::gzip::zlib_error == e.error()
             ? (" zlib error:" + std::to_string(e.zlib_error_code()))
             : std::string(""))));
  }
}

void DualFastq2SamWorkflow::alignDualFastqBlock(
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
    const sam::SamGenerator&               sam)
{
  align::PairBuilder pairBuilder(
      similarity,
      options_.alnMinScore_,
      options_.alignerUnpairedPen_,
      options_.alignerXsPairPen_,
      options_.alignerSecAligns_,
      options_.alignerSecScoreDelta_,
      options_.alignerSecPhredDelta_,
      options_.alignerSecAlignsHard_,
      options_.alignerMapqMinLen_,
      options_.alignerSampleMapq0_);

  // aligner is not stateless, make sure each thread uses its own.
  align::Aligner aligner(
      refSeq_,
      htConfig_,
      hashtable_,
      options_.mapOnly_,
      options_.swAll_,
      similarity,
      options_.gapInitPenalty_,
      options_.gapExtendPenalty_,
      options_.unclipScore_,
      options_.alnMinScore_,
      options_.alignerMapqMinLen_,
      options_.alignerUnpairedPen_,
      options_.mapperFilterLenRatio_,
      !options_.methodSmithWaterman_.compare("mengyao"));

  std::vector<char> r1Block;
  // arbitrary preallocation to avoid unnecessary copy/paste
  r1Block.reserve(RECORDS_AT_A_TIME_ * 1024);
  std::vector<char> r2Block;
  r2Block.reserve(RECORDS_AT_A_TIME_ * 1024);

  // records in output format
  std::vector<char> tmpBuffer;
  tmpBuffer.reserve(RECORDS_AT_A_TIME_ * 1024);
  boost::iostreams::filtering_ostream ostrm;
  ostrm.push(boost::iostreams::back_insert_device<std::vector<char>>(tmpBuffer));
  // minimum data required for insert size calculation
  std::vector<char> insBuffer;
  insBuffer.reserve(RECORDS_AT_A_TIME_ * 1024);

  ReadGroupAlignmentCounts& mappingMetricsLocal = mappingMetricsVector[threadID];
  threadID++;

  do {
    const int ourBlock = blockToStart_++;

    // keep track of which read is reading
    int r1Records = -1;
    int r2Records = -1;
    // no more than 2 threads are allowed in
    common::CPU_THREADS().execute(
        lock,
        [&](common::ThreadPool::lock_type& lock) {
          readBlockThread(lock, ourBlock, r1Records, r1Block, r1Reader, r2Records, r2Block, r2Reader);
        },
        2);

    //      readBlockThread(lock, ourBlock, r1Records, r1Block, r1Reader, r2Records, r2Block, r2Reader);
    // if our block, go see if we've read something
    if (blockToRead_ == ourBlock) {
      // could be -1 or 0 if by the time we decide to read our block we run out of data
      if (0 < r1Records) {
        if (r1Records != r2Records) {
          throw std::logic_error(std::string("fastq files have different number of records "));
        }

        // both reads completed, let others move on
        ++blockToRead_;
        common::CPU_THREADS().notify_all();

        while (blockToGetInsertSizes_ != ourBlock) {
          common::CPU_THREADS().waitForChange(lock);
        }

        align::InsertSizeParameters insertSizeParameters;
        {
          common::unlock_guard<common::ThreadPool::lock_type> unlock(lock);
          boost::iostreams::filtering_istream                 inputR1;
          inputR1.push(boost::iostreams::basic_array_source<char>{&r1Block.front(),
                                                                  &r1Block.front() + r1Block.size()});
          boost::iostreams::filtering_istream inputR2;
          inputR2.push(boost::iostreams::basic_array_source<char>{&r2Block.front(),
                                                                  &r2Block.front() + r2Block.size()});
          insertSizeParameters = requestInsertSizeInfo(insertSizeDistribution, inputR1, inputR2);
        }
        assert(blockToGetInsertSizes_ == ourBlock);
        ++blockToGetInsertSizes_;
        common::CPU_THREADS().notify_all();

        while (options_.mapperNumThreads_ == int(cpuThreads) || blockToAlign_ != ourBlock) {
          common::CPU_THREADS().waitForChange(lock);
        }
        ++cpuThreads;
        ++blockToAlign_;
        common::CPU_THREADS().notify_all();
        {
          common::unlock_guard<common::ThreadPool::lock_type> unlock(lock);
          boost::iostreams::filtering_istream                 inputR1;
          inputR1.push(boost::iostreams::basic_array_source<char>{&r1Block.front(),
                                                                  &r1Block.front() + r1Block.size()});
          boost::iostreams::filtering_istream inputR2;
          inputR2.push(boost::iostreams::basic_array_source<char>{&r2Block.front(),
                                                                  &r2Block.front() + r2Block.size()});
          insBuffer.clear();
          tmpBuffer.clear();

          alignDualFastq(
              insertSizeParameters,
              inputR1,
              inputR2,
              aligner,
              singlePicker,
              pairBuilder,
              [&](const sequences::Read& r, const align::Alignment& a) {
                sam.generateRecord(ostrm, r, a, options_.rgid_) << "\n";

                const auto before = insBuffer.size();
                insBuffer.resize(before + sequences::SerializedRead::getByteSize(r));
                const auto before2 = insBuffer.size();
                insBuffer.resize(before2 + align::SerializedAlignment::getByteSize(a));

                // resize can invalidate references...
                sequences::SerializedRead& sr =
                    *reinterpret_cast<sequences::SerializedRead*>(&insBuffer.front() + before);
                sr << r;

                align::SerializedAlignment& sa =
                    *reinterpret_cast<align::SerializedAlignment*>(&insBuffer.front() + before2);
                sa << a;

                mappingMetricsLocal.addRecord(sa, sr);
              });
          ostrm.flush();
        }

        --cpuThreads;
        common::CPU_THREADS().notify_all();

        if (options_.preserveMapAlignOrder_) {
          while (blockToStore_ != ourBlock) {
            common::CPU_THREADS().waitForChange(lock);
          }
        } else {
          while (-1 != blockToStore_) {
            common::CPU_THREADS().waitForChange(lock);
          }
          blockToStore_ = ourBlock;
        }

        {
          common::unlock_guard<common::ThreadPool::lock_type> unlock(lock);
          for (auto it = insBuffer.begin(); insBuffer.end() != it;) {
            char*                            p     = &*it;
            const sequences::SerializedRead* pRead = reinterpret_cast<const sequences::SerializedRead*>(p);
            it += pRead->getByteSize();
            p = &*it;
            const align::SerializedAlignment* pAlignment =
                reinterpret_cast<const align::SerializedAlignment*>(p);
            it += pAlignment->getByteSize();
            // sam.generateRecord(os, *pRead, *pAlignment, options_.rgid_) << "\n";
            insertSizeDistribution.add(*pAlignment, *pRead);
          }
          if (!os.write(&tmpBuffer.front(), tmpBuffer.size())) {
            throw std::logic_error(std::string("Error writing output stream. Error: ") + strerror(errno));
          }
        }
        assert(blockToStore_ == ourBlock);
        if (options_.preserveMapAlignOrder_) {
          ++blockToStore_;
        } else {
          blockToStore_ = -1;
        }
        common::CPU_THREADS().notify_all();
      }
    }
    // else we're a thread that waited to read its block until it learned that there will be no more
    // input data. just quietly exit
  } while (!common::CPU_THREADS().checkThreadFailed() && !r1Eof_ && !r2Eof_);
}

void DualFastq2SamWorkflow::parseDualFastq(
    std::istream& r1Stream,
    std::istream& r2Stream,
    std::ostream& os,
    std::ostream& insertSizeDistributionLogStream,
    std::ostream& mappingMetricsLogStream)
{
  std::chrono::system_clock::time_point timeStart = std::chrono::system_clock::now();

  std::cerr << "Running dual fastq workflow on " << options_.mapperNumThreads_ << " threads. System supports "
            << std::thread::hardware_concurrency() << " threads." << std::endl;

  align::InsertSizeDistribution insertSizeDistribution(
      options_.samplingEnabled_,
      options_.alignerPeq25_,
      options_.alignerPeq50_,
      options_.alignerPeq75_,
      options_.alignerPeMeanInsert_,
      options_.alignerPeStddevInsert_,
      options_.alignerPeMeanReadLen_,
      options_.peStatsIntervalSize_,
      options_.peStatsSampleSize_,
      options_.peStatsIntervalMemory_,
      options_.peStatsIntervalDelay_,
      options_.peStatsContinuousUpdate_,
      options_.peStatsUpdateLogOnly_,
      options_.alignerMapqMax_,
      options_.alignerPeOrientation_,
      options_.alignerResqueSigmas_,
      options_.alignerResqueCeilFactor_,
      options_.alignerResqueMinIns_,
      options_.alignerResqueMaxIns_,
      insertSizeDistributionLogStream);

  // idle threads needed to hold results that arrive out of order
  const int                             poolThreadCount = options_.mapperNumThreads_ * 2;
  ReadGroupAlignmentCounts              mappingMetricsGlobal(mappingMetricsLogStream);
  std::vector<ReadGroupAlignmentCounts> mappingMetricsVector(
      poolThreadCount, ReadGroupAlignmentCounts(mappingMetricsLogStream));

  const align::SimilarityScores similarity(options_.matchScore_, options_.mismatchScore_);
  align::SinglePicker           singlePicker(
      similarity,
      options_.alnMinScore_,
      options_.suppMinScoreAdj_,
      options_.alignerSecAligns_,
      options_.alignerSecScoreDelta_,
      options_.alignerSecPhredDelta_,
      options_.alignerSecAlignsHard_,
      options_.alignerMapqMinLen_,
      options_.alignerSampleMapq0_);

  fastq::FastqNRecordReader r1Reader(r1Stream);
  fastq::FastqNRecordReader r2Reader(r2Stream);

  const sam::SamGenerator sam(htConfig_);

  std::size_t cpuThreads = 0;
  std::size_t threadID   = 0;
  blockToStore_          = options_.preserveMapAlignOrder_ ? 0 : -1;
  // let all threads do the job have twice the hardware to make sure there are threads to
  // align while others are stuck in the save queue by one that takes
  // unexpectedly long time
  common::CPU_THREADS(poolThreadCount)
      .execute(
          [&](common::ThreadPool::lock_type& lock) {
            alignDualFastqBlock(
                lock,
                cpuThreads,
                threadID,
                mappingMetricsVector,
                insertSizeDistribution,
                r1Reader,
                r2Reader,
                os,
                singlePicker,
                similarity,
                sam);
          },
          options_.mapperNumThreads_);

  // aggregate and print mapping metrics
  for (std::size_t ii = 0; ii < mappingMetricsVector.size(); ii++) {
    mappingMetricsGlobal.add(mappingMetricsVector[ii]);
  }

  mappingMetricsGlobal.printStats(std::chrono::system_clock::now() - timeStart);

  insertSizeDistribution.forceInitDoneSending();
  std::cerr << insertSizeDistribution << std::endl;
}

}  // namespace workflow
}  // namespace dragenos
