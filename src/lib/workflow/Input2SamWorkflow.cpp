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

#include <cerrno>
#include <cstring>
#include <fstream>
#include <limits>

#include "boost/iostreams/filter/gzip.hpp"

#include <boost/filesystem.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/filtering_stream.hpp>

#include "align/Aligner.hpp"
#include "align/SinglePicker.hpp"
#include "bam/BamBlockReader.hpp"
#include "bam/Tokenizer.hpp"
#include "common/Debug.hpp"
#include "common/Threads.hpp"
#include "fastq/FastqBlockReader.hpp"
#include "fastq/Tokenizer.hpp"
#include "io/Bam2ReadTransformer.hpp"
#include "io/Fastq2ReadTransformer.hpp"
#include "mapping_stats.hpp"
#include "options/DragenOsOptions.hpp"
#include "reference/ReferenceDir.hpp"
#include "sam/SamGenerator.hpp"

#include "workflow/DualFastq2SamWorkflow.hpp"
#include "workflow/Input2SamWorkflow.hpp"

#include "workflow/alignment/AlignmentUtils.hpp"

namespace dragenos {
namespace workflow {

std::ostream& operator<<(std::ostream& os, const align::Aligner::ReadPair& pair)
{
  if (2 == pair.size()) {
    return os << "RedPair(" << pair.front() << "-" << pair.back() << ")";
  }

  assert(1 == pair.size());
  return os << "RedPair(" << pair.front();
}

template <typename Tokenizer>
align::InsertSizeParameters requestInsertSizeInfo(
    const options::DragenOsOptions& options,
    align::InsertSizeDistribution&  insertSizeDistribution,
    std::istream&                   input)
{
  Tokenizer                  tokenizer(input);
  align::Aligner::Read::Name lastName;

  align::InsertSizeParameters ret;
  bool                        retDone = false;
  while (tokenizer.next()) {
    if (!insertSizeDistribution.notGoingToBlock()) {
      // Our thread supplied readgroup_insert_stats with last of the input data. Next call to
      // getInsertSizeParameters will block, unless we go and add some alignment results to the
      // readgroup_insert_stats. So, break this fake loop and have everything aligned with dummy insert size
      // parameters so that readgroup_insert_stats gets into the correct state of mind.
      break;
    }
    const auto& token = tokenizer.token();

    const auto& name = token.getName(options.inputQnameSuffixDelim_);
    // only request insert sizes if pair of reads with identical names arrives
    if (align::Aligner::Read::Name(name.first, name.second) == lastName) {
      if (!retDone) {
        ret     = insertSizeDistribution.getInsertSizeParameters(token.readLength());
        retDone = true;
      } else {
        // keep calling this to make sure the insert stats
        // counts match properly
        insertSizeDistribution.getInsertSizeParameters(token.readLength());
      }
    }
    lastName.assign(name.first, name.second);
  }

  return ret;
}

template <typename ReadTransformer>
ReadTransformer makeReadTransformer(const options::DragenOsOptions& options);

template <>
io::FastqToReadTransformer makeReadTransformer<io::FastqToReadTransformer>(
    const options::DragenOsOptions& options)
{
  return io::FastqToReadTransformer(options.inputQnameSuffixDelim_, options.fastqOffset_);
}

template <>
io::BamToReadTransformer makeReadTransformer<io::BamToReadTransformer>(
    const options::DragenOsOptions& options)
{
  return io::BamToReadTransformer(options.inputQnameSuffixDelim_);
}

/**
 * \brief parses single fastq file. Detects interleaved paired records onthe fly.
 *        Should support mixed single/paired data too (not tested)
 */
template <typename ReadTransformer, typename Tokenizer, typename StoreOp>
void alignSingleInput(
    align::InsertSizeParameters&    insertSizeParameters,
    const options::DragenOsOptions& options,
    std::istream&                   input,
    align::Aligner&                 aligner,
    const align::SinglePicker&      singlePicker,
    const align::PairBuilder&       pairBuilder,
    StoreOp                         store)
{
  Tokenizer tokenizer(input);

  align::AlignmentPairs      alignmentPairs;
  align::Aligner::Alignments alignments;

  ReadTransformer          input2Read = makeReadTransformer<ReadTransformer>(options);
  align::Aligner::ReadPair pair;

  align::Aligner::Read::Name lastName;
  uint64_t                   fragmentId = 0;
  while (tokenizer.next()) {
    const auto& token = tokenizer.token();

    const auto& name = token.getName(options.inputQnameSuffixDelim_);
    if (options.interleaved_ && align::Aligner::Read::Name(name.first, name.second) == lastName) {
      // interleaved fastq case, just treat it as paired
      input2Read(token, 1, fragmentId - 1, pair[1]);
      alignment::alignAndStorePair(
          insertSizeParameters, pair, aligner, singlePicker, pairBuilder, alignmentPairs, store);
      lastName.clear();
    } else {
      if (!lastName.empty()) {
        alignment::alignAndStoreSingle(pair.at(0), aligner, singlePicker, alignments, store);
      }
      input2Read(token, 0, fragmentId, pair[0]);
      ++fragmentId;
      lastName.assign(name.first, name.second);
    }
  }

  if (!lastName.empty()) {
    // last unprocesses single-ended read
    alignment::alignAndStoreSingle(pair.at(0), aligner, singlePicker, alignments, store);
  }

  assert(input.eof());
}

template <typename BlockReader>
BlockReader makeBlockReader(std::istream& is, const options::DragenOsOptions& options);

template <>
bam::BamBlockReader makeBlockReader<bam::BamBlockReader>(
    std::istream& is, const options::DragenOsOptions& options)
{
  return bam::BamBlockReader(is, options.inputQnameSuffixDelim_);
}

template <>
fastq::FastqBlockReader makeBlockReader<fastq::FastqBlockReader>(
    std::istream& is, const options::DragenOsOptions& options)
{
  return fastq::FastqBlockReader(is);
}

template <typename ReadTransformer, typename Tokenizer, typename BlockReader>
void parseSingleInput(
    std::istream&                       is,
    std::ostream&                       os,
    const options::DragenOsOptions&     options,
    const reference::ReferenceSequence& refSeq,
    const reference::HashtableConfig&   htConfig,
    const reference::Hashtable&         hashtable,
    std::ostream&                       mappingMetricsLogStream)
{
  std::chrono::system_clock::time_point timeStart = std::chrono::system_clock::now();

  align::InsertSizeDistribution insertSizeDistribution(
      options.samplingEnabled_,
      options.alignerPeq25_,
      options.alignerPeq50_,
      options.alignerPeq75_,
      options.alignerPeMeanInsert_,
      options.alignerPeStddevInsert_,
      options.alignerPeMeanReadLen_,
      options.peStatsIntervalSize_,
      options.peStatsSampleSize_,
      options.peStatsIntervalMemory_,
      options.peStatsIntervalDelay_,
      options.peStatsContinuousUpdate_,
      options.peStatsUpdateLogOnly_,
      options.alignerMapqMax_,
      options.alignerPeOrientation_,
      options.alignerResqueSigmas_,
      options.alignerResqueCeilFactor_,
      options.alignerResqueMinIns_,
      options.alignerResqueMaxIns_,
      std::cerr);

  const align::SimilarityScores similarity(options.matchScore_, options.mismatchScore_);
  align::SinglePicker           singlePicker(
      similarity,
      options.alnMinScore_,
      options.suppMinScoreAdj_,
      options.alignerSecAligns_,
      options.alignerSecScoreDelta_,
      options.alignerSecPhredDelta_,
      options.alignerSecAlignsHard_,
      options.alignerMapqMinLen_,
      options.alignerSampleMapq0_);

  const sam::SamGenerator sam(htConfig);

  // idle threads needed to hold results that arrive out of order
  const int                             poolThreadCount = options.mapperNumThreads_ * 2;
  ReadGroupAlignmentCounts              mappingMetricsGlobal(mappingMetricsLogStream);
  std::vector<ReadGroupAlignmentCounts> mappingMetricsVector(
      poolThreadCount, ReadGroupAlignmentCounts(mappingMetricsLogStream));

  BlockReader reader = makeBlockReader<BlockReader>(is, options);

  std::size_t threadID              = 0;
  std::size_t cpuThreads            = 0;
  int         blockToRead           = 0;
  int         blockToGetInsertSizes = 0;
  int         blockToAlign          = 0;
  int         blockToStore          = options.preserveMapAlignOrder_ ? 0 : -1;
  // let all threads do the job have twice the hardware to make sure there are threads to
  // align while others are stuck in the save queue by one that takes
  // unexpectedly long time
  common::CPU_THREADS(poolThreadCount)
      .execute(
          [&](common::ThreadPool::lock_type& lock) {
            align::PairBuilder pairBuilder(
                similarity,
                options.alnMinScore_,
                options.alignerUnpairedPen_,
                options.alignerXsPairPen_,
                options.alignerSecAligns_,
                options.alignerSecScoreDelta_,
                options.alignerSecPhredDelta_,
                options.alignerSecAlignsHard_,
                options.alignerMapqMinLen_,
                options.alignerSampleMapq0_);

            align::Aligner aligner(
                refSeq,
                htConfig,
                hashtable,
                options.mapOnly_,
                options.swAll_,
                similarity,
                options.gapInitPenalty_,
                options.gapExtendPenalty_,
                options.unclipScore_,
                options.alnMinScore_,
                options.alignerMapqMinLen_,
                options.alignerUnpairedPen_,
                options.mapperFilterLenRatio_,
                !options.methodSmithWaterman_.compare("mengyao"));

            static const std::size_t BUFFER_SIZE = 1024 * 256;

            // records in output format
            std::vector<char> tmpBuffer;
            tmpBuffer.reserve(BUFFER_SIZE * 2);
            boost::iostreams::filtering_ostream ostrm;
            ostrm.push(boost::iostreams::back_insert_device<std::vector<char>>(tmpBuffer));

            //    char              inBuffer[BUFFER_SIZE];
            std::vector<char> inBuffer(BUFFER_SIZE);
            std::vector<char> outBuffer;

            //    auto                      lock                = common::CPU_THREADS().lock();
            ReadGroupAlignmentCounts& mappingMetricsLocal = mappingMetricsVector[threadID];
            threadID++;

            // common::CPU_THREADS().checkThreadFailed() is needed to prevent lucky threads
            // from spending eternity trying to complete the processing of broken data
            while (!common::CPU_THREADS().checkThreadFailed() && !reader.eof()) {
              const std::size_t n = reader.read(&inBuffer[0], BUFFER_SIZE);

              const int ourBlock = blockToRead;
              ++blockToRead;

              while (blockToGetInsertSizes != ourBlock) {
                common::CPU_THREADS().waitForChange(lock);
              }

              align::InsertSizeParameters insertSizeParameters;
              if (options.interleaved_) {
                // sending paired data to readgroup_insert_stats is only allowed if it is treated as paired
                // data. Else, the sent and received counts will mismatch and the whole thing gets stuck
                common::unlock_guard<common::ThreadPool::lock_type> unlock(lock);
                boost::iostreams::filtering_istream                 istrm;
                istrm.push(boost::iostreams::basic_array_source<char>{&inBuffer[0], &inBuffer[0] + n});
                insertSizeParameters =
                    requestInsertSizeInfo<Tokenizer>(options, insertSizeDistribution, istrm);
              }
              assert(blockToGetInsertSizes == ourBlock);
              ++blockToGetInsertSizes;
              common::CPU_THREADS().notify_all();

              while (options.mapperNumThreads_ == int(cpuThreads) || blockToAlign != ourBlock) {
                common::CPU_THREADS().waitForChange(lock);
              }
              ++cpuThreads;
              ++blockToAlign;
              common::CPU_THREADS().notify_all();
              {
                common::unlock_guard<common::ThreadPool::lock_type> unlock(lock);
                //    std::cerr << "read n:" << n << " end: " << std::string(buffer, buffer + n) << std::endl;
                boost::iostreams::filtering_istream istrm;
                istrm.push(boost::iostreams::basic_array_source<char>{&inBuffer[0], &inBuffer[0] + n});
                outBuffer.clear();
                tmpBuffer.clear();

                alignSingleInput<ReadTransformer, Tokenizer>(
                    insertSizeParameters,
                    options,
                    istrm,
                    aligner,
                    singlePicker,
                    pairBuilder,
                    [&](const sequences::Read& r, const align::Alignment& a) {
                      sam.generateRecord(ostrm, r, a, options.rgid_) << "\n";

                      const auto before = outBuffer.size();
                      outBuffer.resize(before + sequences::SerializedRead::getByteSize(r));
                      const auto before2 = outBuffer.size();
                      outBuffer.resize(before2 + align::SerializedAlignment::getByteSize(a));

                      // resize can invalidate references...
                      sequences::SerializedRead& sr =
                          *reinterpret_cast<sequences::SerializedRead*>(&outBuffer.front() + before);
                      sr << r;

                      align::SerializedAlignment& sa =
                          *reinterpret_cast<align::SerializedAlignment*>(&outBuffer.front() + before2);
                      sa << a;

                      mappingMetricsLocal.addRecord(sa, sr);
                    });
                ostrm.flush();
              }

              --cpuThreads;
              common::CPU_THREADS().notify_all();

              if (options.preserveMapAlignOrder_) {
                while (blockToStore != ourBlock) {
                  common::CPU_THREADS().waitForChange(lock);
                }
              } else {
                while (-1 != blockToStore) {
                  common::CPU_THREADS().waitForChange(lock);
                }
                blockToStore = ourBlock;
              }

              {
                common::unlock_guard<common::ThreadPool::lock_type> unlock(lock);
                for (auto it = outBuffer.begin(); outBuffer.end() != it;) {
                  char*                            p = &*it;
                  const sequences::SerializedRead* pRead =
                      reinterpret_cast<const sequences::SerializedRead*>(p);
                  it += pRead->getByteSize();
                  p = &*it;
                  const align::SerializedAlignment* pAlignment =
                      reinterpret_cast<const align::SerializedAlignment*>(p);
                  it += pAlignment->getByteSize();
                  // sam.generateRecord(os, *pRead, *pAlignment, options.rgid_) << "\n";
                  insertSizeDistribution.add(*pAlignment, *pRead);
                }
                if (!os.write(&tmpBuffer.front(), tmpBuffer.size())) {
                  throw std::logic_error(
                      std::string("Error writing output stream. Error: ") + strerror(errno));
                }
              }

              assert(ourBlock == blockToStore);
              if (options.preserveMapAlignOrder_) {
                ++blockToStore;
              } else {
                blockToStore = -1;
              }
              common::CPU_THREADS().notify_all();
            }
          },
          options.mapperNumThreads_);

  // aggregate and print mapping metrics
  for (std::size_t ii = 0; ii < mappingMetricsVector.size(); ii++) {
    mappingMetricsGlobal.add(mappingMetricsVector[ii]);
  }

  mappingMetricsGlobal.printStats(std::chrono::system_clock::now() - timeStart);

  insertSizeDistribution.forceInitDoneSending();
  if (options.interleaved_) {
    std::cerr << insertSizeDistribution << std::endl;
  }
}

bool isBam(const std::string& fileName)
{
  return fileName.length() > 4 && ".bam" == fileName.substr(fileName.length() - 4);
}

bool isGzip(const std::string& fileName)
{
  return (fileName.length() > 3 && ".gz" == fileName.substr(fileName.length() - 3));
}

void parseSingleInput(
    std::ostream&                       os,
    const options::DragenOsOptions&     options,
    const reference::ReferenceSequence& refSeq,
    const reference::HashtableConfig&   htConfig,
    const reference::Hashtable&         hashtable,
    std::ostream&                       mappingMetricsLogStream)
{
  std::cerr << "Running fastq workflow on " << options.mapperNumThreads_ << " threads. System supports "
            << std::thread::hardware_concurrency() << " threads." << std::endl;

  std::ifstream                       file(options.inputFile1_, std::ios_base::in | std::ios_base::binary);
  boost::iostreams::filtering_istream input;
  if (isGzip(options.inputFile1_) || isBam(options.inputFile1_)) {
    input.push(boost::iostreams::gzip_decompressor());
  }

  input.push(file);
  input.exceptions(std::ios_base::badbit);
  try {
    if (isBam(options.inputFile1_)) {
      parseSingleInput<io::BamToReadTransformer, bam::Tokenizer, bam::BamBlockReader>(
          input, os, options, refSeq, htConfig, hashtable, mappingMetricsLogStream);
    } else {
      parseSingleInput<io::FastqToReadTransformer, fastq::Tokenizer, fastq::FastqBlockReader>(
          input, os, options, refSeq, htConfig, hashtable, mappingMetricsLogStream);
    }
  } catch (boost::iostreams::gzip_error& e) {
    BOOST_THROW_EXCEPTION(std::runtime_error(
        e.what() + std::string(" ") + std::to_string(e.error()) +
        (boost::iostreams::gzip::zlib_error == e.error()
             ? (" zlib error:" + std::to_string(e.zlib_error_code()))
             : std::string(""))));
  }
}

void input2Sam(const dragenos::options::DragenOsOptions& options)
{
  if (options.buildHashTable_ || options.htUncompress_) {
    return;
  }

  DRAGEN_OS_THREAD_CERR << "Version: " << common::Version::string() << std::endl;
  DRAGEN_OS_THREAD_CERR << "argc: " << options.argc() << " argv: " << options.getCommandLine() << std::endl;

  const reference::ReferenceDir7 referenceDir(
      options.refDir_, options.mmapReference_, options.loadReference_);

  /**
   ** \brief memory mapped hashtable data
   **
   ** Note: the custom destructor would be munmap, but munmap needs to know
   ** the size of the memory segment (hashtableConfig_.getHashtableBytes())
   ** which requires using a lambda as a constructor, which in turn requires
   ** declaring the type of the destructor as "std::function<void(void*)>>"
   ** instead of using the type of the function pointer "void(*)(uint64_t*)"
   **
   ** Note: the type can't be const because of munmap signature.
   **/
  //const std::unique_ptr<uint64_t, std::function<void(uint64_t*)>> hashtableData_;
  const reference::Hashtable hashtable(
      &referenceDir.getHashtableConfig(), referenceDir.getHashtableData(), referenceDir.getExtendTableData());

  std::ofstream os;
  namespace bfs = boost::filesystem;
  if (!options.outputDirectory_.empty()) {
    if (!exists(bfs::path(options.outputDirectory_))) {
      BOOST_THROW_EXCEPTION(common::IoException(
          ENOENT, std::string("Output directory does not exist: ") + options.outputDirectory_));
    }
    const auto filePath = bfs::path(options.outputDirectory_) / (options.outputFilePrefix_ + ".sam");
    os.open(filePath.c_str());
    if (!os) {
      BOOST_THROW_EXCEPTION(common::IoException(
          errno, std::string("Failed to create SAM file: ") + filePath.string() + ": " + strerror(errno)));
    }
    if (options.verbose_) {
      std::cerr << "INFO: writing SAM file to " << filePath << std::endl;
    }
  }
  std::ostream& samFile = os.is_open() ? os : std::cout;
  sam::SamGenerator::generateHeader(
      samFile, referenceDir.getHashtableConfig(), options.getCommandLine(), options.rgid_, options.rgsm_);

  std::ofstream mappingMetricsLogStream;

  if (!options.outputDirectory_.empty()) {
    namespace bfs = boost::filesystem;
    if (!exists(bfs::path(options.outputDirectory_))) {
      BOOST_THROW_EXCEPTION(common::IoException(
          ENOENT, std::string("Output directory does not exist: ") + options.outputDirectory_));
    }

    const auto filePathMetrics =
        bfs::path(options.outputDirectory_) / (options.outputFilePrefix_ + ".mapping_metrics.csv");
    mappingMetricsLogStream.open(filePathMetrics.c_str());
    if (!mappingMetricsLogStream) {
      BOOST_THROW_EXCEPTION(common::IoException(
          errno,
          std::string("Failed to create mapping metrics file: ") + filePathMetrics.string() + ": " +
              strerror(errno)));
    }

    if (options.verbose_) {
      std::cerr << "INFO: writing mapping metrics stats into " << filePathMetrics << std::endl;
    }
  }

  if (options.inputFile2_.empty()) {
    parseSingleInput(
        samFile,
        options,
        referenceDir.getReferenceSequence(),
        referenceDir.getHashtableConfig(),
        hashtable,
        mappingMetricsLogStream.is_open() ? mappingMetricsLogStream : std::cerr);
  } else {
    DualFastq2SamWorkflow workflow(
        options, referenceDir.getReferenceSequence(), referenceDir.getHashtableConfig(), hashtable);
    std::ofstream insertSizeDistributionLogStream;

    if (!options.outputDirectory_.empty()) {
      namespace bfs = boost::filesystem;
      if (!exists(bfs::path(options.outputDirectory_))) {
        BOOST_THROW_EXCEPTION(common::IoException(
            ENOENT, std::string("Output directory does not exist: ") + options.outputDirectory_));
      }
      const auto filePath =
          bfs::path(options.outputDirectory_) / (options.outputFilePrefix_ + ".insert-stats.tab");
      insertSizeDistributionLogStream.open(filePath.c_str());
      if (!insertSizeDistributionLogStream) {
        BOOST_THROW_EXCEPTION(common::IoException(
            errno,
            std::string("Failed to create insert-stats file: ") + filePath.string() + ": " +
                strerror(errno)));
      }

      if (options.verbose_) {
        std::cerr << "INFO: writing insert stats into " << filePath << std::endl;
      }
      insertSizeDistributionLogStream
          << "RG\tID\tQ25\tQ50\tQ75\tMEAN\tSTD\tMIN-INS\tMAX-INS\tRESCUE-MIN-INS\tRESCUE-MAX-INS\tNPAIRS\tINTSIZE\tSTALLS\tDELAY"
          << std::endl;
    }
    workflow.parseDualFastq(
        samFile,
        insertSizeDistributionLogStream.is_open() ? insertSizeDistributionLogStream : std::cerr,
        mappingMetricsLogStream.is_open() ? mappingMetricsLogStream : std::cerr);
  }
}

}  // namespace workflow
}  // namespace dragenos
