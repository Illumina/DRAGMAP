#include "gtest/gtest.h"

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <cstring>
#include <cerrno>
#include <stdexcept>
#include <memory>

#include <boost/filesystem.hpp>

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>

#include "sequences/Read.hpp"
#include "sequences/Seed.hpp"
#include "reference/Hashtable.hpp"
#include "map/Mapper.hpp"
#include "reference/ReferenceSequence.hpp"
#include "common/Exceptions.hpp"


// generate refSeq as it is done in hash_table.c
// encoding 4 bases per byte - ACGT -> 0123
//std::vector<uint8_t> generateRefSeq(std::vector<char> bases);
std::vector<uint8_t> generateRefSeq(std::string bases);

template <int START, int COUNT>
uint64_t getBits(uint64_t v)
{
  if (START) v = (v >> START);
  constexpr uint64_t MASK = ~((~uint64_t(0)) << (COUNT));
  return v & MASK;
}

class MapperFixture: public ::testing::Test
{
public:
  MapperFixture()
  {
  }
  // The name was changed in gtest 1.8:   static void SetUpTestSuite()
  static void SetUpTestCase();
  static void TearDownTestCase()
  {
  }
  typedef dragenos::sequences::Read Read;
  typedef Read::Name Name;
  typedef Read::Bases Bases;
  typedef Read::Qualities Qualities;
  Bases getBasesFromReference(const size_t referencePosition, const unsigned readLength) const;
  std::string basesToString(const Bases bases) const
  {
    std::string result;
    for (auto b: bases) result.push_back(b);
    return result;
  }
  std::string bases2bpbToString(const Bases bases) const
  {
    static const char tr[] = {'A', 'C', 'G', 'T'};
    std::string result;
    for (auto b: bases) result.push_back(b > 3 ? 'N' : tr[b]);
    return result;
  }
protected:
  typedef dragenos::reference::HashtableConfig HashtableConfig;
  typedef dragenos::reference::HashtableTraits HashtableTraits;
  typedef dragenos::reference::Hashtable Hashtable;
  typedef dragenos::map::Mapper Mapper;
  typedef dragenos::reference::ReferenceSequence ReferenceSequence;
  static HashtableConfig *hashtableConfig;
  static Hashtable *hashtable;
  static ReferenceSequence referenceSequence;
  static std::unique_ptr<Mapper> mapper;
  static Read read;
  // static constexpr size_t READ_REFERENCE_POSITION = 9990016; // in chr1
  static constexpr size_t READ_REFERENCE_POSITION = 300000000; // in chr2
};


class Environment : public ::testing::Environment {
public:
  typedef dragenos::reference::HashtableConfig HashtableConfig;
  typedef dragenos::reference::HashtableTraits HashtableTraits;
  typedef dragenos::reference::Hashtable Hashtable;
  typedef dragenos::reference::ReferenceSequence ReferenceSequence;
  Environment() : hashtableFd(-1), extendTableFd(-1), table(nullptr), extendTable(nullptr)
  {
  }
  virtual ~Environment()
  {
    if (-1 != hashtableFd) close(hashtableFd);
    if ((nullptr != table) && (MAP_FAILED != table))
    {
      const auto tableSize = hashtableConfig->getHashtableBytes();
      munmap(table, tableSize);
      table = nullptr;
    }
  }
  static unsigned char *getReferenceData() {return referenceData;}
  static size_t getReferenceSize() {return referenceFileSize;}
  static Hashtable *getHashtable() {return hashtable.get();}
  static HashtableConfig *getHashtableConfig() {return hashtableConfig.get();}
  static const std::vector<char> &getHashtableConfigText() {return hashtableConfigText;}
  // Override this to define how to set up the environment.
  void SetUp() override;

  // Override this to define how to tear down the environment.
  void TearDown() override;
  static const std::vector<uint8_t> sequence;
  static const std::vector<uint8_t> qualities;
private:
  static std::vector<char> hashtableConfigText;
  static std::unique_ptr<HashtableConfig> hashtableConfig;
  static std::unique_ptr<Hashtable> hashtable;
  int hashtableFd;
  int extendTableFd;
public:
  uint64_t* table;
  uint64_t* extendTable;
private:
  int referenceFd;
  static size_t referenceFileSize;
  static unsigned char *referenceData;
};

TEST_F(MapperFixture, explore)
{
  using dragenos::sequences::Seed;
  using dragenos::reference::HashRecord;
  using dragenos::reference::ExtendTableInterval;
  using dragenos::map::ChainBuilder;
  const size_t readLength = 151;
  const size_t seedLength = hashtableConfig->getPrimarySeedBases();
  const auto SEED_PERIOD = Seed::DEFAULT_PERIOD;
  const auto SEED_PATTERN = Seed::DEFAULT_PATTERN;
  const auto FORCE_LAST_N_SEEDS = Seed::DEFAULT_FORCE_LAST_N;
  const Bases bases = getBasesFromReference(READ_REFERENCE_POSITION, readLength);
  const Qualities qualities(bases.size(), 30);
  const std::string nameString("name");
  const Read::Name name(nameString.begin(), nameString.end());
  Read read;
  read.init(Name(name), Bases(bases), Qualities(qualities), 0, 0);
  std::cerr << "\n    -------------------\n        position: " << READ_REFERENCE_POSITION << " (" << std::hex << READ_REFERENCE_POSITION << std::dec << "): " << basesToString(bases) << std::endl;
  //std::cerr <<  "bases copy: " << basesToString(Bases(bases)) << std::endl;
  //std::cerr <<  "read: " << bases2bpbToString(read.getBases()) << std::endl;
  //std::cerr << "Generating seed offsets for " << readLength << " " << seedLength << " " << SEED_PERIOD << " " << SEED_PATTERN << " " << FORCE_LAST_N_SEEDS << std::endl;
  const auto seedOffsets = Seed::getSeedOffsets(readLength, seedLength, SEED_PERIOD, SEED_PATTERN, FORCE_LAST_N_SEEDS);
  //for (auto offset: seedOffsets) std::cerr << " " << offset;
  //std::cerr << std::endl;
  ChainBuilder chainBuilder;
#if 1
  std::cerr << "\n      ++++++++++++++++++++++++++ Seed offsets... +++++++++++ " << std::endl;
  for (auto offset: seedOffsets)
  {
    unsigned fromHalfExtension = 0;
    const bool isExtended = false;
    const Seed seed (&read, offset, seedLength);
    const auto forwardData = seed.getPrimaryData(false);
    const auto reverseData = seed.getPrimaryData(true);
    const bool seedIsReverseComplement = (reverseData < forwardData);
    const auto hash = hashtable->getPrimaryHasher()->getHash64(seedIsReverseComplement ? reverseData : forwardData);
    const auto matchBits = hashtable->getMatchBits(hash, isExtended);
    const auto virtualByteAddress = hashtable->getVirtualByteAddress(hash);
    const auto bucketIndex = hashtable->getBucketIndex(virtualByteAddress);
    const auto hashThreadId = hashtable->getThreadIdFromVirtualByteAddress(virtualByteAddress);

    std::cerr << "  hits for offset " << std::setw(3) << offset << ": " << std::setfill('0') << (seedIsReverseComplement ? "reverse" : "forward") << " data: "
              << std::setw(14) << std::hex << (seedIsReverseComplement? reverseData : forwardData)
              << " hash: " << std::setw(14) << hash
              << " thread Id: " << std::setw(2) << (unsigned)hashThreadId
              << " bucketIndex: " << std::setw(8) << bucketIndex
              << " matchBits: " << std::setw(8) << matchBits
              << " orientation: " << (seedIsReverseComplement ? "REVERSE" : "FORWARD")
              << std::dec << std::setfill(' ') << std::endl;
    std::vector<HashRecord> hashRecords;
    std::vector<ExtendTableInterval> extendTableIntervals;
    hashtable->getHits(hash, isExtended, hashRecords, extendTableIntervals);
    std::cerr << "     intervals count: " << extendTableIntervals.size() << ":";
    for (const auto &interval: extendTableIntervals)
    {
      std::cerr << "  " << std::dec << std::setw(10) << interval.getStart() << ":" << std::setw(10) << interval.getLength();
    }
    std::cerr << std::endl;
    std::cerr << "     hits count for offset " << std::setw(3) << offset << ": " << std::dec << hashRecords.size();
    unsigned extendCount = 0;
    for (const auto &hit: hashRecords)
    {
      std::cerr << "  " << std::hex << std::setw(8) << (hit.getValue() >> 32) << ":" << std::setw(8) << (hit.getValue() & 0xFFFFFFFF);
    }
    std::cerr << std::endl;
    for (const auto &hit: hashRecords)
    {
      if(HashRecord::EXTEND == hit.getType())
      {
        ++extendCount;
        while (!hashRecords.empty() && (HashRecord::EXTEND == hashRecords.front().getType()))
        {
          const uint64_t primaryMask = (static_cast<uint64_t>(1) << hashtable->getPrimaryCrcBits()) - 1;
          const uint64_t secondaryMask = (static_cast<uint64_t>(1) << hashtable->getSecondaryCrcBits()) - 1;
          const uint64_t addressSegmentMask = primaryMask ^ secondaryMask;
          const uint64_t addressSegment = hash & addressSegmentMask;
          const auto extendRecord = hashRecords.front();
          hashRecords.clear();
          const unsigned halfExtension = extendRecord.getExtensionLength() / 2;
          const unsigned toHalfExtension = fromHalfExtension + halfExtension;
          std::cerr << "        " << std::hex << std::setw(8) << (hit.getValue() >> 32) << ":" << std::setw(8) << (hit.getValue() & 0xFFFFFFFF) << ": Trying extension from " << std::dec << fromHalfExtension << " to half-extension " << toHalfExtension << "..." << std::endl;
          if (seed.isValid(fromHalfExtension + halfExtension))
          {
            const unsigned halfPaddingBits = (6 - halfExtension) * 2;
            const uint64_t data = seed.getExtendedData(fromHalfExtension, toHalfExtension, seedIsReverseComplement);
            //const auto extendBases = seed.getExtendedData(fromHalfExtension, toHalfExtension, seedIsReverseComplement) << halfPaddingBits;
            const uint64_t extendBases = data << halfPaddingBits;
            const uint64_t shiftedExtensionId = extendRecord.getExtensionId() << 24;
            const uint64_t shiftedExtensionIdBin = ((hash & 0x7f) << (24 + 18));
            const auto extendedKey = shiftedExtensionIdBin | shiftedExtensionId | extendBases;
            const auto extendedHash = addressSegment | hashtable->getSecondaryHasher()->getHash64(extendedKey);
            hashtable->getHits(extendedHash, true, hashRecords, extendTableIntervals, true);
            std::cerr << "             extendedKey: " << std::hex << extendedKey << ": extendedHash: " << std::hex << extendedHash << std::dec << ": Found " << hashRecords.size() << " records:";
            for (const auto &hit: hashRecords)
            {
              std::cerr << "  " << std::hex << std::setw(8) << (hit.getValue() >> 32) << ":" << std::setw(8) << (hit.getValue() & 0xFFFFFFFF);
            }
            std::cerr << std::endl;
            std::cerr << "             Found " << extendTableIntervals.size() << " intervals:";
            for (const auto &interval: extendTableIntervals)
            {
              std::cerr << "  " << std::dec << std::setw(10) << interval.getStart() << ":" << std::setw(10) << interval.getLength();
            }
            std::cerr << std::endl;
            fromHalfExtension = toHalfExtension;
          }
        }
      }
    }
    if (extendCount)
    {
      std::cerr << "\n          >>>>>> EXTEND records found: " << extendCount;
    }
    std::cerr << std::dec << std::endl;
  }
  std::cerr << "\n      ++++++++++++++++++++++++++ Seed offsets...done +++++++++++ " << std::endl;
#endif
  mapper.get()->getPositionChains(read, chainBuilder);
  const auto compare = [] (const dragenos::map::SeedChain &lhs, const dragenos::map::SeedChain &rhs) -> bool {return rhs.size() < lhs.size();};
  chainBuilder.sort(compare);
  for (const auto &seedChain: chainBuilder)
  {
    const bool reverseComplement = seedChain.isReverseComplement();
    std::cerr << "SeedChain: RC: " << seedChain.isReverseComplement() << " RS: " << seedChain.hasOnlyRandomSamples() << " size: " << seedChain.size() << std::endl;
    for (const auto &seedPosition: seedChain)
    {
      const auto &seed = seedPosition.getSeed();
      std::cerr << "    " << std::setw(10) << std::dec << seedPosition.getReferencePosition()
                << ": seed (" << std::setw(3) << seed.getReadPosition() <<  ": " << seed.getPrimaryLength() << " + 2*" << std::setw(2) << seedPosition.getHalfExtension() << ")"
                << ": from " <<  std::setw(10) << seedPosition.getFirstProjection(reverseComplement) << " to " << std::setw(10) << seedPosition.getLastProjection(reverseComplement) << std::endl;
    }
  }
}

//::testing::Environment* globalTestEnvironment = nullptr;
Environment* globalTestEnvironment = nullptr;

std::pair<std::string, std::string> keyAndValue(const std::string &line)
{
  const auto equal = line.find(" = ");
  if (std::string::npos == equal)
  {
    BOOST_THROW_EXCEPTION(std::invalid_argument(std::string("unexpected line: ") + line));
  }
  auto key = line.substr(0, equal);
  while (key.back() == ' ') key.pop_back();
  const auto value = line.substr(equal + 3);
  return std::pair<std::string, std::string>(key, value);
}

std::vector<char> Environment::hashtableConfigText;
std::unique_ptr<Environment::HashtableConfig> Environment::hashtableConfig(nullptr);
std::unique_ptr<Environment::Hashtable> Environment::hashtable(nullptr);
unsigned char *Environment::referenceData = nullptr;
size_t Environment::referenceFileSize = 0;

std::vector<char> slurp(const boost::filesystem::path &filePath)
{
  const ssize_t fileSize = file_size(filePath);
  if (0 == fileSize)
  {
    return std::vector<char>();
  }
  std::vector<char> ret;
  ret.resize(fileSize);
  std::ifstream is(filePath.string());
  if (is && is.read(ret.data(), fileSize) && (fileSize == is.gcount()))
  {
    return ret;
  }
  BOOST_THROW_EXCEPTION(dragenos::common::IoException(errno, "Failed to read file"));
}

void Environment::SetUp()
{
  const auto argv = testing::internal::GetArgvs();
  namespace bfs = boost::filesystem;
  ASSERT_TRUE(argv.size() > 1 || (nullptr != getenv("REFDIR"))) << "Checking for reference-directory on the command line or in the environment variable REFDIR";
  const bfs::path referenceDir(argv.size() > 1 ? argv[1] : getenv("REFDIR"));
  std::cerr << "\n" << argv[0] << ": using rederence directory: " << referenceDir << "\n" << std::endl;
  //const bfs::path readsPath(argv[2]);
  //const bfs::path mappingsPath(argv[3]);
  ASSERT_TRUE(exists(referenceDir)) << "checking the existence of the reference-directory: " << referenceDir;
  const bfs::path hashtableConfigTextFile = referenceDir / "hash_table.cfg";
  ASSERT_TRUE(exists(hashtableConfigTextFile)) << "checking the existence of the hashtable config (text): " << hashtableConfigTextFile;
  hashtableConfigText = slurp(hashtableConfigTextFile);
  const bfs::path hashtableConfigFile = referenceDir / "hash_table.cfg.bin";
  ASSERT_TRUE(exists(hashtableConfigFile)) << "checking the existence of the hashtable config: " << hashtableConfigFile;
  const std::vector<char> config = slurp(hashtableConfigFile);
  hashtableConfig.reset(new HashtableConfig(config.data(), config.size()));
  const bfs::path hashtableFile = referenceDir / "hash_table.bin";
  const bfs::path extendTableFile = referenceDir / "extend_table.bin";
  ASSERT_TRUE(exists(hashtableFile)) << "checking the existence of the uncompressed hashtable: " << hashtableFile;
  hashtableFd = open(hashtableFile.c_str(), O_RDONLY, 0);
  ASSERT_LT(-1, hashtableFd) << "failed to open hashtable: " << hashtableFile << ": " << strerror(errno);
  const auto tableSize = boost::filesystem::file_size(hashtableFile);
  ASSERT_EQ(tableSize, hashtableConfig->getHashtableBytes());
  const int prot = PROT_READ;
  const int flags = MAP_PRIVATE | MAP_NORESERVE;
  const int offset = 0;
  table = static_cast<uint64_t*> (mmap(NULL, tableSize, prot, flags, hashtableFd, offset));
  ASSERT_NE(MAP_FAILED, table) << "failed to map hashtable file: " << strerror(errno);
  if (8 <= hashtableConfig.get()->getHashtableVersion())
  {
    extendTableFd = open(extendTableFile.c_str(), O_RDONLY, 0);
    ASSERT_LT(-1, extendTableFd) << "failed to open extendTable: " << extendTableFile << ": " << strerror(errno);
    const auto extendTableSize = boost::filesystem::file_size(extendTableFile);
    ASSERT_EQ(extendTableSize, hashtableConfig->getExtendTableBytes());
    extendTable = static_cast<uint64_t*> (mmap(NULL, extendTableSize, prot, flags, extendTableFd, offset));
    ASSERT_NE(MAP_FAILED, extendTable) << "failed to map extendTable file: " << strerror(errno);
  }
  // As Hashtable can't be initialized after construction, it has to be built dynamically
  hashtable.reset(new Hashtable(hashtableConfig.get(), table, extendTable));
  // We need the actual reference to check that the result of the queries is consistent with the actual reference
  const bfs::path referenceFile = referenceDir / "reference.bin";
  ASSERT_TRUE(exists(referenceFile)) << "checking the existence of the packed reference sequence: " << referenceFile;
  referenceFileSize = file_size(referenceFile);
  referenceFd = open(referenceFile.c_str(), O_RDONLY, 0);
  ASSERT_NE(-1, referenceFd) << "failed to open reference file " << referenceFile << ": " << strerror(errno);
  referenceData = static_cast<unsigned char *>(mmap(NULL, referenceFileSize, prot, flags, referenceFd, offset));
  ASSERT_NE(MAP_FAILED, referenceData) << "failed to mmap reference file " << referenceFile << ": " << strerror(errno);
}

void Environment::TearDown()
{
  if (-1 != hashtableFd)
  {
    close(hashtableFd);
    hashtableFd = -1;
  }
  if (-1 != extendTableFd)
  {
    close(extendTableFd);
    extendTableFd = -1;
  }
  if ((nullptr != table) && (MAP_FAILED != table))
  {
    const auto tableSize =  hashtableConfig->getHashtableBytes();
    munmap(table, tableSize);
    table = nullptr;
  }
  hashtable.reset(nullptr);
  if (-1 != referenceFd) close(referenceFd);
  if ((nullptr != referenceData) && (MAP_FAILED != referenceData))
  {
    munmap(referenceData, referenceFileSize);
    referenceData = nullptr;
  }
}

MapperFixture::HashtableConfig *MapperFixture::hashtableConfig;
MapperFixture::Hashtable *MapperFixture::hashtable;
MapperFixture::ReferenceSequence MapperFixture::referenceSequence;
std::unique_ptr<MapperFixture::Mapper> MapperFixture::mapper;
MapperFixture::Read MapperFixture::read;

MapperFixture::Read::Bases MapperFixture::getBasesFromReference(const size_t referencePosition, const unsigned readLength) const
{
  typedef MapperFixture::Read Read;
  typedef Read::Bases Bases;
  Bases bases;
  for (unsigned i = 0; readLength > i; ++i)
  {
    const unsigned char base = referenceSequence.getBase(i + READ_REFERENCE_POSITION);
    bases.push_back(referenceSequence.decodeBase(base));
    //bases.push_back(referenceSequence.translateTo2bpb(base));
  }
  return bases;
}

void MapperFixture::SetUpTestCase()
{
  hashtableConfig = Environment::getHashtableConfig();
  hashtable = Environment::getHashtable();
  referenceSequence.reset(hashtableConfig->getTrimmedRegions(), Environment::getReferenceData(), Environment::getReferenceSize());
  mapper.reset(new Mapper(hashtable));
  constexpr unsigned READ_LENGTH = 151;
  
  read.init(Read::Name(3, 'A'), Read::Bases(READ_LENGTH, 0), Read::Qualities(READ_LENGTH, 30), 0, 0);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  //::testing::FLAGS_gtest_throw_on_failure = true;
  ::testing::FLAGS_gtest_throw_on_failure = false;
  globalTestEnvironment = new Environment;
  /* ::testing::Environment* const env = */ ::testing::AddGlobalTestEnvironment(globalTestEnvironment);
  return RUN_ALL_TESTS();
  delete globalTestEnvironment;
}

