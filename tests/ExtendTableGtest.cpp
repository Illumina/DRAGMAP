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

#include "common/Bits.hpp"
#include "sequences/Read.hpp"
#include "sequences/Seed.hpp"
#include "reference/Hashtable.hpp"
#include "reference/ReferenceSequence.hpp"
#include "common/Exceptions.hpp"

class ExtendTableFixture: public ::testing::Test
{
public:
  ExtendTableFixture()
  {
  }
  static void SetUpTestCase();
  static void TearDownTestCase()
  {
  }
  void showInterval(const size_t start, const size_t length, const unsigned seedLength = 21, const unsigned halfWing = 0) const;
protected:
  typedef dragenos::reference::HashtableConfig HashtableConfig;
  typedef dragenos::reference::HashtableTraits HashtableTraits;
  typedef dragenos::reference::Hashtable Hashtable;
  typedef dragenos::reference::ReferenceSequence ReferenceSequence;
  static HashtableConfig *hashtableConfig;
  static Hashtable *hashtable;
  static ReferenceSequence referenceSequence;
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
  }
  static unsigned char *getReferenceData() {return referenceData;}
  static size_t getReferenceSize() {return referenceFileSize;}
  static Hashtable *getHashtable() {return hashtable.get();}
  static HashtableConfig *getHashtableConfig() {return hashtableConfig.get();}
  // Override this to define how to set up the environment.
  void SetUp() override;
  // Override this to define how to tear down the environment.
  void TearDown() override;
  static const std::vector<uint8_t> sequence;
  static const std::vector<uint8_t> qualities;
private:
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

Environment* globalTestEnvironment = nullptr;

std::string intervalToString(const dragenos::reference::HashRecord hashRecord)
{
  std::ostringstream os;
  using dragenos::reference::HashRecord;
  const auto type = hashRecord.getType();
  const auto value = hashRecord.getValue();
  switch (type)
  {
    case HashRecord::INTERVAL_SL:
      os << "SL";
      if (!hashRecord.isReverseComplement())
      {
        os << "0";
        os << " length: " <<  ((value >> 15) & 0x1ffUL);
        os << " start: " << (value & 0x7fffUL);
      }
      else
      {
        os << "1";
        os << " length: " <<  ((value >> 8) & 0xffffUL);
        os << " start: " << (value & 0xffUL);
      }
      break;
    case HashRecord::INTERVAL_SLE:
      os << "SLE";
      os << " exlift: " <<  ((value >> 16) & 0xffUL);
      os << " length: " <<  ((value >> 8) & 0xffUL);
      os << " start: " << (value & 0xffUL);
      break;
    case HashRecord::INTERVAL_S:
      os << "S  ";
      os << " start: " << (value & 0xffffffUL);
      break;
    case HashRecord::INTERVAL_L:
      os << "L  ";
      os << " length: " <<  (value & 0xffffffUL);
      break;
    default: os << "unknown"; break;
  }
  return os.str();
}

TEST_F(ExtendTableFixture, DISABLED_FindExtendTableIntervals)
{
  using dragenos::reference::HashRecord;
  using dragenos::reference::Bucket;

  const Bucket * const buckets = reinterpret_cast<Bucket *>(globalTestEnvironment->table);
  const size_t bucketCount = hashtableConfig->getHashtableBucketCount();
  size_t count = 0;
  for (size_t bucketId = 0; bucketCount > bucketId; ++ bucketId)
  {
    //bool inSet = false;
    //const auto begin = reinterpret_cast<const HashRecord *>(buckets + bucketId);
    //const auto end = begin + 8;
    bool hasS = false;
    bool hasL = false;
    for (const auto &hashRecord: buckets[bucketId])
    {
      const auto type = hashRecord.getType();
      if ((HashRecord::CHAIN_CON_MASK == type) || (HashRecord::CHAIN_CON_LIST == type))
      {
        break;
      }
      //if ((HashRecord::INTERVAL_SL == type) || (HashRecord::INTERVAL_SLE == type) ||(HashRecord::INTERVAL_S == type) ||(HashRecord::INTERVAL_L == type))
      //{
      //  if (!inSet) std::cerr << "  Bucket id: " << bucketId << std::endl;
      //  inSet = true;
      //  std::cerr << "    record: " << std::hex << std::setfill('0')
      //            << std::setw(8) << (hashRecord.getValue() >> 32) << ":" << (hashRecord.getValue() & 0xffffffff)
      //            << ": " << intervalToString(hashRecord)
      //            << std::setfill(' ') << std::dec << std::endl;
      //}
      hasS |= (HashRecord::INTERVAL_S == type);
      hasL |= (HashRecord::INTERVAL_L == type);
    }
    if (hasS && hasL)
    {
      std::cerr << "  Bucket id: " << bucketId << ":" << std::endl;
      for (const auto &hashRecord: buckets[bucketId])
      {
        const auto type = hashRecord.getType();
        const auto value = hashRecord.getValue();
        std::cerr << "    record: " << std::hex << std::setfill('0')
                  << std::setw(8) << (value >> 32) << ":" << (value & 0xffffffff);
        if ((HashRecord::INTERVAL_SL == type) || (HashRecord::INTERVAL_SLE == type) ||(HashRecord::INTERVAL_S == type) ||(HashRecord::INTERVAL_L == type))
        {
          std::cerr << ": " << intervalToString(hashRecord);
        }
        std::cerr << std::setfill(' ') << std::dec << std::endl;
      }
      ++count;
    }
    //count += inSet;
    //if (count > 100)
    //{
    //  break;
    //}
  }
}

unsigned char getBase(const unsigned char *referenceData, size_t position)
{
  const auto value = referenceData[position / 2];
  return (position % 2) ? (value >> 4) : (value & 0xF);
}

void ExtendTableFixture::showInterval(const size_t start, const size_t length, const unsigned seedLength, const unsigned halfWing) const
{
  const uint64_t * const extendTable = globalTestEnvironment->extendTable;
  const unsigned char *referenceData = globalTestEnvironment->getReferenceData();
  std::cerr << "start: " << start << " length: " << length << std::endl;
  for (size_t i = 0; 12 > i; ++i)
  { 
    if (2 == i) std::cerr << std::endl;
    const uint64_t extendRecord = extendTable[start + i - 2];
    const size_t position = (extendRecord & 0xffffffff);
    std::cerr << std::hex << std::setfill('0') << std::setw(8) << (extendRecord >> 32) << ":" << std::setw(8) << position << ":";
    for (unsigned j = 0; j < seedLength; ++j) std::cerr << " " << (unsigned)getBase(referenceData, position + j);
    std::cerr << std::setfill(' ') << std::dec << std::endl;
  }
  std::cerr << std::endl; 
  for (size_t i = 0; 12 > i; ++i)
  { 
    if (10 == i) std::cerr << std::endl;
    const uint64_t extendRecord = extendTable[start + length - 10 + i];
    const size_t position = (extendRecord & 0xffffffff);
    std::cerr << std::hex << std::setfill('0') << std::setw(8) << (extendRecord >> 32) << ":" << std::setw(8) << position << ":";
    for (unsigned j = 0; j < seedLength; ++j) std::cerr << " " << (unsigned)getBase(referenceData, position + j);
    std::cerr << std::setfill(' ') << std::dec << std::endl;
  }
}

TEST_F(ExtendTableFixture, ExploreExtendTable)
{
  // First example: EXTEND followed by and interval:
  // record: 7bf9df58:f2300726
  // record: 7bf9df59:f9000000: SLE exlift: 0 length: 0 start: 0
  // record: 7bf9df58:fa582838: S   start: 5777464
  // record: 7bf9df5a:fb0168c4: L   length: 92356
  //
  // Second example: 
  // record: 044d5968:f2280559
  // record: 044d5969:f9000000: SLE exlift: 0 length: 0 start: 0
  // record: 044d5968:fa4092ac: S   start: 4231852
  // record: 044d596a:fb013a82: L   length: 80514
  //
  // Third example:
  // record: b33ac730:f2301041
  // record: b33ac731:f9000001: SLE exlift: 0 length: 0 start: 1
  // record: b33ac731:fa7c4fc2: S   start: 0x27c4fc2
  // record: b33ac732:fb01a0a2: L   length: 106658
  //
  // Fourth example:
  // record: 9c50a1d8:f230116d
  // record: 9c50a1d9:f900001a: SLE exlift: 0 length: 0 start: 26
  // record: 9c50a1d9:fa46391e: S   start: 4602142
  // record: 9c50a1da:fb017940: L   length: 96576
  //
  // Fifth example:
  // record: c851b148:f2300595
  // record: c851b149:f900001a: SLE exlift: 0 length: 0 start: 26
  // record: c851b148:fabb10ae: S   start: 12259502
  // record: c851b14a:fb01f3b0: L   length: 127920
  //
  // First:
  //const size_t start = 5777464;
  //const size_t length = 92356;
  // Second:
  //const size_t start = 4231852;
  //const size_t length = 80514;
  // Third:
  //const size_t start = 0x27c4fc2; // both MSB on SLE and Carry on S
  //const size_t length = 106658;
  // Fourth:
  //const size_t start = 0x1b46391e; // both MSB on SLE and Carry on S
  //const size_t length = 96576;
  // Fifth:
  const size_t start = 0x1abb10ae; // MSB but no Carry
  const size_t length = 127920;
  showInterval(start, length);
}

TEST_F(ExtendTableFixture, CheckExtendInterval)
{
  const size_t start = 257097262; 
  const size_t length = 24;
  showInterval(start, length);
}

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
  ASSERT_TRUE(exists(referenceDir)) << "checking the existence of the reference-directory: " << referenceDir;
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
  if ((nullptr != extendTable) && (MAP_FAILED != extendTable))
  { 
    const auto extendTableSize = hashtableConfig->getExtendTableBytes();
    munmap(extendTable, extendTableSize);
    extendTable = nullptr;
  }
  hashtable.reset(nullptr);
  if (-1 != referenceFd) close(referenceFd);
  if ((nullptr != referenceData) && (MAP_FAILED != referenceData))
  {
    munmap(referenceData, referenceFileSize);
    referenceData = nullptr;
  }
}

ExtendTableFixture::HashtableConfig *ExtendTableFixture::hashtableConfig;
ExtendTableFixture::Hashtable *ExtendTableFixture::hashtable;
ExtendTableFixture::ReferenceSequence ExtendTableFixture::referenceSequence;
void ExtendTableFixture::SetUpTestCase()
{
  hashtableConfig = Environment::getHashtableConfig();
  hashtable = Environment::getHashtable();
  referenceSequence.reset(hashtableConfig->getTrimmedRegions(), Environment::getReferenceData(), Environment::getReferenceSize());
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

// Hashes a data buffer using a CRC polynomial of the same length, bit-by-bit.
// The same-length requirement makes the hash reversible.  No information is lost,
// and different data words are guaranteed to have different hashes.  This means
// that a key match can be confirmed by comparing the hash alone, or rather the
// portion not implicitly matched by use as address bits.
void* crcHashSlow(int bits, void const* poly, void const* data, void* hash)
{
  int bytes = (bits + 7) >> 3, topByte = bytes - 1;
  int topBitMask = (1 << ((bits + 7) % 8)), topByteMask = ((topBitMask << 1) - 1);
  int i, j, subtract;

#define POLY ((unsigned char*)poly)
#define HASH ((unsigned char*)hash)
  // Since data and polynomial are the same length, copy in all the data bytes immediately.
  // This data order doesn't match normal CRC computation, which by processing byte zero first,
  // effectively treats it as the most significant byte of the dividend.  But with odd bit
  // lengths, this works better.
  memcpy(hash, data, bytes);
  // Loop through the bits
  for (i = 0; i < bits; i++) {
    // Plan to subtract the polynomial if the MSB is 1
    subtract = (HASH[topByte] & topBitMask);

    // Left-shift the remainder (corresponds to right-shifting the polynomial position)
    for (j = topByte; j > 0; j--) HASH[j] = (HASH[j] << 1) | (HASH[j - 1] >> 7);
    HASH[0] <<= 1;

    // Subtract the polynomial if required to cancel the MSB shifted out
    if (subtract)
      for (j = 0; j < bytes; j++) HASH[j] ^= POLY[j];
  }
  // Mask off unused positions in the top byte
  HASH[topByte] &= topByteMask;
  // Return the hash pointer
  return hash;
#undef POLY
#undef HASH
}

// Optimized 64-bit version
void* crcHash64Init(int bits, void const* poly)
{
  int       bytes = (bits + 7) >> 3, i, j;
  int       bufQw = (1 + 256 * bytes);
  uint64_t *init = (uint64_t*)calloc(8, bufQw), *p = init;
  uint64_t  data;

  if (!init) return NULL;
  // Store the byte count in the init buffer
  *p++ = bytes;
  // Store the slow-hash of each byte value 0-255 in each byte position in a data buffer
  for (i = 0; i < bytes; i++) {
    for (j = 0; j < 256; j++) {
      data = (uint64_t)j << (i << 3);
      crcHashSlow(bits, poly, &data, p++);
    }
  }
  return init;
}

// Optimized 64-bit version
void* crcHash64(uint64_t* init, uint8_t const* data, uint64_t* hash)
{
  uint64_t h     = 0;
  int      bytes = *init++;

  while (bytes--) {
    h ^= init[*data++];
    init += 256;
  }
  *hash = h;
  return hash;
}

uint8_t encodebaseTo2Bits(const char c)
{
  switch (c)
  {
    case 'A': return 0; break;
    case 'C': return 1; break;
    case 'G': return 2; break;
    case 'T': return 3; break;
    default: throw std::invalid_argument(std::string("base must be ACGT: ") + c);
  };
}

//std::vector<uint8_t> generateRefSeq(std::vector<char> bases)
std::vector<uint8_t> generateRefSeq(std::string bases)
{
  assert(0 == (bases.size() % 4));
  std::vector<uint8_t> refSeq;
  refSeq.resize((bases.size() + 3) / 4);
  uint8_t seqByte = 0;
  for (unsigned i = 0; bases.size() > i; ++i)
  {
    const auto c = bases[i];
    const uint8_t encoded = encodebaseTo2Bits(c);
    seqByte |= ((encoded & 3) << ((i & 3) << 2));
    if (3 == (i & 3))
    {
      refSeq[i/4] = seqByte;
      seqByte = 0;
    }
  }
  return refSeq;
}
