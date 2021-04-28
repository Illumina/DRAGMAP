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

// copied from git/dragen/sec/common/hash_generation/crc_hash.c
// implementation at end of file
void* crcHashSlow(int bits, void const* poly, void const* data, void* hash);
void* crcHash64Init(int bits, void const* poly);
void* crcHash64(uint64_t* init, uint8_t const* data, uint64_t* hash);

// generate refSeq as it is done in hash_table.c
// encoding 4 bases per byte - ACGT -> 0123
//std::vector<uint8_t> generateRefSeq(std::vector<char> bases);
std::vector<uint8_t> generateRefSeq(std::string bases);

class HashtableFixture: public ::testing::Test
{
public:
  HashtableFixture()
  {
  }
  // The name was changed in gtest 1.8:   static void SetUpTestSuite()
  static void SetUpTestCase();
  static void TearDownTestCase()
  {
  }
protected:
  typedef dragenos::reference::HashtableConfig HashtableConfig;
  typedef dragenos::reference::HashtableTraits HashtableTraits;
  typedef dragenos::reference::Hashtable Hashtable;
  typedef dragenos::reference::ReferenceSequence ReferenceSequence;
  static HashtableConfig *hashtableConfig;
  static Hashtable *hashtable;
  static ReferenceSequence referenceSequence;
  //static const uint64_t fastqOffset = 64; // base quality offset, +64 or +33
  //static const bool     convertCtoT = false; // whether to convert C->T for methylation
  //static const bool     convertGtoA = false; // whether to convert G->A for methylation
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


struct TestRead
{
  std::string name_;
  std::string bases_;
  std::string qualities_;

  dragenos::sequences::Read::Name name() const {return dragenos::sequences::Read::Name(name_.begin(), name_.end());}
  dragenos::sequences::Read::Bases bases() const {return dragenos::sequences::Read::Bases(bases_.begin(), bases_.end());}
  dragenos::sequences::Read::Qualities qualities() const {return dragenos::sequences::Read::Qualities(qualities_.begin(), qualities_.end());}
};

TEST_F(HashtableFixture, getBits)
{
  using namespace dragenos::common::bits;
  const uint64_t v = 0xFEDCBA9876543210;
  ASSERT_EQ(v, (getBits<0, 64>(v)));
  ASSERT_EQ(0, (getBits<0, 4>(v)));
  ASSERT_EQ(1, (getBits<4, 4>(v)));
  ASSERT_EQ(2, (getBits<8, 4>(v)));
  ASSERT_EQ(3, (getBits<12, 4>(v)));
  ASSERT_EQ(4, (getBits<16, 4>(v)));
  ASSERT_EQ(5, (getBits<20, 4>(v)));
  ASSERT_EQ(6, (getBits<24, 4>(v)));
  ASSERT_EQ(7, (getBits<28, 4>(v)));
  ASSERT_EQ(8, (getBits<32, 4>(v)));
  ASSERT_EQ(9, (getBits<36, 4>(v)));
  ASSERT_EQ(10, (getBits<40, 4>(v)));
  ASSERT_EQ(11, (getBits<44, 4>(v)));
  ASSERT_EQ(12, (getBits<48, 4>(v)));
  ASSERT_EQ(13, (getBits<52, 4>(v)));
  ASSERT_EQ(14, (getBits<56, 4>(v)));
  ASSERT_EQ(15, (getBits<60, 4>(v)));
}

TEST_F(HashtableFixture, DISABLED_ExploreHhashtablev8)
{
  
}

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

TEST_F(HashtableFixture, DISABLED_FindExtendTableIntervals)
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

TEST_F(HashtableFixture, ExploreExtendTable)
{
  const uint64_t * const extendTable = globalTestEnvironment->extendTable;
  // First example: EXTEND followed by and interval:
  // record: 7bf9df58:f2300726
  // record: 7bf9df59:f9000000: SLE exlift: 0 length: 0 start: 0
  // record: 7bf9df58:fa582838: S   start: 5777464
  // record: 7bf9df5a:fb0168c4: L   length: 92356
  //
  // Secod example: 
  // record: 044d5968:f2280559
  // record: 044d5969:f9000000: SLE exlift: 0 length: 0 start: 0
  // record: 044d5968:fa4092ac: S   start: 4231852
  // record: 044d596a:fb013a82: L   length: 80514
  //
  //const size_t start = 4231852;
  const size_t start = 0;
  for (size_t i = 0; 10 > i; ++i)
  {
    std::cerr << std::hex << std::setfill('0') << std::setw(8) << (extendTable[start + i] >> 32) << ":" << std::setw(8) << (extendTable[start + i] & 0xffffffff) 
              << std::setfill(' ') << std::dec << std::endl;
  }
}

TEST_F(HashtableFixture, DISABLED_FindSeedsInTheBeginning)
//TEST_F(HashtableFixture, FindSeedsInTheBeginning)
{
  using dragenos::reference::HashRecord;
  ASSERT_EQ(sizeof(uint64_t), sizeof(HashRecord));
  const size_t hashtableBytes = hashtableConfig->getHashtableBytes();
  const size_t hashtableRecords = hashtableBytes / 8;
  const HashRecord * const table = reinterpret_cast<HashRecord *>(globalTestEnvironment->table);
  std::vector<HashRecord> records;
  records.reserve(100000);
  uint32_t minPosition = 999999999;
  bool printNext = false;
  auto printRecord = [](size_t i, const HashRecord &record)
  {
    std::cerr << "i: " << std::setw(10) << i << ", position: " << std::setw(10) << record.getPosition() << std::hex << ", bucket: " << std::setw(8) << (i/8) 
              << ", record: " <<  std::setw(8) << (record.getValue() >> 32) << " " << std::setfill('0') << std::setw(2) << (unsigned)record.getThreadId() 
              << " " << std::setfill('0') << std::setw(6) << record.getHashBits() << std::setfill(' ') << " " << record.isExtendedSeed() << " " << record.isLastInThread() << " " << record.isReverseComplement() <<  std::dec << std::endl;
  };
  for (size_t i = 0; hashtableRecords > i; ++i)
  {
    const auto &record = table[i];
    if (printNext)
    {
      printNext = record.isHit() && !record.isLastInThread();
      printRecord(i, record);
      continue;
    }
    if(record.isHit() && !record.isExtendedSeed() /* && !record.isReverseComplement() */ && (0 != record.getPosition()))
    {
      if (record.getPosition() < 180000)
      {
        records.push_back(record);
        //if (record.getPosition() <=163970) 
        if (record.getPosition() <=164015 + 20) 
        {
          printRecord(i, record);
          printNext = !record.isLastInThread();
        }
      }
      minPosition = std::min(minPosition, record.getPosition());
    }
  }
  std::cerr << "size: " << records.size() << ", min: " << minPosition << std::endl;
  // TODO: check that we have the good encoding and Crc hash for the seed at position 163840
/*

164012: 0035de20855c4d8e: 005937e5dcd: 0164df977: 05c4d8e: 
164013: 002643a64efec778: 003f600b72c: 00fd802dc: 07ec778: 
164014: 0024e8498f420b6f: 003d20b9d54: 00f482e75: 0420b6f: 
164015: 003f4027ea80c9a2: 0068c2421c6: 01a309087: 000c9a2: 
164016: 00299c8789f45cc2: 0044eb407c7: 0113ad01f: 0745cc2: 
164017: 002c2baf9129799a: 0049285ac86: 0124a16b2: 029799a: 
164018: 00003ec58be38a7d: 000067f72fa: 00019fdcb: 0638a7d: 
164019: 0006515647af6f7b: 000a76b6e6a: 0029dadb9: 02f6f7b: 
164020: 0001863cc7e83890: 00028654ab1: 000a1952a: 0683890: 
164021: 001b7173f2541ca2: 002d73e8095: 00b5cfa02: 0541ca2: 
164022: 002090d28fc16982: 0035efdcbe1: 00d7bf72f: 0416982: 
164023: 002e68bad0a434ca: 004cdd75698: 013375d5a: 02434ca: 
164024: 0000ae00db80d929: 000120316b8: 000480c5a: 000d929: 
164025: 0006756713b7bbae: 000ab272b89: 002ac9cae: 037bbae: 
164026: 001aa977ab443af0: 002c28ae33a: 00b0a2b8c: 0443af0: 
164027: 001dbaa1297f1c3a: 00313d1aeca: 00c4f46bb: 07f1c3a: 
164028: 00212226390ba9a4: 0036e08f4e7: 00db823d3: 00ba9a4: 
164029: 001848898e42ea69: 00283823d39: 00a0e08f4: 042ea69: 
164030: 00204cc546c7377e: 00357f26bd2: 00d5fc9af: 047377e: 
164031: 0023271f3e1819c4: 003a38cbbed: 00e8e32ef: 01819c4: 
164032: 0025b167d37bbcc0: 003e6dd3f63: 00f9b74fd: 07bbcc0: 
164033: 00127ca53770fdb6: 001e9e71a3d: 007a79c68: 070fdb6: 
164034: 0022d3a73e8851c7: 0039ae8cff9: 00e6ba33f: 00851c7: 
164035: 0035de20855c4d8e: 005937e5dcd: 0164df977: 05c4d8e: 
164036: 002643a64efec778: 003f600b72c: 00fd802dc: 07ec778: 
164037: 0024e8498f420b6f: 003d20b9d54: 00f482e75: 0420b6f: 
164038: 0009287b35d36195: 000f2b0c112: 003cac304: 0536195: 
164039: 00127ca53770fdb6: 001e9e71a3d: 007a79c68: 070fdb6: 
164040: 0022d3a73e8851c7: 0039ae8cff9: 00e6ba33f: 00851c7: 
164041: 0035de20855c4d8e: 005937e5dcd: 0164df977: 05c4d8e: 
164042: 002643a64efec778: 003f600b72c: 00fd802dc: 07ec778: 
164043: 0024e8498f420b6f: 003d20b9d54: 00f482e75: 0420b6f: 

 *********** Hashtable V8 *************************
i:   85800915, position:     164020, bucket:   a3a6fa, record: 5b41c480 16 683890 0 0 0
i:   85800916, position:  491275458, bucket:   a3a6fa, record: 5b41c480 16 683890 0 0 0
i:   85800917, position: 2361018399, bucket:   a3a6fa, record: 5b41c480 16 683890 0 0 0
i:   85800918, position: 1387614846, bucket:   a3a6fa, record: 5b41c481 16 683890 0 0 1
i:   85800919, position: 2727213102, bucket:   a3a6fa, record: 5b41c483 16 683890 0 1 1
i:  351104461, position:     164019, bucket:  29dadb9, record: 357b7bd8 0d 2f6f7b 0 0 0
i:  351104462, position: 2133813782, bucket:  29dadb9, record: 225d5733 08 4baae6 0 1 1
i:  358933876, position:     164025, bucket:  2ac9cae, record: c5bddd70 31 37bbae 0 0 0
i:  358933877, position:  955934793, bucket:  2ac9cae, record: dd60e332 37 2c1c66 0 1 0
i:  405708885, position:     164030, bucket:  305d40a, record: e167ee01 38 2cfdc0 0 0 1
i:  405708886, position:  491275468, bucket:  305d40a, record: e167ee01 38 2cfdc0 0 0 1
i:  405708887, position: 2361018409, bucket:  305d40a, record: e167ee03 38 2cfdc0 0 1 1
i: 1012682444, position:     164024, bucket:  78b89d9, record: 2c7107a9 0b 0e20f5 0 0 1
i: 1012682445, position:  491275462, bucket:  78b89d9, record: 2c7107a9 0b 0e20f5 0 0 1
i: 1012682446, position: 2271026583, bucket:  78b89d9, record: 2c7107a9 0b 0e20f5 0 0 1
i: 1012682447, position: 4110514254, bucket:  78b89d9, record: f5f5771a 3d 3eaee3 0 1 0
i: 1275560092, position:     164027, bucket:  980f013, record: 457b5321 11 2f6a64 0 0 1
i: 1275560093, position:  491275465, bucket:  980f013, record: 457b5321 11 2f6a64 0 0 1
i: 1275560094, position: 1068032020, bucket:  980f013, record: 457b5321 11 2f6a64 0 0 1
i: 1275560095, position: 2361018406, bucket:  980f013, record: 457b5323 11 2f6a64 0 1 1
i: 1276756702, position:     164018, bucket:  983385b, record: 5a55cd29 16 4ab9a5 0 0 1
i: 1276756703, position:  491275456, bucket:  983385b, record: 5a55cd2b 16 4ab9a5 0 1 1
i: 1349535650, position:     164029, bucket:  a0e08f4, record: 9e175348 27 42ea69 0 0 0
i: 1349535651, position:  491275467, bucket:  a0e08f4, record: 9e175348 27 42ea69 0 0 0
i: 1349535652, position: 2361018408, bucket:  a0e08f4, record: 9e175348 27 42ea69 0 0 0
i: 1349535653, position: 1717893350, bucket:  a0e08f4, record: 9e175349 27 42ea69 0 0 1
i: 1349535654, position: 1810367446, bucket:  a0e08f4, record: 9e175349 27 42ea69 0 0 1
i: 1349535655, position: 4110471557, bucket:  a0e08f4, record: a5a5a569 29 34b4ad 0 0 1
i: 1481727078, position:     164026, bucket:  b0a2b8c, record: 9e21d780 27 443af0 0 0 0
i: 1481727079, position: 4110671682, bucket:  b0a2b8c, record: f0f0e6c6 3c 1e1cd8 1 1 0
i: 1809824126, position:     164022, bucket:  d7bf72f, record: f20b4c10 3c 416982 0 0 0
i: 1809824127, position: 4110582041, bucket:  d7bf72f, record: 82828282 20 505050 0 1 0
i: 1841372825, position:     164028, bucket:  db823d3, record: 705d4d20 1c 0ba9a4 0 0 0
i: 1841372826, position:  491275466, bucket:  db823d3, record: 705d4d20 1c 0ba9a4 0 0 0
i: 1841372827, position: 2361018407, bucket:  db823d3, record: 705d4d20 1c 0ba9a4 0 0 0
i: 1841372828, position: 1173208795, bucket:  db823d3, record: 705d4d21 1c 0ba9a4 0 0 1
i: 1841372829, position: 1810367447, bucket:  db823d3, record: 705d4d21 1c 0ba9a4 0 0 1
i: 1841372830, position: 2524182161, bucket:  db823d3, record: 705d4d21 1c 0ba9a4 0 0 1
i: 1841372831, position: 4110596772, bucket:  db823d3, record: a4a4a4a4 29 149494 1 0 0
i: 1854408244, position:     164021, bucket:  dd100c6, record: da996679 36 532ccf 0 0 1
i: 1854408245, position:  491275459, bucket:  dd100c6, record: da996679 36 532ccf 0 0 1
i: 1854408246, position: 2271026580, bucket:  dd100c6, record: da996679 36 532ccf 0 0 1
i: 1854408247, position: 2361018400, bucket:  dd100c6, record: da99667b 36 532ccf 0 1 1
i: 2312536315, position:     164016, bucket: 113ad01f, record: e3a2e610 38 745cc2 0 0 0
i: 2312536316, position:  181960267, bucket: 113ad01f, record: f37301c2 3c 6e6038 0 1 0
i: 2454762901, position:     164017, bucket: 124a16b2, record: 414bccd0 10 29799a 0 0 0
i: 2454762902, position:  491275455, bucket: 124a16b2, record: 414bccd0 10 29799a 0 0 0
i: 2454762903, position: 4110529313, bucket: 124a16b2, record: 9a9a9a9a 26 535353 0 1 0
i: 2579163857, position:     164023, bucket: 13375d5a, record: 4d21a650 13 2434ca 0 0 0
i: 2579163858, position:  491275461, bucket: 13375d5a, record: 4d21a650 13 2434ca 0 0 0
i: 2579163859, position: 2271026582, bucket: 13375d5a, record: 4d21a650 13 2434ca 0 0 0
i: 2579163860, position: 2289834295, bucket: 13375d5a, record: 4d21a650 13 2434ca 0 0 0
i: 2579163861, position: 2361018402, bucket: 13375d5a, record: 4d21a650 13 2434ca 0 0 0
i: 2579163862, position: 2727213099, bucket: 13375d5a, record: 4d21a651 13 2434ca 0 0 1
i: 2579163863, position: 4110431460, bucket: 13375d5a, record: e3c5524c 38 78aa49 1 0 0
i: 3246648837, position:     164015, bucket: 18307dc0, record:  1f04ff9 00 3e09ff 0 0 1
i: 3246648838, position:   18408227, bucket: 18307dc0, record:  1f04ff9 00 3e09ff 0 0 1
i: 3246648839, position: 4110543918, bucket: 18307dc0, record: ffffffff 3f 7fffff 1 1 1
size: 15492, min: 164015
*/
/*
 ********* Hashtable V7 ************
i:   88104482, position:     163857, bucket:   a80bc4, record: 8738b4d0 21 67169a 0 0 0
i:   88104483, position: 1170195791, bucket:   a80bc4, record: 9950b5bd 26 2a16b7 1 0 1
i:   88104484, position: 2527274721, bucket:   a80bc4, record: 87940810 21 728102 0 0 0
i:   88104485, position: 1881874082, bucket:   a80bc4, record: 879aa2b3 21 735456 0 1 1
i:  326266458, position:     163868, bucket:  26e4dcb, record: 73786250 1c 6f0c4a 0 0 0
i:  326266459, position: 2651802530, bucket:  26e4dcb, record: 73a6f3a2 1c 74de74 0 1 0
i:  585378490, position:     163843, bucket:  45c8557, record: ec122c3a 3b 024587 0 1 0
i:  719018456, position:     163858, bucket:  55b6b3b, record: 69808dea 1a 3011bd 0 1 0
i:  945034968, position:     163854, bucket:  70a82db, record: 661a89d2 19 43513a 0 1 0
i: 1175941496, position:     163866, bucket:  8c2ee2f, record: e3fa4cda 38 7f499b 0 1 0
i: 1246235720, position:     163863, bucket:  9490189, record: 20c7db5a 08 18fb6b 0 1 0
i: 1251782843, position:     163855, bucket:  9539617, record: feb8969a 3f 5712d3 0 1 0
i: 1305065840, position:     163867, bucket:  9b9372e, record: c5e18942 31 3c3128 0 1 0
i: 1317662010, position:     163850, bucket:  9d13da7, record: fa917410 3e 522e82 0 0 0
i: 1317662011, position: 1866893293, bucket:  9d13da7, record: faa3077c 3e 5460ef 1 0 0
i: 1317662012, position: 2461380759, bucket:  9d13da7, record: faa3077d 3e 5460ef 1 0 1
i: 1317662013, position: 1936430087, bucket:  9d13da7, record: faa3077d 3e 5460ef 1 0 1
i: 1317662014, position:  655167895, bucket:  9d13da7, record: faa3077d 3e 5460ef 1 0 1
i: 1317662015, position: 4110659670, bucket:  9d13da7, record: efefefef 3b 7dfdfd 1 1 1
i: 1319795192, position:     163845, bucket:  9d54f3f, record: e93331e2 3a 26663c 0 1 0
i: 1528050906, position:     163848, bucket:  b62869b, record: 7b4ea3ea 1e 69d47d 0 1 0
i: 1620992618, position:     163860, bucket:  c13cc4d, record: ac3d96f2 2b 07b2de 0 1 0
i: 2229760632, position:     163851, bucket: 109cee4f, record: f8072852 3e 00e50a 0 1 0
i: 2259045819, position:     163852, bucket: 10d4c9b7, record: f34f6aca 3c 69ed59 0 1 0
i: 2461650632, position:     163840, bucket: 125739d9, record: 20f75d8a 08 1eebb1 0 1 0
i: 2615785281, position:     163856, bucket: 137d36e8, record: 195cea5a 06 2b9d4b 0 1 0
i: 2617211309, position:     163864, bucket: 137fef35, record: b0c33928 2c 186725 0 0 0
i: 2617211310, position: 2212288375, bucket: 137fef35, record: b6394d98 2d 4729b3 0 0 0
i: 2617211311, position: 4110489025, bucket: 137fef35, record: b8b8b8ad 2e 171715 1 0 1
i: 3540465321, position:     163846, bucket: 1a60e6d5, record: b5a119f0 2d 34233e 0 0 0
i: 3540465322, position: 2705840826, bucket: 1a60e6d5, record: a6a1c890 29 543912 0 0 0
i: 3540465323, position: 2421611295, bucket: 1a60e6d5, record: a6da2182 29 5b4430 0 1 0
i: 3546073317, position:     163847, bucket: 1a6b991c, record: 93cb332a 24 796665 0 1 0
i: 3977797539, position:     163844, bucket: 1da30bf4, record: 89f644f0 22 3ec89e 0 0 0
i: 3977797540, position:  821505477, bucket: 1da30bf4, record: 924ec0dc 24 49d81b 1 0 0
i: 3977797541, position: 1432481257, bucket: 1da30bf4, record: 95861ce0 25 30c39c 0 0 0
i: 3977797542, position:   81847960, bucket: 1da30bf4, record: 8a5f0a2a 22 4be145 0 1 0
*/
}

TEST_F(HashtableFixture, CheckHashtableConfig)
{
std::cerr << "\n---------------------------" << std::endl;
  const auto hashtableConfig = Environment::getHashtableConfig();
  using dragenos::sequences::CrcPolynomial;
  using dragenos::sequences::CrcHasher;
  CrcPolynomial priCrcPolynomial(hashtableConfig->getPrimaryCrcBits(), hashtableConfig->getPriCrcPoly());
  ASSERT_TRUE(priCrcPolynomial == "2C991CE6A8DD55"); // the value from hash_table.cfg
  CrcPolynomial secCrcPolynomial(hashtableConfig->getSecondaryCrcBits(), hashtableConfig->getSecCrcPoly());
  ASSERT_TRUE(secCrcPolynomial == "1524CA66E8D39"); // the value from hash_table.cfg
  CrcHasher priCrcHasher(priCrcPolynomial);
  CrcHasher secCrcHasher(secCrcPolynomial);
  ASSERT_EQ(21, hashtableConfig->getPrimarySeedBases());
}

std::ostream &printBases(std::ostream &os, const dragenos::sequences::Read &read)
{
  using dragenos::sequences::Read;
  for(auto b: read.getBases())  os << (unsigned)b << " ";
  for(auto b: read.getBases())  os << Read::decodeBase(b);
  return os;
}

// Encode a sequence with 2 bases per byte into a sequence with 4 bases per byte
std::vector<uint8_t> encode4bpbTo2bpb(const std::vector<unsigned char> sequence)
{
  auto encode = [] (unsigned char b)
  {
    switch (b)
    {
      case 2: return 1;
      case 4: return 2;
      case 8: return 3;
      default: return 0;
    }
  };
  std::vector<uint8_t> encoded;
  bool even = true;
  for (const uint8_t b4x2: sequence)
  {
    // convert the 2x 4 bits into 2x 2bits
    const uint8_t b2x2 = encode(b4x2 & 0xF) | (encode(b4x2 >> 4) << 2);
    if (even)
    {
      encoded.push_back(b2x2);
    }
    else
    {
      encoded.back() |= (b2x2 << 4);
    }
    even = ! even;
  }
  return encoded; 
}

TEST_F(HashtableFixture, CheckReferenceBasesHashtableV8)
{
  //Chr1 starts at position 163840 - 9984 first bases (Ns) are trimmed. 
  ASSERT_EQ(163840, hashtableConfig->getSequences()[0].seqStart);
  ASSERT_EQ(9984, hashtableConfig->getSequences()[0].begTrim);
  // Beginning of the untrimmed sequence from FASTA:
  //                                             NNNNNNNNNNNNNNNNTAACCCTAAC
  // CCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTA
  // ACCCTAACCCTAACCCTAACCCTAACCCAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCC
  // TAACCCTAACCCTAACCCTAACCCTAACCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACC
  // CTAACCCTAACCCCTAACCCTAACCCTAAACCCTAAACCCTAACCCTAACCCTAACCCTAACCCTAACCC
  // CAACCCCAACCCCAACCCCAACCCCAACCCCAACCCTAACCCCTAACCCTAACCCTAACCCTACCCTAAC
  // CCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCCCTAACCCTAACCCTAACCCTAACCCTAACCC
  // TAACCCTAACCCCTAACCCTAACCCTAACCCTAACCCTCGCGGTACCCTCAGCCGGCCCGCCCGCCCGGG
  //
  // corresponding hexdump from reference.bin:
  const std::vector<unsigned char> sequence {
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,  0x18, 0x21, 0x22, 0x18, 0x21, 0x22, 0x18, 0x21,
    0x22, 0x18, 0x21, 0x22, 0x18, 0x21, 0x22, 0x18,  0x21, 0x22, 0x18, 0x21, 0x22, 0x18, 0x21, 0x22,
    0x18, 0x21, 0x22, 0x18, 0x21, 0x22, 0x18, 0x21,  0x22, 0x18, 0x21, 0x22, 0x18, 0x21, 0x22, 0x18,
    0x21, 0x22, 0x18, 0x21, 0x22, 0x18, 0x21, 0x22,  0x18, 0x21, 0x22, 0x18, 0x21, 0x22, 0x11, 0x22,
    0x82, 0x11, 0x22, 0x82, 0x11, 0x22, 0x82, 0x11,  0x22, 0x82, 0x11, 0x22, 0x82, 0x11, 0x22, 0x82,
    0x11, 0x22, 0x22, 0x18, 0x21, 0x22, 0x18, 0x21,  0x22, 0x18, 0x21, 0x22, 0x18, 0x21, 0x22, 0x18,
    0x21, 0x82, 0x11, 0x22, 0x82, 0x11, 0x22, 0x82,  0x11, 0x22, 0x82, 0x11, 0x22, 0x82, 0x11, 0x22,
    0x82, 0x11, 0x22, 0x82, 0x11, 0x22, 0x82, 0x11,  0x22, 0x82, 0x11, 0x22, 0x22, 0x18, 0x21, 0x22};

  const std::vector<uint8_t> encoded2bpb = encode4bpbTo2bpb(sequence);
  ASSERT_EQ(sequence.size(), encoded2bpb.size() * 2);
  ASSERT_EQ(0x00, encoded2bpb[0]);
  ASSERT_EQ(0x00, encoded2bpb[1]);
  ASSERT_EQ(0x00, encoded2bpb[2]);
  ASSERT_EQ(0x00, encoded2bpb[3]);
  ASSERT_EQ(0x43, encoded2bpb[4]); // 18, 21: TAAC: 01000011
  ASSERT_EQ(0x35, encoded2bpb[5]); // 22, 18: CCTA: 00110101
  ASSERT_EQ(0x54, encoded2bpb[6]); // 21, 22: ACCC: 01010100

  ASSERT_EQ(128, sequence.size());
  ASSERT_EQ(0x00014000 * 2, 163840);
  for (size_t i = 0; 0x14000 > i; ++i)
  {
    ASSERT_EQ(0, referenceSequence.getData()[i]) << "i: " << i;
    ASSERT_EQ(0, referenceSequence.getBase(2*i)) << "i: " << i;
    ASSERT_EQ(0, referenceSequence.getBase(2*i + 1)) << "i: " << i;
  }
  // check that all the bases trimmed at the beginning are N
  for (size_t i = 0; 9984 > i; ++i)
  {
    ASSERT_EQ(15, referenceSequence.getBase(163840 + i));
  }
  for (size_t i = 0x14000; 0x14008 > i; ++i)
  {
    ASSERT_EQ(0xFF, globalTestEnvironment->getReferenceData()[i]) << "i: " << i << ": " << (i - 0x14000);
    ASSERT_EQ(0xFF, referenceSequence.getData()[i]) << "i: " << i << ": " << (i - 0x14000);
    ASSERT_EQ(0xF, referenceSequence.getBase(9984 + 2*i)) << "i: " << i;
    ASSERT_EQ(0xF, referenceSequence.getBase(9984 + 2*i + 1)) << "i: " << i;
  }
  ASSERT_EQ(0x18, referenceSequence.getData()[0x14000 + 8]);
  ASSERT_EQ(0x8, referenceSequence.getBase(163840 + 9984 + 8 * 2));
  ASSERT_EQ(0x1, referenceSequence.getBase(163840 + 9984 + 8 * 2 + 1));
  ASSERT_EQ(0x21, referenceSequence.getData()[0x14000 + 9]);
  ASSERT_EQ(0x1, referenceSequence.getBase(163840 + 9984 + 9 * 2));
  ASSERT_EQ(0x2, referenceSequence.getBase(163840 + 9984 + 9 * 2 + 1));

  // The first primary seed that has a forward hit is at position 164015 = 163840 + 175
  // Generate the seeds from a few bases before to a few bases after and hash them
  // start at offset ((175/4) * 4)==172 and do 32 primary seeds
  const unsigned seedLength = hashtableConfig->getPrimarySeedBases();
  ASSERT_EQ(21, seedLength);
  const uint64_t seedMask = (((uint64_t)1) << (seedLength * 2)) - 1;
  // check that the mask has 2 * seedLength LSB set
  for (unsigned i = 0; 2 * seedLength > i; ++i)
  {
    ASSERT_EQ(1, (seedMask >> i) & 1) << "i: " << i;
  }
  // check that the mask bits after that al 0
  ASSERT_EQ(0, seedMask >> (2 * seedLength));
  const unsigned newBaseShift = 2 * seedLength - 2;
  using dragenos::sequences::CrcPolynomial;
  CrcPolynomial priCrcPolynomial(hashtableConfig->getPrimaryCrcBits(), hashtableConfig->getPriCrcPoly());
  ASSERT_TRUE(priCrcPolynomial == "2C991CE6A8DD55"); // the value from hash_table.cfg
  ASSERT_EQ(54, hashtableConfig->getPrimaryCrcBits());
  // void* crcHash64Init(int bits, void const* poly);
  const uint64_t poly = 0x2C991CE6A8DD55;
  uint64_t *init = reinterpret_cast<uint64_t*>(crcHash64Init(54, &poly));
  ASSERT_EQ(encoded2bpb.size(), 64);
  const size_t begin = 175 / 4;
  ASSERT_GT(64, begin + 8 + 5); // enough space to get 32 deeds of 21 bases
  uint64_t seedValue = (*reinterpret_cast<const uint64_t *>(encoded2bpb.data() + begin)) & seedMask;
  const uint8_t *bytePtr = encoded2bpb.data() + begin + 5; // pointer to the byte that contains the next bases
  uint8_t byte = (*bytePtr) >> 2;
  //std::cerr << "\n     ------------ printing hash values for region " << (163840 + 4*begin) << " to " << (163840 + 4*begin + 32) << std::endl;
  const unsigned table_size_64ths = 53;
  ASSERT_EQ(table_size_64ths, hashtableConfig->getTableSize64Ths()); 
  for (size_t i = 0; 32 > i; ++i)
  {
    // verify that the seed value matches the reference sequence
    auto tmp = seedValue;
    for (size_t j = 0; seedLength > j; ++j)
    {
// TODO FIX the introduction of the new bases
      ASSERT_EQ(1U << (tmp & 3), referenceSequence.getBase(163840 + 9984 + 4*begin + i + j)) << "i: " << i << ", j: " << j;
      tmp = (tmp >> 2);
    }
    uint64_t hashValue;
    crcHash64(init, reinterpret_cast<uint8_t const*>(&seedValue), &hashValue);
    // get the same value from the hashtable
    const auto primaryHasher = hashtable->getPrimaryHasher();
    const uint64_t hashValueFromHashtable = primaryHasher->getHash64(seedValue);
    ASSERT_EQ(hashValue, hashValueFromHashtable);
    //std::cerr << "Seed Value : Hash Value : " << std::hex << seedValue << " " << hashValueFromHashtable << std::dec << std::endl;
    // print the hashValue, virtual address, bucket and hash bits
    //const uint64_t addressBits = dragenos::common::bits::getBits<19, 35>(hashValue);
    //const uint64_t virtualAddress = (addressBits * table_size_64ths) / 64;
    //const uint64_t bucket = virtualAddress >> 6;
    //const uint64_t hashBits = dragenos::common::bits::getBits<0, 23>(hashValue);
    //std::cerr << (163840 + 4*begin + i) << ": " 
    //          << std::hex << std::setfill('0') << std::setw(16) << hashValue << std::setfill(' ') << ": "
    //          << std::hex << std::setfill('0') << std::setw(11) << virtualAddress << std::setfill(' ') << ": "
    //          << std::hex << std::setfill('0') << std::setw(9) << bucket<< std::setfill(' ') << ": "
    //          << std::hex << std::setfill('0') << std::setw(7) << hashBits<< std::setfill(' ') << ": "
    //          << std::dec << std::endl;
    // move to the next seed
    seedValue = seedValue >> 2;
    uint64_t newBase = (byte & 3);
    if (2 == (i % 4))
    {
      ++bytePtr;
      byte = *bytePtr;
    }
    else
    {
      byte = (byte >> 2);
    }
    seedValue |= (newBase << newBaseShift);
  }
}

//TEST_F(HashtableFixture, CheckReferenceBasesHashtableV7)
TEST_F(HashtableFixture, DISABLED_CheckReferenceBasesHashtableV7)
{
  using dragenos::sequences::Read;
  //ChrM starts at position 163840
  // i: 2461650632, position:     163840, bucket: 125739d9, record: 20f75d8a 08 1eebb1 0 1 0
  std::cerr << "\n  Expected seed at start of ChrM: 2461650632, position:     163840, bucket: 125739d9, record: 20f75d8a 08 1eebb1 0 1 0" << std::endl;
  size_t position = 163840;
  // 00014000  14 28 21 41 84 82 81 12  22 82 81 18 21 12 82 12 
  // 00014010  42 44 41 82 82 22 81 24  81 88 44 18 88 88 42 28

  const std::vector<unsigned char> sequence {1,4, 2,8, 2,1, 4,1, 8,4, 8,2, 8,1, 1,2,   2,2, 8,2, 8,1, 1,8, 2,1, 1,2, 8,2, 1,2,
                                             4,2, 4,4, 4,1, 8,2, 8,2, 2,2, 8,1, 2,4,   8,1, 8,8, 4,4, 1,8, 8,8, 8,8, 4,2, 2,8 };
  ASSERT_EQ(64, sequence.size());
  for(auto s: sequence)
  {
    ASSERT_EQ(s, referenceSequence.getBase(position));
    ++position;
  }
  auto translate = [] (unsigned char c) -> char
  {
    switch (c) {
      case 1: return 'A';
      case 2: return 'C';
      case 4: return 'G';
      case 8: return 'T';
      default: return 'N';
    }
  };
  TestRead r;
  for(auto s: sequence)
  {
    r.bases_.push_back(translate(s));
  }
  r.qualities_.resize(r.bases_.size(), 'B');
  ASSERT_EQ(r.bases_.size(), sequence.size());
  ASSERT_EQ(r.qualities_.size(), sequence.size());
  std::vector<uint8_t> refSeq = generateRefSeq(r.bases_);
  //uint64_t *init = (uint64_t*)crcHash64Init(hashtableConfig->getPrimaryCrcBits(), hashtableConfig->getPriCrcPoly());
  uint64_t *init = (uint64_t*)crcHash64Init(42, hashtableConfig->getPriCrcPoly());
  auto data = refSeq;
  //constexpr unsigned seedLen = 21;
  //const uint64_t priSeedMask = ((uint64_t)1 << (seedLen * 2)) - 1;

  // TODO: sort out the reference bin
  //
  //uint64_t refBinMask  = binBits ? ((0ULL - 1) << binBits) : 0;
  //hashKey |= (pos & refBinMask) << KEY_ANCHOR_OFFSET;

  uint64_t hash1;
  crcHash64(init, data.data(), &hash1);
  auto virtualByteAddress1 = dragenos::common::bits::getBits<19, 35>(hash1);
  std::cerr << "hash with all bases: " << std::hex << hash1 << " : " << dragenos::common::bits::getBits<0, 23>(hash1) << std::dec << std::endl; 
  std::cerr << "Virtual address from hashtable: " << std::hex << virtualByteAddress1 << std::dec << std::endl;
  std::cerr << "Bucket index from hashtable: " << std::hex << hashtable->getBucketIndex(virtualByteAddress1) << std::dec << std::endl;
  std::cerr << "thread Id from hashtable: " << std::hex << hashtable->getThreadIdFromVirtualByteAddress(virtualByteAddress1) << std::dec << std::endl;
  data[5] &= 3;
  data[6] = 0;
  data[7] = 0;
  crcHash64(init, data.data(), &hash1);
  virtualByteAddress1 = dragenos::common::bits::getBits<19, 35>(hash1);
  std::cerr << "hash with 21 bases: " << std::hex << hash1 << " : " << dragenos::common::bits::getBits<0, 23>(hash1) << std::dec << std::endl; 
  std::cerr << "Virtual address from hashtable: " << std::hex << virtualByteAddress1 << std::dec << std::endl;
  std::cerr << "Bucket index from hashtable: " << std::hex << hashtable->getBucketIndex(virtualByteAddress1) << std::dec << std::endl;
  std::cerr << "thread Id from hashtable: " << std::hex << hashtable->getThreadIdFromVirtualByteAddress(virtualByteAddress1) << std::dec << std::endl;


  free(init);
  init = NULL;

  r.name_ = "dummy";
  Read read;
  read.init(r.name(), r.bases(), r.qualities(), 0, 0);
  std::cerr << "\n    Read: ";  printBases(std::cerr, read); std::cerr << std::endl;
  const dragenos::sequences::Seed seed(&read, 0U, hashtableConfig->getPrimarySeedBases());

  const auto hashtableConfig = Environment::getHashtableConfig();
  using dragenos::sequences::CrcPolynomial;
  using dragenos::sequences::CrcHasher;
  CrcPolynomial priCrcPolynomial(hashtableConfig->getPrimaryCrcBits(), hashtableConfig->getPriCrcPoly());
  ASSERT_TRUE(priCrcPolynomial == "2C991CE6A8DD55"); // the value from hash_table.cfg
  //CrcPolynomial secCrcPolynomial(hashtableConfig->getSecCrcPoly(), hashtableConfig->getSecondaryCrcBits());
  //ASSERT_TRUE(secCrcPolynomial == "1524CA66E8D39"); // the value from hash_table.cfg
  CrcHasher priCrcHasher(priCrcPolynomial);
  //CrcHasher secCrcHasher(hashtableConfig->getSecondaryCrcBits(), secCrcPolynomial);
  const auto seedData = seed.getPrimaryData(false);
  const auto hash = priCrcHasher.getHash64(seedData);
  std::cerr << "Hash value from crc hasher: " << std::hex << hash << std::dec << std::endl;
  std::cerr << "Hash value from hashtable: " << std::hex << hashtable->getPrimaryHasher()->getHash64(seedData) << std::dec << std::endl;
  std::cerr << "expected Virtual address: " << std::hex << (hash * 60 / 64) << std::dec << std::endl;
  const auto virtualByteAddress = hashtable->getVirtualByteAddress(hash);
  std::cerr << "Virtual address from hashtable: " << std::hex << virtualByteAddress << std::dec << std::endl;
  std::cerr << "Bucket index from hashtable: " << std::hex << hashtable->getBucketIndex(virtualByteAddress) << std::dec << std::endl;
  std::cerr << "thread Id from hashtable: " << std::hex << hashtable->getThreadIdFromVirtualByteAddress(virtualByteAddress) << std::dec << std::endl;
#if 0
  typedef dragenos::reference::HashRecord HashRecord;
// TODO: add some tests for getHits
  std::vector<HashRecord> hashRecords;
  hashtable->getHits(seed, hashRecords);
  for (const auto &hashRecord: hashRecords)
  {
    std::cerr << hashRecord.getType() << std::endl;
  }
  uint64_t slowData = 0x375437b21d8;
  uint64_t slowHash = 0;
  uint64_t slowRev = 0;
  crcHashSlow(hashtableConfig->getPrimaryCrcBits(), hashtableConfig->getPriCrcPoly(), &slowData, &slowHash);
  crcHashSlow(hashtableConfig->getPrimaryCrcBits(), hashtableConfig->getPriCrcPoly(), &slowHash, &slowRev);
  std::cerr << "      crcHashSlow: slowData: " << std::hex << slowData << std::dec << std::endl;
  std::cerr << "      crcHashSlow: slowHash: " << std::hex << slowHash << std::dec << std::endl;
  std::cerr << "      crcHashSlow: slowRev: " << std::hex << slowRev << std::dec << std::endl;
  const uint64_t MASK_RECORD = ~((~0UL) << 23);
  uint64_t recordHash = (slowHash & MASK_RECORD);
  uint64_t addressHash = slowHash >> 19;
  std::cerr << "      crcHashSlow: recordHash: " << std::hex << recordHash << std::dec << std::endl;
  std::cerr << "      crcHashSlow: addressHash: " << std::hex << addressHash << std::dec << std::endl;
#endif
}

TEST_F(HashtableFixture, ExploreExtend)
{
  using dragenos::reference::Bucket;
  using dragenos::reference::HashRecord;
  const auto table = globalTestEnvironment->table;
  const auto bucketCount = hashtableConfig->getHashtableBytes() / sizeof(Bucket);
  auto buckets = reinterpret_cast<Bucket*>(table);
  bool found = false;
  size_t count = 0;
  uint64_t maxExtensionLength = 0;
  for (size_t i = 0; i < bucketCount && !found; ++i)
  {
    const auto &bucket = buckets[i];
    for (const auto &hashRecord: bucket)
    {
      //if (HashRecord::EXTEND == bucket[0].getType())
      if (HashRecord::EXTEND == hashRecord.getType())
      {
        ++count;
        maxExtensionLength = std::max(maxExtensionLength, hashRecord.getExtensionLength());
        if (1 == (hashRecord.getExtensionLength() % 2)) // as the extension length is supposed to be the total it should be even
        {
          std::cerr << "  QQQQQQQQQQQQ: " << (unsigned)hashRecord.getExtensionLength() <<  std::endl;
        }
#if 0
        std::cerr << "\n  ------------\ni: " << i << std::endl;
        for (const auto &hashRecord: bucket)
        {
          std::cerr << "  type: " << hashRecord.getType();
          //if (hashRecord.isExtendedSeed())
          {
            std::cerr << " value: " << std::hex << hashRecord.getValue() << " " << std::oct << hashRecord.getValue() << std::hex << " ";
            std::cerr << "thread Id: " << (unsigned)hashRecord.getThreadId() << " ";
            std::cerr << "hash bits: " << hashRecord.getHashBits() << " ";
            std::cerr << hashRecord.isExtendedSeed() << " ";
            std::cerr << hashRecord.isLastInThread() << " ";
            std::cerr << hashRecord.hasRandomSamples() << " ";
            std::cerr << 0xf << " " << HashRecord::EXTEND << " ";
            std::cerr << hashRecord.hasRepairRecords() << " ";
            std::cerr << hashRecord.isAltLiftover() << " ";
            std::cerr << (uint32_t)hashRecord.getExtensionLength() << " ";
            std::cerr << (uint32_t)hashRecord.getExtensionId() << " ";
            std::cerr << std::dec << std::endl;
            found |= (count > 10);
          }
          std::cerr << std::endl;
        }
#endif

      }
    }
  }
  std::cerr << "\n ------------------- \n    Found " << count << " EXTEND records in " << bucketCount << " buckets. Max Extension Length: " << (unsigned)maxExtensionLength << std::endl;
  //for (const auto p: {787527575UL, 2796142304UL, 213253491UL, 1333620284UL})
  //{
    
  //}
  //const auto table = Environment::table;
}

TEST_F(HashtableFixture, DISABLED_ExploreHits)
{
  using dragenos::reference::Bucket;
  const auto table = globalTestEnvironment->table;
  auto buckets = reinterpret_cast<Bucket*>(table);
  bool found = false;
  for (size_t i = 0; i < 10000000 && !found; ++i)
  {
    const auto &bucket = buckets[i];
    //for (const auto &hashRecord: bucket)
    {
      if (bucket[0].isHit())
      {
        std::cerr << "\n  ------------\ni: " << i << std::endl;
        for (const auto &hashRecord: bucket)
        {
          std::cerr << "  type: " << hashRecord.getType();
          if (hashRecord.isHit())
          {
            std::cerr << " value: " << std::hex << hashRecord.getValue() << " " << std::oct << hashRecord.getValue() << std::hex << " ";
            std::cerr << "thread Id: " << (unsigned)hashRecord.getThreadId() << " ";
            std::cerr << "hash bits: " << hashRecord.getHashBits() << " ";
            std::cerr << hashRecord.isExtendedSeed() << " ";
            std::cerr << hashRecord.isLastInThread() << " ";
            std::cerr << hashRecord.isReverseComplement() << " ";
            const size_t position = hashRecord.getPosition();
            std::cerr << "position: " << std::dec << position << std::endl;
            for (unsigned i = position -5 ; i < position + 25; ++i)
            {
              std::cerr << (unsigned int)referenceSequence.getBase(i) << " ";
            }
            found |= (!hashRecord.isReverseComplement());
          }
          std::cerr << std::endl;
        }
      }
    }
  }
  //for (const auto p: {787527575UL, 2796142304UL, 213253491UL, 1333620284UL})
  //{
    
  //}
  //const auto table = Environment::table;
}

TEST_F(HashtableFixture, DISABLED_Config)
{
  ASSERT_NE(nullptr, hashtableConfig);
  const auto &buffer = Environment::getHashtableConfigText();
  std::istringstream is(std::string(buffer.data(), buffer.size()));
  std::string line;
  getline(is, line); // 1st line is dragen version
  getline(is, line); // 2nd line is command line
  getline(is, line);
  ASSERT_EQ(std::string("#    Hash table version 7"), line);
  ASSERT_LE(7, hashtableConfig->getHashtableVersion());
  getline(is, line); // 4th line is "#"
  getline(is, line); // 5th line is "# Do not modify.
  getline(is, line); // 6th line is empty
  std::vector<std::pair<std::string, std::string>> keyValues;
  while (is && getline(is, line))
  {
    keyValues.push_back(keyAndValue(line));
  }
  const std::map<std::string, std::string> keyValuesMap(keyValues.begin(), keyValues.end());
  ASSERT_EQ(keyValuesMap.size(), keyValues.size()); // check that all keys are unique
  ASSERT_EQ(stol(keyValuesMap.at("reference_len_raw")), hashtableConfig->getReferenceLength());
  ASSERT_EQ(stoi(keyValuesMap.at("reference_sequences")), hashtableConfig->getNumberOfSequences());
  ASSERT_LT(5 * hashtableConfig->getNumberOfSequences(), keyValues.size()); // at least 5 entries per sequence
  for(unsigned i = 0; hashtableConfig->getNumberOfSequences() > i; ++i)
  {
    const std::string indexString = std::to_string(i);
    const auto numberOfSequences = hashtableConfig->getNumberOfSequences();
    const auto numberOfKeys = keyValues.size();
    const size_t index = numberOfKeys - 5 * (numberOfSequences - i);
    // check that we have the expected keys in the text config file
    ASSERT_EQ(std::string("reference_sequence") + indexString, keyValues[index].first);
    ASSERT_EQ(std::string("reference_start") + indexString, keyValues[index + 1].first);
    ASSERT_EQ(std::string("reference_beg_trim") + indexString, keyValues[index+ 2 ].first);
    ASSERT_EQ(std::string("reference_end_trim") + indexString, keyValues[index + 3].first);
    ASSERT_EQ(std::string("reference_len") + indexString, keyValues[index + 4].first);
    // check that the hashtable config values match the values in the text file for this sequence
    const std::string expectedName = keyValues[index].second.substr(1, keyValues[index].second.size() - 2); // enclosed in ''
    ASSERT_EQ(expectedName, hashtableConfig->getSequenceName(i));
    const auto &sequence = hashtableConfig->getSequences()[i];
    ASSERT_EQ(keyValues[index + 1].second, std::to_string(sequence.seqStart));
    ASSERT_EQ(keyValues[index + 2].second, std::to_string(sequence.begTrim));
    ASSERT_EQ(keyValues[index + 3].second, std::to_string(sequence.endTrim));
    ASSERT_EQ(keyValues[index + 4].second, std::to_string(sequence.seqLen));
  }
}

TEST_F(HashtableFixture, HashtableContent)
{
  ASSERT_NE(nullptr, hashtable);
}

TEST_F(HashtableFixture, DISABLED_HashSeed)
{
  using namespace dragenos::sequences;
  const CrcHasher primaryHasher(hashtableConfig->getPrimaryPolynomial());
  const CrcHasher secondaryHasher(hashtableConfig->getSecondaryPolynomial());
  TestRead tr = {"@dummy", "ACGTAACCGGTTAAACCCGGGTTTAAAACCCCGGGGTTTT", "BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB"};
  dragenos::sequences::Read read;
  read.init(tr.name(), tr.bases(), tr.qualities(), 0, 0);
  ASSERT_LT(hashtableConfig->getPrimarySeedBases(), read.getLength());
  const dragenos::sequences::Seed primarySeed(&read, 0U, hashtableConfig->getPrimarySeedBases());
  const dragenos::sequences::Seed extendedSeed(&read, 7U, hashtableConfig->getPrimarySeedBases());
  ASSERT_EQ(primaryHasher.getHash64(primarySeed.getPrimaryData(false)), hashtable->getPrimaryHasher()->getHash64(primarySeed.getPrimaryData(false)));
  ASSERT_EQ(secondaryHasher.getHash64(extendedSeed.getPrimaryData(false)), hashtable->getSecondaryHasher()->getHash64(extendedSeed.getPrimaryData(false)));
}

TEST_F(HashtableFixture, DISABLED_GetExtendRecord)
{
#if 0
  using namespace dragenos::sequences;
  const CrcHasher primaryHasher(hashtableConfig->getPrimaryPolynomial());
  const CrcHasher secondaryHasher(hashtableConfig->getSecondaryPolynomial());
  // sequence taken from ChrM
  TestRead tr = {"@dummy", "TTGGTATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTG", "BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB"};
  dragenos::sequences::Read read;
  read.init(tr.name(), tr.bases(), tr.qualities(), 0, 0);
  ASSERT_LT(hashtableConfig->getPrimarySeedBases(), read.getLength());
  const dragenos::sequences::Seed seed(&read, 0U, hashtableConfig->getPrimarySeedBases());
  typedef dragenos::reference::HashRecord HashRecord;
  std::vector<HashRecord> hashRecords;
  hashtable->getHits(seed, hashRecords);
  ASSERT_EQ(2, hashRecords.size());
  // first record is an EXTEND
  ASSERT_FALSE(hashRecords[0].isHit());
  ASSERT_EQ(HashRecord::EXTEND, hashRecords[0].getOpCode());
  ASSERT_EQ(HashRecord::EXTEND, hashRecords[0].getType());
  ASSERT_FALSE(hashRecords[0].isLastInThread());
  ASSERT_TRUE(hashRecords[0].hasRandomSamples());
  ASSERT_FALSE(hashRecords[0].isLastInThread());
  ASSERT_FALSE(hashRecords[0].isExtendedSeed());
  ASSERT_EQ(2621441, hashRecords[0].getFrequency());
  // NOT checking for AL (Alt Liftover) or RF (Repair Flag) as they are deprecated
  // second sample is a HIT (random sample) and the last entry in thread
  ASSERT_TRUE(hashRecords[1].isHit());
  ASSERT_EQ(HashRecord::HIT, hashRecords[1].getType());
  ASSERT_FALSE(hashRecords[1].isReverseComplement());
  ASSERT_TRUE(hashRecords[1].isLastInThread());
  ASSERT_FALSE(hashRecords[1].isExtendedSeed());
  ASSERT_EQ(1686365997, hashRecords[1].getPosition());
  std::cerr << "\n\n ----------------- \n";
  for (unsigned i = hashRecords[1].getPosition(); i < hashRecords[1].getPosition() + 25; ++i)
  {
    std::cerr << (unsigned int)referenceSequence.getBase(i) << " ";
  }
  std::cerr << "\n ---------------- \n" << std::endl;
#endif
}

TEST_F(HashtableFixture, DoItAll)
{
  //ASSERT_TRUE(false);
#if 0
  Mapper mapper;
  mapper.load(directory);
  HashtableConfig* config = mapper.get_hashtable_config();
  ReferenceGenome *reference = mapper.get_reference_genome();

  //  ResultChecker results;
  const int NUM_TEST_SEQUENCES = 10;
  for ( int i = 0; i < NUM_TEST_SEQUENCES; ++i ) {
    const bool insert_random_snp = false;

    // fetch_seed is a test hook -- remember that this isn't really part of
    // hashtable functionality
    const int seed_length = config->get_primary_seed_bases();
    seed_t seed = reference->fetch_seed(seed_length, i, insert_random_snp);

    HitList hits;
    mapper.get_hits(seed, hits);//, results);
#endif
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

HashtableFixture::HashtableConfig *HashtableFixture::hashtableConfig;
HashtableFixture::Hashtable *HashtableFixture::hashtable;
HashtableFixture::ReferenceSequence HashtableFixture::referenceSequence;
void HashtableFixture::SetUpTestCase()
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
