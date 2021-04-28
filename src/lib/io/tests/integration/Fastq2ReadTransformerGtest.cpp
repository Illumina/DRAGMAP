#include "gtest/gtest.h"

#include <algorithm>

#include "align/Aligner.hpp"
#include "io/Fastq2ReadTransformer.hpp"

const std::string THREE_RECORDS(
    "@blah\n"
    "ACGT\n"
    "+\n"
    "#AAA\n"

    "@NB551322:14:HFVLLBGX9:4:11401:24054:1050 2:N:0:CGGCTATG+CCGTCGCC\n"
    "NTGTCGGGGCAGGCAGGGCTCCTCGGGCAGCGGCTCATGAGAGAAGACGGAATCCTCCCCTGAGGAGCACGTAGAGCTCCGGGTGTCGGGAAAGCTGGGGG\n"
    "+\n"
    "#AAAAEEEEEEEEEEEEEEEEEEEEEEE6EEEEEEEEEEEEEEEEEAEAEEEEEEEEEEEEAEEE/EEEEAEEEAEEEEEEEEEEEEEEEEEEEEE/EEEE\n"

    "@NB551322:14:HFVLLBGX9:4:11401:9125:1052 2:N:0:CGGCTATG+CCGTCGCC\n"
    "NTGCCTCGCACATAGCACCTGCCCCTGTGCAGCCCCCCATGATTCGGCG\n"
    "+\n"
    "#AAAAEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE");

namespace dragenos {
namespace io {

struct TagNameBufferDump {
  std::string name_;
};

template <typename DumpT>
DumpT& dump(DumpT& dump, const io::FastqToReadTransformer& transformer)
{
  dump.name_ = std::string(transformer.tmpName_.begin(), transformer.tmpName_.end());
  return dump;
}

}  // namespace io
}  // namespace dragenos

using namespace dragenos;

TEST(Fastq2ReadTransformer, ProperMemoryHandover)
{
  std::stringstream stm(THREE_RECORDS);

  fastq::Tokenizer           t(stm, 1024);
  io::FastqToReadTransformer fastq2Read;
  align::Aligner::Read       read;

  ASSERT_EQ(true, t.next());
  fastq2Read(t.token(), 0, 1, read);

  ASSERT_EQ(true, t.next());
  fastq2Read(t.token(), 0, 1, read);

  io::TagNameBufferDump name;
  io::dump(name, fastq2Read);
  // The implementation ensures that memory buffers are being passed arround
  // without being destroyed.
  // After second conversion, the initial buffer must come back to the token
  ASSERT_EQ(std::string("blah"), name.name_);
}
