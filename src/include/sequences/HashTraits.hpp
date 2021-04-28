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

#ifndef ALIGN_HASH_TRAITS_HPP
#define ALIGN_HASH_TRAITS_HPP

#include <algorithm>

namespace dragenos {
namespace sequences {

//--------------------------------------------------------------------------------craczy
// We parameterize hashes by the following attribute class
//
class HashTraits {
public:
  //  enum { HASH_BUCKET_BYTES_LOG2 = 6 };
  //  enum { HASH_BUCKET_BYTES = (1 << HASH_BUCKET_BYTES_LOG2) };
  //  enum { HASH_RECORD_BYTES_LOG2 = 3 };
  //  enum { HASH_RECORD_BYTES = (1 << HASH_RECORD_BYTES_LOG2) };
  //  enum { HASH_MIN_EXTRA_BITS = 3 };
  //  enum { HASH_MAX_EXTRA_BITS = 3 };
  //  enum { MAX_PROBES = (1 << HASH_MAX_EXTRA_BITS) };
  //  enum { HASH_RECORDS_PER_BUCKET_LOG2 = (HASH_BUCKET_BYTES_LOG2 - HASH_RECORD_BYTES_LOG2) };
  //  enum { HASH_RECORDS_PER_BUCKET = (HASH_BUCKET_BYTES / HASH_RECORD_BYTES) };
  enum { HASH_RECORD_HASH_BITS = 23 };
  enum { HASH_RECORD_EXT_ID_BITS = 8 };
  //  enum { HASH_RECORD_THREAD_BITS = 6 };
  //  enum { HASH_ADDR_MINUS_ID_BITS = 16 };
  //  enum { HASH_THREADS_LOG2 = HASH_RECORDS_PER_BUCKET_LOG2 };
  //  enum { HASH_THREADS = HASH_RECORDS_PER_BUCKET };
  enum { SEC_CRC_BITS_MINUS_EXT_ID_HASH_BITS = 36 };
};

//class MAX_HASH_THREADS {
// public:
//  uint64_t operator() () {
//    return static_cast<uint64_t>(HashTraits::HASH_THREADS);
//  }
//};

}  // namespace sequences
}  // namespace dragenos

#endif  // #ifndef ALIGN_HASH_TRAITS_HPP
