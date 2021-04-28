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

#ifndef REFERENCE_HASHTABLE_TRAITS_HPP
#define REFERENCE_HASHTABLE_TRAITS_HPP

#include <algorithm>

#include "sequences/HashTraits.hpp"

namespace dragenos {
namespace reference {

//--------------------------------------------------------------------------------adam
// We parameterize hashtables by the following attribute class
//
class HashtableTraits {
public:
  constexpr static unsigned BUCKET_RECORDS_LOG2          = 3;  // 8=2^3 Hash Records per Bucket
  constexpr static unsigned HASH_BUCKET_BYTES_LOG2       = 6;
  constexpr static unsigned HASH_BUCKET_BYTES            = (1 << HASH_BUCKET_BYTES_LOG2);
  constexpr static unsigned HASH_RECORD_BYTES_LOG2       = 3;
  constexpr static unsigned HASH_RECORD_BYTES            = (1 << HASH_RECORD_BYTES_LOG2);
  constexpr static unsigned HASH_MIN_EXTRA_BITS          = 3;
  constexpr static unsigned HASH_MAX_EXTRA_BITS          = 3;
  constexpr static unsigned MAX_PROBES                   = (1 << HASH_MAX_EXTRA_BITS);
  constexpr static unsigned HASH_RECORDS_PER_BUCKET_LOG2 = (HASH_BUCKET_BYTES_LOG2 - HASH_RECORD_BYTES_LOG2);
  constexpr static unsigned HASH_RECORDS_PER_BUCKET      = (HASH_BUCKET_BYTES / HASH_RECORD_BYTES);
  constexpr static unsigned HASH_RECORD_HASH_BITS        = sequences::HashTraits::HASH_RECORD_HASH_BITS;
  //  constexpr static unsigned HASH_RECORD_EXT_ID_BITS = 8;
  constexpr static unsigned HASH_RECORD_THREAD_BITS = 6;
  constexpr static unsigned HASH_ADDR_MINUS_ID_BITS = 16;
  constexpr static unsigned HASH_THREADS_LOG2       = HASH_RECORDS_PER_BUCKET_LOG2;
  constexpr static unsigned HASH_THREADS            = HASH_RECORDS_PER_BUCKET;
  //  constexpr static unsigned SEC_CRC_BITS_MINUS_EXT_ID_HASH_BITS = 36;
  constexpr static unsigned MAX_WRAP_BYTES_LOG2 = 15;
};

class MAX_HASH_THREADS {
public:
  uint64_t operator()() { return static_cast<uint64_t>(HashtableTraits::HASH_THREADS); }
};

}  // namespace reference
}  // namespace dragenos

#endif  // #ifndef REFERENCE_HASHTABLE_TRAITS_HPP
