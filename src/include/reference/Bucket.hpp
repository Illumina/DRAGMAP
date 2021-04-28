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

#ifndef REFERENCE_BUCKET_HPP
#define REFERENCE_BUCKET_HPP

#include <array>
#include "reference/HashRecord.hpp"
#include "reference/HashtableTraits.hpp"

namespace dragenos {
namespace reference {

template <size_t N_LOG2>
class BucketT : public std::array<HashRecord, 1 << N_LOG2> {
public:
  static constexpr size_t hashRecordCount = (1 << N_LOG2);
};

enum { HASH_BUCKET_BYTES_LOG2 = 6 };
enum { HASH_BUCKET_BYTES = (1 << HASH_BUCKET_BYTES_LOG2) };
enum { HASH_RECORD_BYTES_LOG2 = 3 };

typedef BucketT<HashtableTraits::BUCKET_RECORDS_LOG2> Bucket;

}  // namespace  reference
}  // namespace dragenos

#endif  // #ifndef REFERENCE_BUCKET_HPP
