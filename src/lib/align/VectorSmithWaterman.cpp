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

#include "align/VectorSmithWaterman.hpp"
#include <sstream>
#include <vector>
#include "ssw/ssw.hpp"

namespace dragenos {
namespace align {

void VectorSmithWaterman::destroyReadContext(int readIdx)
{
#ifdef __AVX2__
  init_destroy_avx2(profile_[readIdx]);
  init_destroy_avx2(profileRev_[readIdx]);
#else
  init_destroy_sse2(profile_[readIdx]);
  init_destroy_sse2(profileRev_[readIdx]);
#endif

  profile_[readIdx]    = NULL;
  profileRev_[readIdx] = NULL;
}

void VectorSmithWaterman::initReadContext(
    const unsigned char* queryBegin, const unsigned char* queryEnd, int readIdx)
{
  querySize_[readIdx] = std::distance(queryBegin, queryEnd);
  query_[readIdx].resize(querySize_[readIdx]);
  queryRev_[readIdx].resize(querySize_[readIdx]);

  std::copy(queryBegin, queryEnd, query_[readIdx].begin());

  const std::reverse_iterator<const unsigned char*> rbegin(queryEnd);
  const std::reverse_iterator<const unsigned char*> rend(queryBegin);
  std::copy(rbegin, rend, queryRev_[readIdx].begin());

  destroyReadContext(readIdx);

  const int8_t* queryBeginInt    = (int8_t*)query_[readIdx].data();
  const int8_t* queryRevBeginInt = (int8_t*)queryRev_[readIdx].data();

#ifdef __AVX2__
  // AVX2 variant initializes profile only for 8-bit scoring (last argument 0),
  // 16-bit scoring is initialized only when needed
  profile_[readIdx] =
      ssw_init_avx2(queryBeginInt, querySize_[readIdx], sswScoringMat_, sswAlphabetSize_, sswBias_, 0);
  profileRev_[readIdx] =
      ssw_init_avx2(queryRevBeginInt, querySize_[readIdx], sswScoringMat_, sswAlphabetSize_, sswBias_, 0);
#else
  profile_[readIdx] =
      ssw_init_sse2(queryBeginInt, querySize_[readIdx], sswScoringMat_, sswAlphabetSize_, sswBias_, 2);
  profileRev_[readIdx] =
      ssw_init_sse2(queryRevBeginInt, querySize_[readIdx], sswScoringMat_, sswAlphabetSize_, sswBias_, 2);
#endif
}

// returns alignment score
// returns operations list in cigar
uint16_t VectorSmithWaterman::align(
    const unsigned char* /*queryBegin*/,
    const unsigned char* /*queryEnd*/,
    const unsigned char* databaseBegin,
    const unsigned char* databaseEnd,
    bool                 reverseQuery,
    std::string&         cigar,
    int                  readIdx)
{
  const int8_t* databaseBeginInt = (int8_t*)databaseBegin;
  const int8_t* databaseEndInt   = (int8_t*)databaseEnd;

  uint16_t score = 0;

  // const int querySize = std::distance(queryBeginInt, queryEndInt);
  const int dbSize = std::distance(databaseBeginInt, databaseEndInt);

#ifdef __AVX2__
  s_profile_avx2* profile;
#else
  s_profile_sse2* profile;
#endif

  // use the already built profile
  int querySize = querySize_[readIdx];
  if (reverseQuery) {
    profile = profileRev_[readIdx];
  } else {
    profile = profile_[readIdx];
  }

  s_align* result;

  uint8_t flag = 0;
  //flag |= 0x08;  // report ref position
  //flag |= 0x0F;  // report cigar
  // 1 - always compute cigar; 1 << 5 - get just the start and end positions
  flag             = 1;
  uint16_t filters = 0;
  int32_t  filterd = 0;
  int32_t  maskLen = querySize / 2;

  result =
#ifdef __AVX2__
      ssw_align_avx2(
#else
      ssw_align_sse2(
#endif
          profile, databaseBeginInt, dbSize, gapInit_, gapExtend_, flag, filters, filterd, maskLen);

  this->getCigarOperations(*result, querySize, cigar);

  int softClipStart = result->read_begin1;
  int softClipEnd   = querySize - result->read_end1 - 1;
#ifdef TRACE_VECTOR_SMITH_WATERMAN

  printf(
      "convert_cigar cig len %i score %i ref %i -- %i \n",
      result->cigarLen,
      result->score1,
      result->ref_begin1,
      result->ref_end1);

  std::string cigarres = convert_cigar(*result, querySize);
  printf("VEC SW cigar %s \n", cigarres.c_str());

#endif

  score = result->score1;

  align_destroy(result);
  uint16_t unclipScoreAdjsutment = (softClipStart ? 0 : unclipScore_) + (softClipEnd ? 0 : unclipScore_);
  unclipScoreAdjsutment          = std::min(unclipScoreAdjsutment, score);
  return score - unclipScoreAdjsutment;
}

// converts to dos-like operations list
//

void VectorSmithWaterman::getCigarOperations(
    const s_align& s_al, const int& query_len, std::string& operations)
{
  operations.clear();

  if (s_al.ref_begin1 > 0) {
    operations.resize(operations.size() + s_al.ref_begin1, 'N');
  }

  if (s_al.cigarLen > 0) {
    if (s_al.read_begin1 > 0) {
      operations.resize(operations.size() + s_al.read_begin1, 'S');
    }

    for (int i = 0; i < s_al.cigarLen; ++i) {
      operations.resize(operations.size() + cigar_int_to_len(s_al.cigar[i]), cigar_int_to_op(s_al.cigar[i]));
    }

    int end = query_len - s_al.read_end1 - 1;
    if (end > 0) {
      operations.resize(operations.size() + end, 'S');
    }
  }
}

// returns traditional human readable cigar string
std::string VectorSmithWaterman::convert_cigar(const s_align& s_al, const int& query_len)
{
  std::string result_cigar_string;
  if (s_al.cigarLen > 0) {
    std::ostringstream cigar_string;
    if (s_al.read_begin1 > 0) {
      //      uint32_t cigar = to_cigar_int(s_al.read_begin1, 'S');
      cigar_string << s_al.read_begin1 << 'S';
    }

    for (int i = 0; i < s_al.cigarLen; ++i) {
      cigar_string << cigar_int_to_len(s_al.cigar[i]) << cigar_int_to_op(s_al.cigar[i]);
    }

    int end = query_len - s_al.read_end1 - 1;
    if (end > 0) {
      //      uint32_t cigar = to_cigar_int(end, 'S');
      cigar_string << end << 'S';
    }

    result_cigar_string = cigar_string.str();
  }  // end if
  return result_cigar_string;
}

}  // namespace align
}  // namespace dragenos
