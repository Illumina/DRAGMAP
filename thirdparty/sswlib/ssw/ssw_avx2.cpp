/**
 ** DRAGEN Open Source Software
 ** Copyright (c) 2019-2020 Illumina, Inc.
 ** All rights reserved.
 **
 ** Based on SSW implementation
 ** https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library
 ** Version 0.1.4
 ** Last revision by Mengyao Zhao on 07/19/16 <zhangmp@bc.edu>
 **
 ** License: MIT
 ** Copyright (c) 2012-2015 Boston College
 ** Copyright (c) 2021 Illumina
 **
 ** Permission is hereby granted, free of charge, to any person obtaining a copy of this
 ** software and associated documentation files (the "Software"), to deal in the Software
 ** without restriction, including without limitation the rights to use, copy, modify,
 ** merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
 ** permit persons to whom the Software is furnished to do so, subject to the following
 ** conditions:
 ** The above copyright notice and this permission notice shall be included in all copies
 ** or substantial portions of the Software.
 ** THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
 ** INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
 ** PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 ** HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 ** OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 ** SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 **/

#include <immintrin.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <boost/assert.hpp>
#include "ssw.hpp"
#include "ssw_internal.hpp"

#ifdef __AVX2__

constexpr int AVX2_BYTE_ELEMS = 32;

static inline __attribute__((always_inline)) void* memalign_local(const size_t alignment, const size_t size)
{
  void* ptr = NULL;
  BOOST_ASSERT(0 == posix_memalign(&ptr, alignment, size));
  BOOST_ASSERT(ptr != nullptr);
  return ptr;
}

static inline int32_t getSegLen(const int32_t readLen)
{
  return (readLen + AVX2_BYTE_ELEMS - 1) / AVX2_BYTE_ELEMS;
}

/* Generate query profile rearrange query sequence & calculate the weight of match/mismatch. */
static __m256i* qP_byte_init (
    const int8_t* read_num,
    const int8_t* mat,
    const int32_t readLen,
    const int32_t n,	/* the edge length of the squre matrix mat */
    uint8_t bias) {

  int32_t segLen = getSegLen(readLen); /* Split the 128 bit register into 16 pieces.
								     Each piece is 8 bit. Split the read into 16 segments.
								     Calculate 16 segments in parallel.
   */
  __m256i* vProfile = (__m256i*)memalign_local(AVX2_BYTE_ELEMS, n * segLen * sizeof(__m256i));
  int8_t* t = (int8_t*)vProfile;
  int32_t nt, i, j, segNum;


  std::vector<int8_t> bonuses(segLen * AVX2_BYTE_ELEMS, bias);
  bonuses[0] = bias + UNCLIP_BONUS;
  bonuses[readLen-1] = bias + UNCLIP_BONUS;

  /* Generate query profile, rearrange query sequence & calculate the weight of match/mismatch */
  for (nt = 0; LIKELY(nt < n); nt ++) {
    for (i = 0; i < segLen; i ++) {
      j = i;

      for (segNum = 0; LIKELY(segNum < AVX2_BYTE_ELEMS) ; segNum ++) {
        *t++ =
          j>= readLen ?
              bias :
              mat[nt * n + read_num[j]] + bonuses[j];
        j += segLen;
      }
    }
  }

  return vProfile;
}

static __m256i* qP_byte_rev (
    const int8_t* read_num,
    const int8_t* mat,
    const int32_t readLen,
    const int32_t fullreadLen,
    const int32_t n,
    uint8_t bias) {
  //printf("qp_byte rev readLen %i fullreadLen %i \n",readLen,fullreadLen);

  constexpr int ELEMS = 32;

  int32_t segLen = getSegLen(readLen);
  __m256i* vProfile = (__m256i*)malloc(n * segLen * sizeof(__m256i));
  int8_t* t = (int8_t*)vProfile;
  int32_t nt, i, j, segNum;

  std::vector<int8_t> bonuses(segLen * AVX2_BYTE_ELEMS, bias);
  if (readLen == fullreadLen) {
    bonuses[0] = bias + UNCLIP_BONUS;
  }
  bonuses[readLen-1] = bias + UNCLIP_BONUS;

  for (nt = 0; LIKELY(nt < n); nt ++) {
    for (i = 0; i < segLen; i ++) {
      j = i;
      for (segNum = 0; LIKELY(segNum < ELEMS) ; segNum ++) {

        *t++ = j>= readLen ? bias : mat[nt * n + read_num[j]] + bonuses[j];
        j += segLen;
      }
    }
  }
  return vProfile;
}


static inline uint8_t max32fun(__m256i v)
{
  uint8_t* s = (uint8_t*)&v;
  uint8_t max = 0;
  // recognized by compiler
  for (int i = 0; i < 32; i++) {
    if (s[i] > max) {
      max = s[i];
    }
  }

  return max;
}

// shift left that handles 128-bit lane crossing
template<unsigned int N>
__m256i _mm256_shift_left(__m256i a)
{
  __m256i mask =  _mm256_srli_si256(
          _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0,0,3,0)), 16-N);
  return _mm256_or_si256(_mm256_slli_si256(a,N),mask);
}

/* Striped Smith-Waterman
   Record the highest score of each reference position.
   Return the alignment score and ending position of the best alignment, 2nd best alignment, etc.
   Gap begin and gap extension are different.
   wight_match > 0, all other weights < 0.
   The returned positions are 0-based.
 */
template<bool segLenIs5, int8_t ref_dir>
static void sw_avx2_byte (const int8_t* ref,
    int32_t refLen,
    int32_t readLen,
    const uint8_t weight_gapO, /* will be used as - */
    const uint8_t weight_gapE, /* will be used as - */
    const __m256i* vProfile,
    uint8_t terminate,	/* the best alignment score: used to terminate
												   the matrix calculation when locating the
												   alignment beginning point. If this score
												   is set to 0, it will not be used */
    uint8_t bias,  /* Shift 0 point to a positive value. */
    int32_t /*maskLen*/,
    alignment_end& best) {

  uint8_t max = 0;		                     /* the max alignment score */
  int32_t end_read = readLen - 1;
  int32_t end_ref = -1; /* 0_based best alignment ending point; Initialized as isn't aligned -1. */
  int32_t segLen;

  if (segLenIs5) {
    segLen = 5;
  }
  else {
    segLen = (readLen + AVX2_BYTE_ELEMS - 1) / AVX2_BYTE_ELEMS; /* number of segment */
  }

  /* Define 32 byte 0 vector. */
  __m256i vZero = _mm256_set1_epi32(0);
  size_t sz = segLen * sizeof(__m256i);
  __m256i* pvHStore = (__m256i*) memalign_local(AVX2_BYTE_ELEMS, sz);
  __m256i* pvHLoad = (__m256i*) memalign_local(AVX2_BYTE_ELEMS, sz);
  __m256i* pvE = (__m256i*) memalign_local(AVX2_BYTE_ELEMS, sz);
  __m256i* pvHmax = (__m256i*) memalign_local(AVX2_BYTE_ELEMS, sz);
  
  // no need to clear pvHLoad and pvHmax
  memset(pvHStore, 0, sz);
  memset(pvE, 0, sz);

  int32_t i, j;
  /* 32 byte insertion begin vector */
  __m256i vGapO = _mm256_set1_epi8(weight_gapO);

  /* 32 byte insertion extension vector */
  __m256i vGapE = _mm256_set1_epi8(weight_gapE);

  /* 32 byte bias vector */
  __m256i vBias = _mm256_set1_epi8(bias);


  __m256i vTerminate = _mm256_set1_epi8(terminate);

  __m256i vMaxScore = vZero; /* Trace the highest score of the whole SW matrix. */
  __m256i vMaxMark = vZero; /* Trace the highest score till the previous column. */
  __m256i vTemp;
  int32_t begin = 0, end = refLen, step = 1;

  /* outer loop to process the reference sequence */
  if (ref_dir == 1) {
    begin = refLen - 1;
    end = -1;
    step = -1;
  }
  // for reference
  for (i = begin; LIKELY(i != end); i += step) {
    int32_t cmp;
    __m256i e, vF = vZero, vMaxColumn = vZero; /* Initialize F value to 0.
							   Any errors to vH values will be corrected in the Lazy_F loop. */

    __m256i vH = _mm256_loadu_si256(&pvHStore[segLen - 1]);
    vH = _mm256_shift_left<1>(vH);

    const __m256i* vP = vProfile + ref[i] * segLen; /* Right part of the vProfile */

    /* Swap the 2 H buffers. */
    __m256i* pv = pvHLoad;
    pvHLoad = pvHStore;
    pvHStore = pv;

    /* inner loop to process the query sequence */
    for (j = 0; LIKELY(j < segLen); ++j) {

      __m256i p = _mm256_loadu_si256(vP + j);

      vH = _mm256_adds_epu8(vH, p);
      vH = _mm256_subs_epu8(vH, vBias); /* vH will be always > 0 */

      /* Get max from vH, vE and vF. */
      e = _mm256_loadu_si256(pvE + j); // originally _mm256_loadu_si256

      vH = _mm256_max_epu8(vH, e);
      vH = _mm256_max_epu8(vH, vF);

      vMaxColumn = _mm256_max_epu8(vMaxColumn, vH);

      /* Save vH values. */
      _mm256_storeu_si256(pvHStore + j, vH);

      /* Update vE value. */
      vH = _mm256_subs_epu8(vH, vGapO); /* saturation arithmetic, result >= 0 */

      e = _mm256_subs_epu8(e, vGapE);
      e = _mm256_max_epu8(e, vH);
      _mm256_storeu_si256(pvE + j, e);

      /* Update vF value. */
      vF = _mm256_subs_epu8(vF, vGapE);
      vF = _mm256_max_epu8(vF, vH);

      /* Load the next vH. */
      vH = _mm256_loadu_si256(pvHLoad + j);
    }

    /* Lazy_F loop: has been revised to disallow adjacent insertion and then deletion, so don't update E(i, j), learn from SWPS3 */
    /* reset pointers to the start of the saved data */
    j = 0;
    vH = _mm256_loadu_si256 (pvHStore + j);

    /*  the computed vF value is for the given column.  since */
    /*  we are at the end, we need to shift the vF value over */
    /*  to the next column. */
    vF = _mm256_shift_left<1>(vF);

    vTemp = _mm256_subs_epu8 (vH, vGapO);
    vTemp = _mm256_subs_epu8 (vF, vTemp);
    vTemp = _mm256_cmpeq_epi8 (vTemp, vZero);
    cmp  = _mm256_movemask_epi8 (vTemp);

    while ((unsigned int)cmp != 0xffffffff)
    {
      vH = _mm256_max_epu8 (vH, vF);
      vMaxColumn = _mm256_max_epu8(vMaxColumn, vH);
      _mm256_storeu_si256 (pvHStore + j, vH);
      vF = _mm256_subs_epu8 (vF, vGapE);
      j++;
      if (j >= segLen)
      {
        j = 0;
        vF = _mm256_shift_left<1>(vF);
      }
      vH = _mm256_loadu_si256 (pvHStore + j);

      vTemp = _mm256_subs_epu8 (vH, vGapO);
      vTemp = _mm256_subs_epu8 (vF, vTemp);
      vTemp = _mm256_cmpeq_epi8 (vTemp, vZero);
      cmp  = _mm256_movemask_epi8 (vTemp);
    }

    vMaxScore = _mm256_max_epu8(vMaxScore, vMaxColumn);
    vTemp = _mm256_cmpeq_epi8(vMaxMark, vMaxScore);
    cmp = _mm256_movemask_epi8(vTemp);
    if ((unsigned int)cmp != 0xffffffff) {
      uint8_t temp;
      vMaxMark = vMaxScore;
      temp = max32fun(vMaxScore);
      vMaxScore = vMaxMark;

      if (LIKELY(temp > max)) {
        max = temp;
        if (max + bias >= 255) break;	//overflow
        end_ref = i;

        /* Store the column with the highest alignment score in order to trace the alignment ending position on read. */
        for (j = 0; LIKELY(j < segLen); ++j) {
          __m256i v = _mm256_loadu_si256(&pvHStore[j]);
          _mm256_storeu_si256(&pvHmax[j], v);
        }
      }
    }

    vTemp = _mm256_cmpeq_epi8(vTerminate, vMaxColumn);
    cmp = _mm256_movemask_epi8(vTemp);
    if ((unsigned int)cmp != 0) {
      break;
    }
  }

  /* Trace the alignment ending position on read. */
  uint8_t *t = (uint8_t*)pvHmax;
  int32_t column_len = segLen * AVX2_BYTE_ELEMS;
  for (i = 0; LIKELY(i < column_len); ++i, ++t) {
    int32_t temp;
    if (*t == max) {
      temp = i / AVX2_BYTE_ELEMS + i % AVX2_BYTE_ELEMS * segLen;
      if (temp < end_read) end_read = temp;
    }
  }

  free(pvHmax);
  free(pvE);
  free(pvHLoad);
  free(pvHStore);

  /* Record the best alignment. */
  best.score = max + bias >= 255 ? 255 : max;
  best.ref = end_ref;
  best.read = end_read;
}

s_profile_avx2* ssw_init_avx2 (
    const int8_t* read, const int32_t readLen,
    const int8_t* mat, const int32_t n, const int32_t bias,
    const int8_t score_size)
{
  BOOST_ASSERT(score_size == 0 &&
      "The AVX2 SSW variant initializes profile only for 8-bit scoring, "
      "16-bit scoring profile is initialized only when needed."
  );

  s_profile_avx2* p = (s_profile_avx2*)calloc(1, sizeof(struct _profile_avx2));
  p->profile_byte = 0;
  p->profile_word = 0;
  p->bias = bias;

  p->profile_byte = qP_byte_init (read, mat, readLen, n, bias);

  p->read = read;
  p->mat = mat;
  p->readLen = readLen;
  p->n = n;
  return p;
}

void init_destroy_avx2 (s_profile_avx2* p) {
  if (p != nullptr) {
    free(p->profile_byte);
    free(p->profile_word);
    free(p);
  }
}

// free
static int8_t* seq_reverse(const int8_t* seq, int32_t end)  /* end is 0-based alignment ending position */
{
  int8_t* reverse = (int8_t*)calloc(end + 1, sizeof(int8_t));
  int32_t start = 0;
  while (LIKELY(start <= end)) {
    reverse[start] = seq[end];
    reverse[end] = seq[start];
    ++ start;
    -- end;
  }
  return reverse;
}


s_align* ssw_align_avx2 (
    const s_profile_avx2* prof,
    const int8_t* ref,
    int32_t refLen,
    const uint8_t weight_gapO,
    const uint8_t weight_gapE,
    const uint8_t flag,	//  (from high to low) bit 5: return the best alignment beginning position; 6: if (ref_end1 - ref_begin1 <= filterd) && (read_end1 - read_begin1 <= filterd), return cigar; 7: if max score >= filters, return cigar; 8: always return cigar; if 6 & 7 are both setted, only return cigar when both filter fulfilled
    const uint16_t filters,
    const int32_t filterd,
    const int32_t maskLen) {
  //printf("----- ssw_align readLen %i   read %p ------\n",prof->readLen,prof->read);
  alignment_end best, best_reverse;
  int32_t band_width = 0, readLen = prof->readLen;
  int8_t* read_reverse = 0;
  cigar* path;
  s_align* r = (s_align*)calloc(1, sizeof(s_align));
  r->ref_begin1 = -1;
  r->read_begin1 = -1;
  r->cigar = 0;
  r->cigarLen = 0;

  int32_t segLen = (readLen + AVX2_BYTE_ELEMS - 1) / AVX2_BYTE_ELEMS;

  // Find the alignment scores and ending positions
  if (prof->profile_byte) {

    if (segLen == 5) {
      sw_avx2_byte<true, 0>(ref, refLen, readLen, weight_gapO, weight_gapE, prof->profile_byte, -1, prof->bias, maskLen, best);
    }
    else {
      sw_avx2_byte<false, 0>(ref, refLen, readLen, weight_gapO, weight_gapE, prof->profile_byte, -1, prof->bias, maskLen, best);
    }

    if (best.score == 255) {
      free(r);

      // use SSW implementation as a fallback when 8-bit scoring overflowed,
      // does not happen with ~150 base reads
      // there is no overflow check in SSW for the 16-bit scoring
      s_profile_sse2* prof_sse2 = ssw_init_sse2 (
          prof->read, prof->readLen,
          prof->mat, prof->n, prof->bias, 2);

      s_align* r16 = ssw_align_sse2(
          prof_sse2, ref, refLen,
          weight_gapO, weight_gapE,
          flag, filters, filterd, maskLen);

      init_destroy_sse2(prof_sse2);
      return r16;
    }
  }else {
    fprintf(stderr, "Please call the function ssw_init before ssw_align.\n");
    free(r);
    return nullptr;
  }
  r->score1 = best.score;
  r->ref_end1 = best.ref;
  r->read_end1 = best.read;
  r->score2 = 0;
  r->ref_end2 = -1;

  if (flag == 0 || (flag == 2 && r->score1 < filters)) {
    return r;
  }

  // Find the beginning position of the best alignment.
  read_reverse = seq_reverse(prof->read, r->read_end1);
  __m256i* vP = qP_byte_rev(read_reverse, prof->mat, r->read_end1 + 1,readLen, prof->n, prof->bias);
  sw_avx2_byte<false, 1>(ref, r->ref_end1 + 1, r->read_end1 + 1, weight_gapO, weight_gapE, vP, r->score1, prof->bias, maskLen, best_reverse);
  free(vP);
  free(read_reverse);

  r->ref_begin1 = best_reverse.ref;
  r->read_begin1 = r->read_end1 - best_reverse.read;

  if ((7&flag) == 0 || ((2&flag) != 0 && r->score1 < filters) || ((4&flag) != 0 && (r->ref_end1 - r->ref_begin1 > filterd || r->read_end1 - r->read_begin1 > filterd))) {
    return r;
  }

  // Generate cigar.
  refLen = r->ref_end1 - r->ref_begin1 + 1;
  readLen = r->read_end1 - r->read_begin1 + 1;
  band_width = abs(refLen - readLen) + 1;

  path = banded_sw(ref + r->ref_begin1, prof->read + r->read_begin1, refLen, readLen, r->score1, weight_gapO, weight_gapE, band_width, prof->mat, prof->n,  r->read_begin1 , prof->readLen);
  if (path == nullptr) {
    free(r);
    r = nullptr;
  }
  else {
    r->cigar = path->seq;
    r->cigarLen = path->length;
    free(path);
  }

  return r;
}

#endif // #ifdef __AVX2__
