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

// internal shared definitions

#pragma once

#include <immintrin.h>

#ifdef __GNUC__
#define LIKELY(x) __builtin_expect((x),1)
#define UNLIKELY(x) __builtin_expect((x),0)
#else
#define LIKELY(x) (x)
#define UNLIKELY(x) (x)
#endif

#define UNCLIP_BONUS 5
/* Convert the coordinate in the scoring matrix into the coordinate in one line of the band. */
#define set_u(u, w, i, j) { int x=(i)-(w); x=x>0?x:0; (u)=(j)-x+1; }

/* Convert the coordinate in the direction matrix into the coordinate in one line of the band. */
#define set_d(u, w, i, j, p) { int x=(i)-(w); x=x>0?x:0; x=(j)-x; (u)=x*3+p; }

/*! @function
  @abstract  Round an integer to the next closest power-2 integer.
  @param  x  integer to be rounded (in place)
  @discussion x will be modified.
 */
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))

typedef struct {
  uint16_t score;
  int32_t ref;   //0-based position
  int32_t read;    //alignment ending position on read, 0-based
} alignment_end;

struct _profile_sse2{
  __m128i* profile_byte;  // 0: none
  __m128i* profile_word;  // 0: none
  const int8_t* read;
  const int8_t* mat;
  int32_t readLen;
  int32_t n;
  uint8_t bias;
};

struct _profile_avx2{
  __m256i* profile_byte;  // 0: none
  __m256i* profile_word;  // 0: none
  const int8_t* read;
  const int8_t* mat;
  int32_t readLen;
  int32_t n;
  uint8_t bias;
};


