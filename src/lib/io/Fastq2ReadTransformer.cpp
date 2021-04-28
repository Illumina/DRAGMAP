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

#include "io/Fastq2ReadTransformer.hpp"

namespace dragenos {
namespace io {

#if VECTOR_REGISTER_WIDTH

#if 0
//__attribute__((noinline))
// vectorized implementation that requires pshufb
void FastqToReadTransformer::convertBases(Base* begin, const std::size_t size)
{
  assert(0 == size % VECTOR_REGISTER_WIDTH);

  const Base* end = begin + size;
  for (Base* b = begin; end != b; b += VECTOR_REGISTER_WIDTH)
  {
    __m128i input = _mm_loadu_si128(reinterpret_cast<const __m128i*>(b));
    // convert to uppercase
//    input = _mm_and_si128(input, _mm_set1_epi8(0xdf));
    for (size_t i = 0; VECTOR_REGISTER_WIDTH != i; ++i)
    {
      (((unsigned char*)&input)[i]) &= 0xdf;
    }
//    std::cout << "caps : " << input << std::endl;

//    const __m128i less = _mm_cmpgt_epi8(input, _mm_set1_epi8('N'));
    __m128i less = input;
    for (size_t i = 0; VECTOR_REGISTER_WIDTH != i; ++i)
    {
      (((unsigned char*)&less)[i]) = (((unsigned char*)&less)[i]) > 'N' ? 0xff : 0x00;
    }
//    std::cout << ">'N' : " << less << std::endl;

    // T needs a special adjustment as it is further than 15 positions away from A
//    const __m128i adj = _mm_and_si128(less, _mm_set1_epi8('T' - 'N' - 1));
    __m128i adj = less;
    for (size_t i = 0; VECTOR_REGISTER_WIDTH != i; ++i)
    {
      (((unsigned char*)&adj)[i]) &= (unsigned char)('T' - 'N' - 1);
    }
//    std::cout << "adj  : " << adj << std::endl;

//    const __m128i tmp = _mm_subs_epu8(input, _mm_set1_epi8('A'));
    __m128i tmp = input;
    for (size_t i = 0; VECTOR_REGISTER_WIDTH != i; ++i)
    {
      (((unsigned char*)&tmp)[i]) -= 'A';
    }
//    std::cout << "-'A' : " << tmp << std::endl;

//    const __m128i idx = _mm_subs_epu8(tmp, adj);
    __m128i idx = tmp;
    for (size_t i = 0; VECTOR_REGISTER_WIDTH != i; ++i)
    {
      (((unsigned char*)&idx)[i]) -= (((unsigned char*)&adj)[i]);
    }
//    std::cout << "idx  : " << idx << std::endl;

    static const char unused = 0;
    // The lookup table.
//    static const __m128i LUT = _mm_setr_epi8(
//      0/*A*/, unused, 1/*C*/, unused, unused, unused, 2/*G*/, unused,
//      unused, unused, unused, unused, unused, 2/*N*/,
//      3/*T*/,
//      unused
//    );

    static const char LUT[VECTOR_REGISTER_WIDTH] = {
      0/*A*/, unused, 1/*C*/, unused, unused, unused, 2/*G*/, unused,
      unused, unused, unused, unused, unused, 2/*N*/, 3/*T*/, unused
    };

    // unfortunately there is no way to autovectorize into pshufb
    const __m128i result = _mm_shuffle_epi8(*reinterpret_cast<const __m128i*>(LUT), idx);
//    __m128i result;
//    for (size_t i = 0; VECTOR_REGISTER_WIDTH != i; ++i)
//    {
//      (((unsigned char*)&result)[i]) = LUT[((unsigned char*)&idx)[i] & 0x0f];
//    }
//    std::cout << "rslt : " << result << std::endl;

    *reinterpret_cast<__m128i*>(b) = result;
  }

}

#endif  // 0

struct VectorOfChars {
  char        v_[VECTOR_REGISTER_WIDTH];
  char&       operator[](std::size_t i) { return v_[i]; }
  const char& operator[](std::size_t i) const { return v_[i]; }
};

//__attribute__((noinline))
void FastqToReadTransformer::convertBases2(sequences::Read::Base* begin, const std::size_t size)
{
  assert(0 == size % VECTOR_REGISTER_WIDTH);

  const sequences::Read::Base* end = begin + size;
  for (sequences::Read::Base* b = begin; end != b; b += VECTOR_REGISTER_WIDTH) {
    //__m128i input = _mm_loadu_si128(reinterpret_cast<const __m128i*>(b));
    VectorOfChars input = *reinterpret_cast<VectorOfChars*>(b);
    // convert to uppercase
    //    input = _mm_and_si128(input, _mm_set1_epi8(0xdf));
    for (size_t i = 0; VECTOR_REGISTER_WIDTH != i; ++i) {
      input[i] &= 0xdf;
    }
    //    std::cout << "caps : " << input << std::endl;

    VectorOfChars n;
    for (size_t i = 0; VECTOR_REGISTER_WIDTH != i; ++i) {
      //(((unsigned char*)&n)[i]) = input[i] == 'N' ? 0xF : 0x00;
      (((unsigned char*)&n)[i]) = input[i] == 'N' ? 0x0 : 0x00;
    }
    //    std::cout << "n : " << n << std::endl;
    VectorOfChars a;
    for (size_t i = 0; VECTOR_REGISTER_WIDTH != i; ++i) {
      a[i] = input[i] == 'A' ? 0x01 : 0x00;
    }
    VectorOfChars c;
    for (size_t i = 0; VECTOR_REGISTER_WIDTH != i; ++i) {
      c[i] = input[i] == 'C' ? 0x02 : 0x00;
    }
    //    std::cout << "c : " << n << std::endl;

    VectorOfChars g;
    for (size_t i = 0; VECTOR_REGISTER_WIDTH != i; ++i) {
      g[i] = input[i] == 'G' ? 0x04 : 0x00;
    }
    //    std::cout << "g : " << n << std::endl;

    VectorOfChars t;
    for (size_t i = 0; VECTOR_REGISTER_WIDTH != i; ++i) {
      t[i] = input[i] == 'T' ? 0x08 : 0x00;
    }
    //    std::cout << "t : " << n << std::endl;

    VectorOfChars result = n;
    //    __m128i result = n;
    for (size_t i = 0; VECTOR_REGISTER_WIDTH != i; ++i) {
      result[i] |= a[i];
    }
    for (size_t i = 0; VECTOR_REGISTER_WIDTH != i; ++i) {
      result[i] |= c[i];
    }
    for (size_t i = 0; VECTOR_REGISTER_WIDTH != i; ++i) {
      result[i] |= g[i];
    }
    for (size_t i = 0; VECTOR_REGISTER_WIDTH != i; ++i) {
      result[i] |= t[i];
    }
    //    std::cout << "rslt : " << result << std::endl;

    *reinterpret_cast<VectorOfChars*>(b) = result;
  }
}

// !!! must be called before convertBases
//__attribute__((noinline))
void FastqToReadTransformer::convertQualities(
    const sequences::Read::Base* b, sequences::Read::Qscore* q, const std::size_t size)
{
  assert(0 == size % VECTOR_REGISTER_WIDTH);
  //  Qscore qmask[VECTOR_REGISTER_WIDTH];
  const sequences::Read::Qscore* end = q + size;
  for (; end != q; q += VECTOR_REGISTER_WIDTH, b += VECTOR_REGISTER_WIDTH) {
    //    __m128i bases = _mm_loadu_si128(reinterpret_cast<const __m128i*>(b));
    const VectorOfChars bases = *reinterpret_cast<const VectorOfChars*>(b);
    //    std::cout << "bazs : " << bases << std::endl;
    //    const __m128i n = _mm_cmpeq_epi8(bases, _mm_set1_epi8('N'));
    VectorOfChars n;
    for (size_t i = 0; VECTOR_REGISTER_WIDTH != i; ++i) {
      (((unsigned char*)&n)[i]) = bases[i] == 'N' ? 0xff : 0x00;
    }
    //    std::cout << "n    : " << n << std::endl;

    //    const __m128i input = _mm_loadu_si128(reinterpret_cast<const __m128i*>(q));
    VectorOfChars input = *reinterpret_cast<const VectorOfChars*>(q);
    //    std::cout << "inpt : " << input << std::endl;

    for (size_t i = 0; VECTOR_REGISTER_WIDTH != i; ++i) {
      input[i] -= q0_[i];
    }

    //    const __m128i result = _mm_andnot_si128(n, input);
    VectorOfChars result;
    for (size_t i = 0; VECTOR_REGISTER_WIDTH != i; ++i) {
      result[i] = ~n[i] & input[i];
    }

    input = *reinterpret_cast<const VectorOfChars*>(q2_);
    for (size_t i = 0; VECTOR_REGISTER_WIDTH != i; ++i) {
      input[i] = n[i] & input[i];
    }

    for (size_t i = 0; VECTOR_REGISTER_WIDTH != i; ++i) {
      result[i] = result[i] | input[i];
    }
    //    std::cout << "rslt : " << result << std::endl;
    *reinterpret_cast<VectorOfChars*>(q) = result;
  }
}

void FastqToReadTransformer::convertBasesAndQualities()
{
  assert(tmpBases_.size() == tmpQscores_.size());
  const std::size_t paddedSize =
      ((tmpQscores_.size() + VECTOR_REGISTER_WIDTH - 1) / VECTOR_REGISTER_WIDTH) * VECTOR_REGISTER_WIDTH;
  convertQualities(&tmpBases_.front(), &tmpQscores_.front(), paddedSize);
  convertBases2(&tmpBases_.front(), paddedSize);
}

#else  // VECTOR_REGISTER_WIDTH
void FastqToReadTransformer::convertBasesAndQualities(const char q0)
{
  assert(bases_.size() == qualities_.size());
  // convert alphabetic to numeric -- obscure implementation to enable auto-vectorization
  for (size_t i = 0; bases_.size() > i; ++i) {
    Base&   c = bases_[i];
    Qscore& q = qualities_[i];
    q -= q0;
    const char isC  = (c == 'C');
    const char isN  = (c == 'N');
    const char notN = (c != 'N');
    const char isG  = (c == 'G') | isN;
    const char isT  = (c == 'T');
    c               = isC + (isG << 1) + (isT << 1) + isT;
    q *= notN;
  }
}

#endif  // VECTOR_REGISTER_WIDTH

}  // namespace io
}  // namespace dragenos
