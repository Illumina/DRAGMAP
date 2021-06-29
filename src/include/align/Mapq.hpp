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

#ifndef ALIGN_MAPQ_HPP
#define ALIGN_MAPQ_HPP

#include <math.h>

#include "align/Score.hpp"
#include "common/DragenLogger.hpp"

namespace dragenos {
namespace align {

typedef int           MapqType;
static const MapqType MAPQ_MAX    = 60;
static const MapqType HW_MAPQ_MAX = 250;

//  static const double MAPQ_COEFF = 20.8;
//  static const double MAPQ_COEFF = 88.28125;
static const double MAPQ_COEFF   = 38912 >> 8;  //152;
static const double MAPQ_COEFF_I = 38912;       //152;
//  static const double MAPQ_COEFF = 139.376;
//  static const double MAPQ_COEFF = 99.27;
//

inline double mapqCoeffScaled(ScoreType snpCost)
{
  return MAPQ_COEFF * (5.0 / snpCost);
}

inline int mapqCoeffScaledI(ScoreType snpCost)
{
  return MAPQ_COEFF_I * 5 / snpCost;
}

inline int log2_approx(const int x)
{
  static const int log2_approx_c[] = {
      0b0000000, 0b0000001, 0b0000011, 0b0000100, 0b0000110, 0b0000111, 0b0001000, 0b0001010, 0b0001011,
      0b0001101, 0b0001110, 0b0001111, 0b0010001, 0b0010010, 0b0010011, 0b0010100, 0b0010110, 0b0010111,
      0b0011000, 0b0011010, 0b0011011, 0b0011100, 0b0011101, 0b0011111, 0b0100000, 0b0100001, 0b0100010,
      0b0100011, 0b0100101, 0b0100110, 0b0100111, 0b0101000, 0b0101001, 0b0101010, 0b0101100, 0b0101101,
      0b0101110, 0b0101111, 0b0110000, 0b0110001, 0b0110010, 0b0110011, 0b0110100, 0b0110101, 0b0110111,
      0b0111000, 0b0111001, 0b0111010, 0b0111011, 0b0111100, 0b0111101, 0b0111110, 0b0111111, 0b1000000,
      0b1000001, 0b1000010, 0b1000011, 0b1000100, 0b1000101, 0b1000110, 0b1000111, 0b1001000, 0b1001001,
      0b1001010, 0b1001011, 0b1001100, 0b1001101, 0b1001110, 0b1001111, 0b1010000, 0b1010001, 0b1010001,
      0b1010010, 0b1010011, 0b1010100, 0b1010101, 0b1010110, 0b1010111, 0b1011000, 0b1011001, 0b1011010,
      0b1011011, 0b1011011, 0b1011100, 0b1011101, 0b1011110, 0b1011111, 0b1100000, 0b1100001, 0b1100001,
      0b1100010, 0b1100011, 0b1100100, 0b1100101, 0b1100110, 0b1100111, 0b1100111, 0b1101000, 0b1101001,
      0b1101010, 0b1101011, 0b1101011, 0b1101100, 0b1101101, 0b1101110, 0b1101111, 0b1101111, 0b1110000,
      0b1110001, 0b1110010, 0b1110011, 0b1110011, 0b1110100, 0b1110101, 0b1110110, 0b1110110, 0b1110111,
      0b1111000, 0b1111001, 0b1111001, 0b1111010, 0b1111011, 0b1111100, 0b1111100, 0b1111101, 0b1111110,
      0b1111111, 0b1111111};

  // Can’t take log of zero.  HW can detect zero input and decide what to do; using a zero result is common.
  //  if (!x) {
  //    return 0;
  //  }

  // Find the integer portion of the log.
  // HW does this by detecting the position of the most significant ‘1’ bit.
  int logInt = 0;
  for (int tmp = x; tmp >>= 1; ++logInt)
    ;

  // Normalize X into [1,2).  HW does this by bit-shifting.
  int norm = (x << 7) >> logInt;

  // Truncate to 7 fractional bits.
  // (Also practical to round to 7 bits, but I don’t think HW does this anywhere.)
  // These 7 fractional bits are the input to HW’s lookup table.
  // High-precision log of the truncated normalized value, within [0,1)
  int logFrac = log2_approx_c[norm & 0x7f];
  // Round to 7 fractional bits.  These 7 bits are the output of HW’s lookup table.
  // Combine integer and fractional portions
  return (logInt << 7) + logFrac;
}

inline int log2_simple(int x)
{
  int logInt = 0;
  for (int tmp = x; tmp >>= 1; ++logInt)
    ;

  //  -- Fractional portion of log2 approximation: 7 bits below the most significant bit
  const int count_log_shift_v = ((x << 7) >> logInt) & 0x7f;
  const int sub_count_log2    = (logInt << 7) | count_log_shift_v;

  return sub_count_log2;
}

inline MapqType aln2mapq(ScoreType snpCost, const double read_len_avg)
{
  const int    log2_length = log2_approx(read_len_avg);
  const double aln2mapq    = mapqCoeffScaled(snpCost) / ((log2_length * log2_length) >> 7);
  return aln2mapq * (1 << 20);
}

inline int mapq2aln(ScoreType snpCost, const double readLength)
{
  const int log2_length = log2_approx(readLength);
  const int mapq2aln    = (log2_length * log2_length / (mapqCoeffScaledI(snpCost) >> 4));
  return mapq2aln;
}

inline MapqType computeMapq(ScoreType snpCost, const ScoreType as, const ScoreType xs, const double n1)
{
  const ScoreType s1 = as;
  const ScoreType s2 = xs;
  //  const MapqType mapq = std::min<MapqType>(
  //    MAPQ_COEFF * log2(n1) * (s1 - s2) / n1, MAPQ_MAX);
  //  const double R = 1.0 - (MaxPossibleScore- as) / (snpCost * n1);
  //  const MapqType mapq = std::min<MapqType>(
  //    (s1 - s2) * MAPQ_COEFF / log2(n1) / log2(n1), MAPQ_MAX);

  //  const MapqType mapq = std::min<MapqType>((s1 - s2) * aln2mapq(snpCost, n1), MAPQ_MAX);
  //  const MapqType mapq = ((s1 - s2) * aln2mapq(snpCost, n1)) >> 7;  // >> 13;
  const int      a2m_scale = aln2mapq(snpCost, n1);
  const MapqType mapq      = ((s1 - s2) * a2m_scale) >> 13;

#ifdef TRACE_SCORING
  std::cerr << "[SCORING]\t"
            << "s1:" << s1 << " s2:" << s2 << " n1:" << n1 << " aln2mapq(n1):" << aln2mapq(snpCost, n1)
            << " a2m_scale:" << a2m_scale << " mapq:" << mapq << std::endl;
#endif  // TRACE_SCORING

  //    std::cerr << "s1:" << s1 << " s2:" << s2 << " n1:" << n1 << " aln2mapq(n1):" << aln2mapq(snpCost, n1) <<
  ////    " a2m_over:" << ((s1 - s2) * aln2mapq(snpCost, n1) > MAPQ_MAX) << -- wrong and confusing
  //      " a2m_scale:" << a2m_scale << " mapq:" << mapq <<  std::endl;

  return mapq;
}

}  // namespace align
}  // namespace dragenos

#endif  // #ifndef ALIGN_MAPQ_HPP
