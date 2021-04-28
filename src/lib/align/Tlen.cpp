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

#include "align/Tlen.hpp"

namespace dragenos {
namespace align {

bool isDefaultLeft(
    const dragenos::align::Alignment& r1Align,
    const dragenos::align::Alignment& r2Align,
    InsertSizeParameters::Orientation prevailingOrientation)
{
  if (r1Align.isReverseComplement() != r2Align.isReverseComplement() &&
      prevailingOrientation != InsertSizeParameters::Orientation::pe_orient_ff_c) {
    return (!r1Align.isReverseComplement()) ^
           (prevailingOrientation == InsertSizeParameters::Orientation::pe_orient_rf_c);
  } else {
    return true;
  }
}

int getTlen(
    const dragenos::align::Alignment& r1Align,
    const dragenos::align::Alignment& r2Align,
    InsertSizeParameters::Orientation prevailingOrientation)
{
  if (r1Align.isUnmapped() || r2Align.isUnmapped() || r1Align.getReference() != r2Align.getReference()) {
    return 0;
  }

  const int tlen_beg = std::min(r1Align.getPosition(), r2Align.getPosition());
  const int tlen_end = std::max(
      r1Align.getPosition() + r1Align.getRefLen() - 1, r2Align.getPosition() + r2Align.getRefLen() - 1);

  const bool default_left = isDefaultLeft(r1Align, r2Align, prevailingOrientation);

  const int midpoint0 = r1Align.getPosition() * 2 + r1Align.getRefLen();
  const int midpoint1 = r2Align.getPosition() * 2 + r2Align.getRefLen();

  const int  abs_tlen_m1 = std::abs(tlen_end - tlen_beg);
  const bool leftmost    = midpoint0 < midpoint1 || (midpoint0 == midpoint1 && default_left);

  //    std::cerr << "cigar0:" << r1Align.getCigar() << " r1Align.getTemplateLength():" << r1Align.getTemplateLength() << std::endl;
  //    std::cerr << "cigar1:" << r2Align.getCigar() << " r2Align.getTemplateLength():" << r2Align.getTemplateLength() << std::endl;
  //    std::cerr << "midpoint0:" << midpoint0 << " midpoint1:" << midpoint1 << " leftmost:" << leftmost << " default_left:"  << default_left << std::endl;

  return leftmost ? (1 + abs_tlen_m1) : (-1 - abs_tlen_m1);
}

}  // namespace align
}  // namespace dragenos
