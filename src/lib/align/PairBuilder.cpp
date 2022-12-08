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

#include "align/PairBuilder.hpp"
#include "align/AlignmentRescue.hpp"
#include "align/SinglePicker.hpp"
#include "align/Tlen.hpp"
#include "common/DragenLogger.hpp"

namespace dragenos {
namespace align {

//-- This ROM contains the phred-scale (-10log10) probability from a normal PDF, with the output scaled to 1
// at input 0,
//-- and the input scaled to saturate 8 bits (0xFF) at the end of a 9-bit range (511).
//--    N = 0..511:  rom(N) = round(-10*log10(standardNormalPDF(N/47.125)/standardNormalPDF(0)))
//-- To use it, given an expected insert-size standard deviation "sigma" we have config register
//--    sigma_factor = min(0xFFFF, round(0x2F200/sigma)), interpreted as 4.12 bit fixed point.
//--    (This saturates when sigma is sligtly less than 3.)
//-- Then look up:
//--    N = int(sigma_factor * abs(insert_size - mean))
//--  constant petab_rom_c : petab_rom_t := (
//--    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00",
//--    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01",
//--    0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x02, 0x02, 0x02, 0x02, 0x02, 0x02, 0x02, 0x02",
//--    0x02, 0x02, 0x02, 0x03, 0x03, 0x03, 0x03, 0x03, 0x03, 0x03, 0x03, 0x03, 0x04, 0x04, 0x04, 0x04",
//--    0x04, 0x04, 0x04, 0x04, 0x05, 0x05, 0x05, 0x05, 0x05, 0x05, 0x05, 0x06, 0x06, 0x06, 0x06, 0x06",
//--    0x06, 0x06, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x09, 0x09",
//--    0x09, 0x09, 0x09, 0x0A, 0x0A, 0x0A, 0x0A, 0x0A, 0x0B, 0x0B, 0x0B, 0x0B, 0x0B, 0x0C, 0x0C, 0x0C",
//--    0x0C, 0x0C, 0x0D, 0x0D, 0x0D, 0x0D, 0x0E, 0x0E, 0x0E, 0x0E, 0x0F, 0x0F, 0x0F, 0x0F, 0x10, 0x10",
//--    0x10, 0x10, 0x11, 0x11, 0x11, 0x11, 0x12, 0x12, 0x12, 0x12, 0x13, 0x13, 0x13, 0x13, 0x14, 0x14",
//--    0x14, 0x15, 0x15, 0x15, 0x15, 0x16, 0x16, 0x16, 0x17, 0x17, 0x17, 0x17, 0x18, 0x18, 0x18, 0x19",
//--    0x19, 0x19, 0x1A, 0x1A, 0x1A, 0x1B, 0x1B, 0x1B, 0x1C, 0x1C, 0x1C, 0x1D, 0x1D, 0x1D, 0x1E, 0x1E",
//--    0x1E, 0x1F, 0x1F, 0x1F, 0x20, 0x20, 0x20, 0x21, 0x21, 0x21, 0x22, 0x22, 0x23, 0x23, 0x23, 0x24",
//--    0x24, 0x24, 0x25, 0x25, 0x26, 0x26, 0x26, 0x27, 0x27, 0x28, 0x28, 0x28, 0x29, 0x29, 0x29, 0x2A",
//--    0x2A, 0x2B, 0x2B, 0x2C, 0x2C, 0x2C, 0x2D, 0x2D, 0x2E, 0x2E, 0x2E, 0x2F, 0x2F, 0x30, 0x30, 0x31",
//--    0x31, 0x32, 0x32, 0x32, 0x33, 0x33, 0x34, 0x34, 0x35, 0x35, 0x36, 0x36, 0x36, 0x37, 0x37, 0x38",
//--    0x38, 0x39, 0x39, 0x3A, 0x3A, 0x3B, 0x3B, 0x3C, 0x3C, 0x3D, 0x3D, 0x3E, 0x3E, 0x3F, 0x3F, 0x40",
//--    0x40, 0x41, 0x41, 0x42, 0x42, 0x43, 0x43, 0x44, 0x44, 0x45, 0x45, 0x46, 0x46, 0x47, 0x47, 0x48",
//--    0x48, 0x49, 0x49, 0x4A, 0x4A, 0x4B, 0x4C, 0x4C, 0x4D, 0x4D, 0x4E, 0x4E, 0x4F, 0x4F, 0x50, 0x51",
//--    0x51, 0x52, 0x52, 0x53, 0x53, 0x54, 0x55, 0x55, 0x56, 0x56, 0x57, 0x57, 0x58, 0x59, 0x59, 0x5A",
//--    0x5A, 0x5B, 0x5C, 0x5C, 0x5D, 0x5D, 0x5E, 0x5F, 0x5F, 0x60, 0x60, 0x61, 0x62, 0x62, 0x63, 0x64",
//--    0x64, 0x65, 0x65, 0x66, 0x67, 0x67, 0x68, 0x69, 0x69, 0x6A, 0x6A, 0x6B, 0x6C, 0x6C, 0x6D, 0x6E",
//--    0x6E, 0x6F, 0x70, 0x70, 0x71, 0x72, 0x72, 0x73, 0x74, 0x74, 0x75, 0x76, 0x76, 0x77, 0x78, 0x78",
//--    0x79, 0x7A, 0x7B, 0x7B, 0x7C, 0x7D, 0x7D, 0x7E, 0x7F, 0x7F, 0x80, 0x81, 0x82, 0x82, 0x83, 0x84",
//--    0x84, 0x85, 0x86, 0x87, 0x87, 0x88, 0x89, 0x8A, 0x8A, 0x8B, 0x8C, 0x8C, 0x8D, 0x8E, 0x8F, 0x8F",
//--    0x90, 0x91, 0x92, 0x92, 0x93, 0x94, 0x95, 0x95, 0x96, 0x97, 0x98, 0x99, 0x99, 0x9A, 0x9B, 0x9C",
//--    0x9C, 0x9D, 0x9E, 0x9F, 0xA0, 0xA0, 0xA1, 0xA2, 0xA3, 0xA4, 0xA4, 0xA5, 0xA6, 0xA7, 0xA8, 0xA8",
//--    0xA9, 0xAA, 0xAB, 0xAC, 0xAC, 0xAD, 0xAE, 0xAF, 0xB0, 0xB1, 0xB1, 0xB2, 0xB3, 0xB4, 0xB5, 0xB6",
//--    0xB6, 0xB7, 0xB8, 0xB9, 0xBA, 0xBB, 0xBC, 0xBC, 0xBD, 0xBE, 0xBF, 0xC0, 0xC1, 0xC2, 0xC3, 0xC3",
//--    0xC4, 0xC5, 0xC6, 0xC7, 0xC8, 0xC9, 0xCA, 0xCA, 0xCB, 0xCC, 0xCD, 0xCE, 0xCF, 0xD0, 0xD1, 0xD2",
//--    0xD3, 0xD3, 0xD4, 0xD5, 0xD6, 0xD7, 0xD8, 0xD9, 0xDA, 0xDB, 0xDC, 0xDD, 0xDE, 0xDE, 0xDF, 0xE0",
//--    0xE1, 0xE2, 0xE3, 0xE4, 0xE5, 0xE6, 0xE7, 0xE8, 0xE9, 0xEA, 0xEB, 0xEC, 0xED, 0xEE, 0xEF, 0xF0",
//--    0xF1, 0xF2, 0xF2, 0xF3, 0xF4, 0xF5, 0xF6, 0xF7, 0xF8, 0xF9, 0xFA, 0xFB, 0xFC, 0xFD, 0xFE, 0xFF");
//-- Alternate version matching legacy DRAGEN penalties defined by Gaussian 2-tail area ("P-value"):
//--    N = 0..511:  rom(N) = round(-10*log10(2*standardNormalCDF(-N/47.125)))
//-- These penalties are higher than the above, e.g.:  0 -> 0, 2 -> 5, 6 -> 10, 15 -> 20, 33 -> 40, 71 -> 80.
//-- They saturate at 0xFF a little sooner.  Although the PDF version above is theoretically more correct
// assuming
//-- the normal distribution is true, we saw some accuracy degradation, so we are rolling back to the legacy
// definition.
static const int petab_rom_c[] = {
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01,
    0x01, 0x02, 0x02, 0x02, 0x02, 0x02, 0x02, 0x02, 0x02, 0x02, 0x02, 0x03, 0x03, 0x03, 0x03, 0x03, 0x03,
    0x03, 0x03, 0x04, 0x04, 0x04, 0x04, 0x04, 0x04, 0x04, 0x04, 0x05, 0x05, 0x05, 0x05, 0x05, 0x05, 0x05,
    0x06, 0x06, 0x06, 0x06, 0x06, 0x06, 0x06, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x08, 0x08, 0x08, 0x08,
    0x08, 0x08, 0x09, 0x09, 0x09, 0x09, 0x09, 0x0A, 0x0A, 0x0A, 0x0A, 0x0A, 0x0A, 0x0B, 0x0B, 0x0B, 0x0B,
    0x0B, 0x0C, 0x0C, 0x0C, 0x0C, 0x0D, 0x0D, 0x0D, 0x0D, 0x0D, 0x0E, 0x0E, 0x0E, 0x0E, 0x0E, 0x0F, 0x0F,
    0x0F, 0x0F, 0x10, 0x10, 0x10, 0x10, 0x11, 0x11, 0x11, 0x11, 0x12, 0x12, 0x12, 0x12, 0x13, 0x13, 0x13,
    0x13, 0x14, 0x14, 0x14, 0x14, 0x15, 0x15, 0x15, 0x16, 0x16, 0x16, 0x16, 0x17, 0x17, 0x17, 0x18, 0x18,
    0x18, 0x18, 0x19, 0x19, 0x19, 0x1A, 0x1A, 0x1A, 0x1A, 0x1B, 0x1B, 0x1B, 0x1C, 0x1C, 0x1C, 0x1D, 0x1D,
    0x1D, 0x1E, 0x1E, 0x1E, 0x1F, 0x1F, 0x1F, 0x20, 0x20, 0x20, 0x21, 0x21, 0x21, 0x22, 0x22, 0x22, 0x23,
    0x23, 0x23, 0x24, 0x24, 0x25, 0x25, 0x25, 0x26, 0x26, 0x26, 0x27, 0x27, 0x27, 0x28, 0x28, 0x29, 0x29,
    0x29, 0x2A, 0x2A, 0x2B, 0x2B, 0x2B, 0x2C, 0x2C, 0x2D, 0x2D, 0x2D, 0x2E, 0x2E, 0x2F, 0x2F, 0x2F, 0x30,
    0x30, 0x31, 0x31, 0x32, 0x32, 0x32, 0x33, 0x33, 0x34, 0x34, 0x35, 0x35, 0x35, 0x36, 0x36, 0x37, 0x37,
    0x38, 0x38, 0x39, 0x39, 0x39, 0x3A, 0x3A, 0x3B, 0x3B, 0x3C, 0x3C, 0x3D, 0x3D, 0x3E, 0x3E, 0x3F, 0x3F,
    0x40, 0x40, 0x41, 0x41, 0x42, 0x42, 0x42, 0x43, 0x43, 0x44, 0x44, 0x45, 0x45, 0x46, 0x46, 0x47, 0x48,
    0x48, 0x49, 0x49, 0x4A, 0x4A, 0x4B, 0x4B, 0x4C, 0x4C, 0x4D, 0x4D, 0x4E, 0x4E, 0x4F, 0x4F, 0x50, 0x51,
    0x51, 0x52, 0x52, 0x53, 0x53, 0x54, 0x54, 0x55, 0x55, 0x56, 0x57, 0x57, 0x58, 0x58, 0x59, 0x59, 0x5A,
    0x5B, 0x5B, 0x5C, 0x5C, 0x5D, 0x5E, 0x5E, 0x5F, 0x5F, 0x60, 0x61, 0x61, 0x62, 0x62, 0x63, 0x64, 0x64,
    0x65, 0x65, 0x66, 0x67, 0x67, 0x68, 0x68, 0x69, 0x6A, 0x6A, 0x6B, 0x6C, 0x6C, 0x6D, 0x6E, 0x6E, 0x6F,
    0x6F, 0x70, 0x71, 0x71, 0x72, 0x73, 0x73, 0x74, 0x75, 0x75, 0x76, 0x77, 0x77, 0x78, 0x79, 0x79, 0x7A,
    0x7B, 0x7B, 0x7C, 0x7D, 0x7D, 0x7E, 0x7F, 0x7F, 0x80, 0x81, 0x82, 0x82, 0x83, 0x84, 0x84, 0x85, 0x86,
    0x86, 0x87, 0x88, 0x89, 0x89, 0x8A, 0x8B, 0x8B, 0x8C, 0x8D, 0x8E, 0x8E, 0x8F, 0x90, 0x91, 0x91, 0x92,
    0x93, 0x94, 0x94, 0x95, 0x96, 0x97, 0x97, 0x98, 0x99, 0x9A, 0x9A, 0x9B, 0x9C, 0x9D, 0x9D, 0x9E, 0x9F,
    0xA0, 0xA0, 0xA1, 0xA2, 0xA3, 0xA4, 0xA4, 0xA5, 0xA6, 0xA7, 0xA8, 0xA8, 0xA9, 0xAA, 0xAB, 0xAC, 0xAC,
    0xAD, 0xAE, 0xAF, 0xB0, 0xB0, 0xB1, 0xB2, 0xB3, 0xB4, 0xB5, 0xB5, 0xB6, 0xB7, 0xB8, 0xB9, 0xBA, 0xBA,
    0xBB, 0xBC, 0xBD, 0xBE, 0xBF, 0xBF, 0xC0, 0xC1, 0xC2, 0xC3, 0xC4, 0xC5, 0xC5, 0xC6, 0xC7, 0xC8, 0xC9,
    0xCA, 0xCB, 0xCC, 0xCC, 0xCD, 0xCE, 0xCF, 0xD0, 0xD1, 0xD2, 0xD3, 0xD4, 0xD4, 0xD5, 0xD6, 0xD7, 0xD8,
    0xD9, 0xDA, 0xDB, 0xDC, 0xDD, 0xDD, 0xDE, 0xDF, 0xE0, 0xE1, 0xE2, 0xE3, 0xE4, 0xE5, 0xE6, 0xE7, 0xE8,
    0xE9, 0xEA, 0xEA, 0xEB, 0xEC, 0xED, 0xEE, 0xEF, 0xF0, 0xF1, 0xF2, 0xF3, 0xF4, 0xF5, 0xF6, 0xF7, 0xF8,
    0xF9, 0xFA, 0xFB, 0xFC, 0xFD, 0xFE, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0xFF};

int PairBuilder::computePairPenalty(
    const InsertSizeParameters& insertSizeParameters,
    const ReadPair&             readPair,
    const Alignment*            a2,
    const Alignment*            a1,
    const bool                  properPair,
    int&                        insert_len,
    int&                        insert_diff) const
{
  int m2a_penalty = -1;
  if (!a1 || !a2 || !properPair) {
    m2a_penalty = aln_cfg_unpaired_pen_;
  } else {
    //      std::cerr << "c1:" << *c1 << " c2:" << *c2 <<std::endl;
    //    const int inp_result_qry_end_gap =
    //        a1->getCigar()
    //            .countEndClips();  //readPair[0].getLength() - c1->lastReadBase() - 1;  //back().getSeed().getReadPosition() - 1;
    //      std::cerr << "inp_result_qry_end_gap:" << inp_result_qry_end_gap << std::endl;
    //    const int result_rrec_qry_end_gap =
    //        a2->getCigar()
    //            .countEndClips();  //readPair[1].getLength() - c2->lastReadBase() - 1;  //back().getSeed().getReadPosition() - 1;
    //      std::cerr << "result_rrec_qry_end_gap:" << result_rrec_qry_end_gap << std::endl;

    const int inp_eff_beg = a1->getUnclippedStartPosition();
    const int inp_eff_end = a1->getUnclippedEndPosition();
    const int rct_eff_beg = a2->getUnclippedStartPosition();
    const int rct_eff_end = a2->getUnclippedEndPosition();

    //      std::cerr << " inp_eff_beg=" << inp_eff_beg << " inp_eff_end=" << inp_eff_end << " rct_eff_beg=" << rct_eff_beg << " rct_eff_end=" << rct_eff_end<< std::endl;

    using Orientation = InsertSizeParameters::Orientation;

    uint32_t insert_beg = 0;
    uint32_t insert_end = 0;
    // -- Forward-forward orientation (same as reverse-reverse): pick the outermost begin and end
    if (insertSizeParameters.orientation_ == Orientation::pe_orient_ff_c or
        insertSizeParameters.orientation_ == Orientation::pe_orient_rr_c) {
      //    std::cerr << "ff or rr" << std::endl;
      insert_beg = std::min(inp_eff_beg, rct_eff_beg);
      insert_end = std::max(inp_eff_end, rct_eff_end);
      //    wrong_beg_v := '0' & inp_result.ref_beg;
      //    wrong_end_v := '0' & inp_result.ref_end;
    }
    //-- Forward-reverse: beginning from the forward one, end from the reverse one
    //-- Reverse-forward: beginning from the reverse one, end from the forward one
    else if ((insertSizeParameters.orientation_ == Orientation::pe_orient_fr_c) ^ a1->isReverseComplement()) {
      //    std::cerr << "fr" << std::endl;
      insert_beg = inp_eff_beg;
      insert_end = rct_eff_end;
      //    wrong_beg_v := '0' & result_rrec.ref_beg;
      //    wrong_end_v := '0' & inp_result.ref_end;
    } else {
      //    std::cerr << "rf?" << std::endl;
      insert_beg = rct_eff_beg;
      insert_end = inp_eff_end;
      //    wrong_beg_v := '0' & inp_result.ref_beg;
      //    wrong_end_v := '0' & result_rrec.ref_end;
    }

    insert_len = insert_end - insert_beg + 1;
    //        std::cerr << "insert_beg:" << insert_beg << " insert_end:" << insert_end << " insert_len:" << insert_len << std::endl;

    insert_diff = std::abs(insert_len - insertSizeParameters.mean_);
    //      std::cerr << "sigma factor:" << insertSizeParameters.getSigmaFactor() << " pe_mean_insert:" << insertSizeParameters.mean_ << " insert_diff:" << insert_diff << std::endl;
    const unsigned   ins_prod                 = insert_diff * insertSizeParameters.getSigmaFactor();
    static const int sigma_factor_frac_bits_c = 12;
    static const int petab_addr_bits_c        = 9;
    const unsigned   ins_prod_q = ((1 << petab_addr_bits_c) - 1) & (ins_prod >> sigma_factor_frac_bits_c);
    //    std::cerr << "ins_prod:" << ins_prod << " ins_prod_q:" << ins_prod_q << std::endl;

    m2a_penalty = ins_prod_q < sizeof(petab_rom_c) / sizeof(*petab_rom_c) ? petab_rom_c[ins_prod_q]
                                                                          : aln_cfg_unpaired_pen_;
  }

  const int m2a_scale =
      mapq2aln(similarity_.getSnpCost(), std::max(aln_cfg_mapq_min_len_, readPair.getLength()));
  //    std::cerr << "m2a_scale:" << m2a_scale << std::endl;
  const int m2a_prod = (m2a_scale * m2a_penalty) >> 10;

  //  const int m2a_prod_frac_bits_c = 10;
  const int m2a_prod_pen = m2a_prod;  // >> m2a_prod_frac_bits_c;

  //      std::cerr << "m2a_prod_pen:" << m2a_prod_pen << std::endl;
  return m2a_prod_pen;
}

/**
 * \brief Find second best pair such that the read at readIdx is not part of it
 * \param nearSub - On return, count of alignments that are withing 1 SNP of suboptimal.
 *        Unchanged if no suboptimal found
 */
AlignmentPairs::const_iterator PairBuilder::findSecondBest(
    const int                            averageReadLength,
    const AlignmentPairs&                pairs,
    const UnpairedAlignments&            unpairedAlignments,
    const AlignmentPairs::const_iterator best,
    int                                  readIdx,
    int&                                 sub_count) const
{
  AlignmentPairs::const_iterator ret = pairs.end();
  for (AlignmentPairs::const_iterator it = pairs.begin(); pairs.end() != it; ++it) {
    //    std::cerr << "trying:" << *it << std::endl;
    if ((pairs.end() == ret || ret->getScore() < it->getScore()) &&
        !best->at(readIdx).isDuplicate(it->at(readIdx))) {
      ret = it;
    }
  }

  if (pairs.end() != ret) {
#ifdef TRACE_SCORING
    std::cerr << "[SCORING]\t"
              << "best:" << *best << std::endl;
    std::cerr << "[SCORING]\t"
              << "second best:" << *ret << std::endl;
#endif  // TRACE_SCORING
    ScoreType other_best_scr = best->at(!readIdx).getScore();
    const int m2a_scale =
        mapq2aln(similarity_.getSnpCost(), std::max(aln_cfg_mapq_min_len_, averageReadLength));
    const ScoreType scaled_max_pen = (m2a_scale * aln_cfg_unpaired_pen_) >> 10;  //27;

    const bool      use_pe_mapq_v = ret->isProperPair();
    const ScoreType list_pe_max =
        use_pe_mapq_v ? ret->getScore() : (ret->at(readIdx).getScore() + other_best_scr - scaled_max_pen);
    const ScoreType list_pe_min = list_pe_max - similarity_.getSnpCost();
    const ScoreType list_se_max =
        use_pe_mapq_v ? (ret->getScore() - other_best_scr + scaled_max_pen) : ret->at(readIdx).getScore();
    const ScoreType list_se_min = list_se_max - similarity_.getSnpCost();
#ifdef TRACE_SCORING
    std::cerr << "[SCORING]\t"
              << "list_pe_min:" << list_pe_min << " list_pe_max:" << list_pe_max << std::endl;
    std::cerr << "[SCORING]\t"
              << "list_se_min:" << list_se_min << " list_se_max:" << list_se_max
              << " sub_pair_score_v=" << ret->getScore()
              << " other_mins_pen_v=" << (other_best_scr - scaled_max_pen)
              << " se_best_score=" << other_best_scr << std::endl;
#endif  // TRACE_SCORING

    if (ret->isProperPair()) {
      for (AlignmentPairs::const_iterator it = pairs.begin(); pairs.end() != it; ++it) {
        if (it->isProperPair() && it->getScore() > list_pe_min && it->getScore() <= list_pe_max) {
          ++sub_count;
#ifdef TRACE_SCORING
          std::cerr << "[SCORING]\t"
                    << "psub_count:" << sub_count << " " << *it << std::endl;
#endif  // TRACE_SCORING
        }
      }
    }

    // std::cerr << "second best: " << *ret << std::endl;
    // std::cerr << "other_best_scr: " << other_best_scr << std::endl;

    for (Alignments::const_iterator it = unpairedAlignments[readIdx].begin();
         unpairedAlignments[readIdx].end() != it;
         ++it) {
      if (!it->isUnmapped() && it->getScore() > list_se_min && it->getScore() <= list_se_max) {
        ++sub_count;
#ifdef TRACE_SCORING
        std::cerr << "[SCORING]\t"
                  << "ssub_count:" << sub_count << " " << *it << std::endl;
#endif  // TRACE_SCORING
      }
    }
  }

  return ret;
}

AlignmentPairs::const_iterator PairBuilder::findSecondBestScore(
    const int                            averageReadLength,
    const AlignmentPairs&                pairs,
    const UnpairedAlignments&            unpairedAlignments,
    const AlignmentPairs::const_iterator best,
    int                                  readIdx,
    int&                                 sub_count,
    ScoreType&                           secondBestPeScore,
    ScoreType&                           secondBestSeScore) const
{
  AlignmentPairs::const_iterator secondBest =
      findSecondBest(averageReadLength, pairs, unpairedAlignments, best, readIdx, sub_count);

  int se_sub_count  = 0;
  secondBestSeScore = std::max(
      alnMinScore_,
      align::findSecondBestScore(
          similarity_.getSnpCost(), unpairedAlignments[readIdx], &best->at(readIdx), se_sub_count));

  if (pairs.end() != secondBest) {
    secondBestPeScore = secondBest->getScore();
#ifdef TRACE_SCORING
    std::cerr << "[SCORING]\t"
              << "findSecondBestScore\tr" << readIdx << " sub_count:" << sub_count
              << " se_sub_count:" << se_sub_count << " secondBestSeScore:" << secondBestSeScore
              << " secondBestPeScore:" << secondBestPeScore << std::endl;
#endif  // TRACE_SCORING
  } else {
#ifdef TRACE_SCORING
    std::cerr << "[SCORING]\t"
              << "findSecondBestScore\tr" << readIdx << " sub_count:" << sub_count
              << " se_sub_count:" << se_sub_count << " secondBestSeScore:" << secondBestSeScore << std::endl;
#endif  // TRACE_SCORING
    sub_count = se_sub_count;
  }
  return secondBest;
}

void PairBuilder::updateEndMapq(
    const int                averageReadLength,
    AlignmentPairs::iterator best,
    const int                readIdx,
    int                      sub_count,
    ScoreType                sub_pair_score_v,
    const ScoreType          secondBestSeScore[],
    const ScoreType          xs[]) const
{
  if (best->at(readIdx).isUnmapped()) {
    best->at(readIdx).setMapq(0);
  } else {
    const int a2m_scale =
        aln2mapq(similarity_.getSnpCost(), std::max(aln_cfg_mapq_min_len_, averageReadLength));
    const ScoreType xs_score_diff = best->getScore() - xs[0] - xs[1];
    const MapqType xs_heur_mapq = std::max(0, ((xs_score_diff * a2m_scale) >> 13) + aln_cfg_xs_pair_penalty_);

    const MapqType a2m_mapq = INVALID_SCORE == sub_pair_score_v
                                  ? computeMapq(
                                        similarity_.getSnpCost(),
                                        best->at(readIdx).getScore(),
                                        std::max(alnMinScore_, secondBestSeScore[readIdx]),
                                        std::max(aln_cfg_mapq_min_len_, averageReadLength))
                                  : computeMapq(
                                        similarity_.getSnpCost(),
                                        best->getScore(),
                                        std::max(alnMinScore_, sub_pair_score_v),
                                        std::max(aln_cfg_mapq_min_len_, averageReadLength));

    const MapqType pe_mapq = a2m_mapq;

    const int      sub_count_log2 = log2_simple(sub_count);
    const MapqType sub_mapq_pen_v = sub_count ? (3 * sub_count_log2) : 0;

    const MapqType mapq_prod_pen = pe_mapq - (sub_mapq_pen_v >> 7);

    const bool mapq0 = (aln_cfg_sample_mapq0_ >= 1 && best->hasOnlyRandomSamples()) ||
                       (aln_cfg_sample_mapq0_ >= 2 && best->isExtra());
    const MapqType mapq = mapq0 ? 0
                                : (INVALID_SCORE != xs_score_diff)
                                      ? std::min(std::max(0, xs_heur_mapq), mapq_prod_pen)
                                      : mapq_prod_pen;

#ifdef TRACE_SCORING
    std::cerr << "[SCORING]\t"
              << "r" << readIdx << " a2m_scale:" << a2m_scale
              << " se_score_diff:" << (best->at(readIdx).getScore() - secondBestSeScore[readIdx])
              << " best_pair_score:" << best->getScore() << " sub_pair_score_v:" << sub_pair_score_v
              << " pe_score_diff:" << (best->getScore() - sub_pair_score_v)
              << " xs_score_diff:" << xs_score_diff << " result_rsub:" << xs[readIdx]
              << " mate_sub_scr:" << xs[!readIdx] << " sub_count:" << sub_count
              << " sub_mapq_pen_v:" << sub_mapq_pen_v << " xs_heur_mapq:" << xs_heur_mapq
              << " pe_mapq:" << pe_mapq << " mapq:" << mapq << std::endl;
#endif  // TRACE_SCORING
#ifdef TRACE_SCORING
    std::cerr << "[SCORING]\t"
              << "a2m_mapq=" << a2m_mapq << " sub_count=" << sub_count << " sub_count_log2=" << sub_count_log2
              << " sub_mapq_pen_v=" << sub_mapq_pen_v << " mapq_prod_pen=" << mapq_prod_pen << std::endl;
#endif  // TRACE_SCORING
    best->at(readIdx).setMapq(mapq);

    best->at(readIdx).setXs(xs[readIdx] > alnMinScore_ ? xs[readIdx] : INVALID_SCORE);
  }
}

void PairBuilder::updateMapq(
    const int                 averageReadLength,
    AlignmentPairs&           pairs,
    const UnpairedAlignments& unpairedAlignments,
    AlignmentPairs::iterator  best) const
{
  int       sub_count[2]         = {0, 0};
  ScoreType sub_pair_score_v[2]  = {INVALID_SCORE, INVALID_SCORE};
  ScoreType secondBestSeScore[2] = {INVALID_SCORE, INVALID_SCORE};

  //  AlignmentPairs::const_iterator secondBestPair[2] = {findSecondBestScore(
  //                                                          averageReadLength,
  //                                                          pairs,
  //                                                          unpairedAlignments,
  //                                                          best,
  //                                                          0,
  //                                                          sub_count[0],
  //                                                          sub_pair_score_v[0],
  //                                                          secondBestSeScore[0]),
  //
  //                                                      findSecondBestScore(
  //                                                          averageReadLength,
  //                                                          pairs,
  //                                                          unpairedAlignments,
  //                                                          best,
  //                                                          1,
  //                                                          sub_count[1],
  //                                                          sub_pair_score_v[1],
  //                                                          secondBestSeScore[1])};

  findSecondBestScore(
      averageReadLength,
      pairs,
      unpairedAlignments,
      best,
      0,
      sub_count[0],
      sub_pair_score_v[0],
      secondBestSeScore[0]);

  findSecondBestScore(
      averageReadLength,
      pairs,
      unpairedAlignments,
      best,
      1,
      sub_count[1],
      sub_pair_score_v[1],
      secondBestSeScore[1]);

  updateEndMapq(
      averageReadLength, best, 0, sub_count[0], sub_pair_score_v[0], secondBestSeScore, secondBestSeScore);

  updateEndMapq(
      averageReadLength, best, 1, sub_count[1], sub_pair_score_v[1], secondBestSeScore, secondBestSeScore);
}

AlignmentPairs::iterator PairBuilder::pickBest(
    const ReadPair&           readPair,
    AlignmentPairs&           alignmentPairs,
    const UnpairedAlignments& unpairedAlignments) const
{
  if (alignmentPairs.empty()) {
    return alignmentPairs.end();
  }

  AlignmentPairs::iterator best = std::max_element(
      alignmentPairs.begin(),
      alignmentPairs.end(),
      [](const AlignmentPair& left, const AlignmentPair& right) {
        return left.getScore() < right.getScore();
      });

  if (true == best->at(0).getIneligibilityStatus() and true == best->at(1).getIneligibilityStatus())
    return alignmentPairs.end();

  updateMapq(readPair.getLength(), alignmentPairs, unpairedAlignments, best);
  //  std::cerr << "best: " << *best << std::endl;

  return best;
}

}  // namespace align
}  // namespace dragenos
