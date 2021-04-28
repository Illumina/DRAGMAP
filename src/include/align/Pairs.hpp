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

#ifndef ALIGN_PAIRS_HPP
#define ALIGN_PAIRS_HPP

#include "align/InsertSizeParameters.hpp"
#include "map/SeedChain.hpp"
#include "sequences/ReadPair.hpp"

namespace dragenos {
namespace align {

bool pairMatch(
    const InsertSizeParameters& insertSizeParameters,
    const sequences::ReadPair&  readPair,
    const map::SeedChain&       first,
    const map::SeedChain&       second);
std::pair<int, int> calculateEffBegEnd(const sequences::Read& read, const map::SeedChain& chain);

bool pairMatch(
    const InsertSizeParameters& insert_stats, const Alignment& inp_result, const Alignment& result_rrec);

}  // namespace align
}  // namespace dragenos

#endif  // #ifndef ALIGN_PAIRS_HPP
