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

#ifndef ALIGN_CALCULATE_REF_START_END_HPP
#define ALIGN_CALCULATE_REF_START_END_HPP

#include <utility>
#include "map/SeedChain.hpp"
#include "sequences/Read.hpp"

namespace dragenos {
namespace align {

std::pair<int64_t, int64_t> calculateRefStartEnd(const sequences::Read& read, const map::SeedChain& chain);

}  // namespace align
}  // namespace dragenos

#endif  // #ifndef ALIGN_CALCULATE_REF_START_END_HPP
