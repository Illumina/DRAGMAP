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

#ifndef ALIGN_SCORE_HPP
#define ALIGN_SCORE_HPP

namespace dragenos {
namespace align {

typedef int            ScoreType;
static const ScoreType INVALID_SCORE = std::numeric_limits<ScoreType>::min();
static const ScoreType MAX_SCORE     = std::numeric_limits<ScoreType>::max();

}  // namespace align
}  // namespace dragenos

#endif  // #ifndef ALIGN_SCORE_HPP
