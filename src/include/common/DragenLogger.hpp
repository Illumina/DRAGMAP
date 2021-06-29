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

#ifndef COMMON_DRAGEN_LOGGER_HPP
#define COMMON_DRAGEN_LOGGER_HPP

#include <fstream>

//#define TRACE_SEED_CHAINS
//#define TRACE_ALIGNMENTS
//#define TRACE_SMITH_WATERMAN
//#define TRACE_SCORING
//#define TRACE_VECTOR_SMITH_WATERMAN

// set it to true to enable generation of mapper log files similar to those produced by the simulator
#define DEBUG_FILES false
//#define DEBUG_FILES true

namespace dragenos {
namespace common {

#define DRAGEN_LOG(os) \
  if (!DEBUG_FILES)    \
    ;                  \
  else                 \
    (os)

std::ofstream& alnResultLog();
std::ofstream& chainLog();
std::ofstream& rescueLog();
std::ofstream& smithWatermanLog();
std::ofstream& sWFetchLog();
std::ofstream& sWSteerLog();

#define DRAGEN_ALN_RESULT_LOG DRAGEN_LOG(::dragenos::common::alnResultLog())
#define DRAGEN_CHAIN_LOG DRAGEN_LOG(::dragenos::common::chainLog())
#define DRAGEN_RESCUE_LOG DRAGEN_LOG(::dragenos::common::rescueLog())
#define DRAGEN_SMITH_WATERMAN_LOG DRAGEN_LOG(::dragenos::common::smithWatermanLog())
#define DRAGEN_S_W_FETCH_LOG DRAGEN_LOG(::dragenos::common::sWFetchLog())
#define DRAGEN_S_W_STEER_LOG DRAGEN_LOG(::dragenos::common::sWSteerLog())

}  // namespace common
}  // namespace dragenos

#endif  // #ifndef COMMON_DRAGEN_LOGGER_HPP
