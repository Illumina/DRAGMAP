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

#include "common/DragenLogger.hpp"

namespace dragenos {
namespace common {

std::ofstream& chainLog()
{
  static std::ofstream log("map_chain_mgr.log");
  return log;
}

std::ofstream& rescueLog()
{
  static std::ofstream log("map_resc_scan.log");
  return log;
}

std::ofstream& smithWatermanLog()
{
  static std::ofstream log("map_sw_scores_0_tog.log");
  return log;
}

std::ofstream& sWSteerLog()
{
  static std::ofstream log("map_sw_steer_0.log");
  return log;
}

}  // namespace common
}  // namespace dragenos
