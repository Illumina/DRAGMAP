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

#ifndef ALIGN_INSERT_SIZE_PARAMETERS_HPP
#define ALIGN_INSERT_SIZE_PARAMETERS_HPP

#include <math.h>
#include <algorithm>
#include <cstdint>
#include <iostream>

namespace dragenos {
namespace align {

struct InsertSizeParameters {
  enum class Orientation {
    pe_orient_fr_c = 0,
    pe_orient_rf_c,  // 01
    pe_orient_ff_c,  // 10
    pe_orient_rr_c   // 11
  };
  int min_  = 0;
  int max_  = 0;
  int mean_ = 0;

  int rescueMin_ = 0;
  int rescueMax_ = 0;
  //  double stddev_;
  uint16_t    sigmaFactor_ = 0;
  Orientation orientation_ = Orientation::pe_orient_fr_c;

  bool     isInitStatDone_ = false;
  uint16_t getSigmaFactor() const
  {
    return sigmaFactor_;
    //    return std::min(static_cast<uint16_t>(0xFFFF), static_cast<uint16_t>(std::round(double(0x2F200) / stddev_)));
  }

  bool isInitDone() const { return isInitStatDone_; }

  InsertSizeParameters() {}

  InsertSizeParameters(
      int         min,
      int         max,
      int         mean,
      int         rescueMin,
      int         rescueMax,
      uint16_t    sigmaFactor,
      Orientation orientation,
      bool        isInitStatDone)
    : min_(min),
      max_(max),
      mean_(mean),
      rescueMin_(rescueMin),
      rescueMax_(rescueMax),
      sigmaFactor_(sigmaFactor),
      orientation_(orientation),
      isInitStatDone_(isInitStatDone)
  {
  }

  InsertSizeParameters(
      int min, int max, int mean, int rescueMin, int rescueMax, double stddev, Orientation orientation)
    : min_(min),
      max_(max),
      mean_(mean),
      rescueMin_(rescueMin),
      rescueMax_(rescueMax),
      sigmaFactor_(std::min(
          static_cast<uint16_t>(0xFFFF), static_cast<uint16_t>(std::round(double(0x2F200) / stddev)))),
      orientation_(orientation)
  {
  }

  friend std::ostream& operator<<(std::ostream& os, const InsertSizeParameters& isp)
  {
    return os << "InsertSizeParameters(" << isp.min_ << "," << isp.max_ << "," << isp.mean_ << ","
              << isp.rescueMin_ << "," << isp.rescueMax_ << "," << isp.sigmaFactor_ << ","
              << int(isp.orientation_) << "," << isp.isInitStatDone_ << ")";
  }
};

}  // namespace align
}  // namespace dragenos

#endif  // #ifndef ALIGN_INSERT_SIZE_PARAMETERS_HPP
