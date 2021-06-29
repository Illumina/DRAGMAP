// Copyright 2018 Edico Genome Corporation. All rights reserved.
// Copyright 2018-2020 Illumina, Inc. All rights reserved.
//
// This file contains confidential and proprietary information of the Edico
// Genome Corporation and is protected under the U.S. and international
// copyright and other intellectual property laws.
//
// $Id$
// $Author$
// $Change$
// $DateTime$
//

#ifndef __ELAPSED_TIMER_HPP__
#define __ELAPSED_TIMER_HPP__

#include "math.h"
#include <chrono>
#include <iomanip>
#include <sstream>
#include <string>

//-------------------------------------------------------------------------------
// Utilities for formatting time (originally in ms) for logging
//
inline std::string timeToFormattedString(uint64_t ms) {
  uint64_t print_ms = ms % 1000;
  uint64_t ss = (ms / 1000) % 60;
  uint64_t mm = (ms / (60 * 1000)) % 60;
  uint64_t hh = (ms / (60 * 60 * 1000));

  std::stringstream ss_time;

  ss_time << std::setfill('0') << std::right;
  ss_time << std::setw(2) << hh << ":";
  ss_time << std::setw(2) << mm << ":";
  ss_time << std::setw(2) << ss << ".";
  ss_time << std::setw(3) << print_ms;

  return ss_time.str();
}

inline std::string timeToSeconds(uint64_t ms) {
  float ss_seconds = round(static_cast<float>(ms) / 10) / 100;

  std::stringstream ss_time;
  ss_time << std::setprecision(2) << std::fixed << ss_seconds;

  return ss_time.str();
}

//-------------------------------------------------------------------------------
// A class that tracks elasped time.
//
class ElapsedTimer {
public:
  typedef std::chrono::high_resolution_clock clock_type;

  ElapsedTimer() : m_start(clock_::now()) {}

  void Reset() { m_start = clock_::now(); }

  double Elapsed() const {
    return std::chrono::duration_cast<second_>(clock_::now() - m_start).count();
  }

private:
  typedef std::chrono::high_resolution_clock clock_;
  typedef std::chrono::duration<double, std::ratio<1>> second_;
  std::chrono::time_point<clock_> m_start;
};

#endif
