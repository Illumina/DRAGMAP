// Copyright 2014 Edico Genome Corporation. All rights reserved.
// This file contains confidential and proprietary information of the Edico Genome
// Corporation and is protected under the U.S. and international copyright and other
// intellectual property laws.
//
//

#include "run_stats.hpp"

RunStats* RunStats::m_instance;
RunStats* RunStats::Instance()
{
  if (!m_instance) {
    m_instance = new RunStats;
  }
  return m_instance;
}
