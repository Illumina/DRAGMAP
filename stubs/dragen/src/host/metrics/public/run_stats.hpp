// Copyright 2014 Edico Genome Corporation. All rights reserved.
// This file contains confidential and proprietary information of the Edico
// Genome Corporation and is protected under the U.S. and international
// copyright and other intellectual property laws.
//
//

#pragma once

#include <string>
#include <vector>
#include <memory>
//
// RP: HA! HA! HA! That's what you get when you write code logging to cout all
// over the place!
#define cout cerr

class RunStats {
public:
  struct RGinsertSizeStats {
    RGinsertSizeStats()
        : m_errMsg(""), m_numInsertsInMean(0), m_q25(0), m_q50(0), m_q75(0),
          m_mean(0), m_stddev(0), m_low(0), m_high(0), m_minInsert(0),
          m_maxInsert(0) {}
    std::string m_errMsg;
    int m_numInsertsInMean;
    int m_q25;
    int m_q50;
    int m_q75;
    double m_mean;
    double m_stddev;
    int m_low;
    int m_high;
    uint32_t m_minInsert;
    uint32_t m_maxInsert;
  };

  void addInsertSizeRG(RGinsertSizeStats &other) {}
  void addInsertSizeHistogram(const std::vector<uint64_t> &){}

  static RunStats *Instance();
  static RunStats *m_instance;
};
