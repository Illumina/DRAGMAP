// Copyright 2016-2018 Edico Genome Corporation. All rights reserved.
//
// This file contains confidential and proprietary information of the Edico Genome
// Corporation and is protected under the U.S. and international copyright and other
// intellectual property laws.
//
// $Id$
// $Author$
// $Change$
// $DateTime$
//

#ifndef KERNEL_DENSITY_HPP
#define KERNEL_DENSITY_HPP

#include <stdint.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <map>
#include <vector>

//--------------------------------------------------------------------------------
// KernelDensity - a class that builds up a epanechnikov kernel from a set of samples,
// to be used for insert stats estimation
//
class KernelDensity {
public:
  KernelDensity(){};

  ~KernelDensity(){};

  void AddSample(double x);

  double GetMin(int x)
  {
    m_curr_var = x;
    DefaultBandwidth();
    return m_min[x];
  };

  double GetMax(int x)
  {
    m_curr_var = x;
    DefaultBandwidth();
    return m_max[x];
  };

  double Pdf(double x);

private:
  void AddSamples(std::vector<double>& x);

  double Pdf(std::vector<double>& data);

  void CalcBandwidth();

  void DefaultBandwidth();

  double EpanechnikovPdf(const double& x, const double& mu, const double& sigma) const;

private:
  std::map<int, double>            m_sum_x;
  std::map<int, double>            m_sum_x2;
  std::map<int, double>            m_count;
  std::map<int, double>            m_min;
  std::map<int, double>            m_max;
  std::map<int, double>            m_default_bandwidth;
  std::map<int, double>            m_bandwidth;
  std::vector<std::vector<double>> m_data;
  unsigned int                     m_curr_var;
};

#endif
