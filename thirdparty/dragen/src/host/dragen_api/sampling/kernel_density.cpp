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

#include <assert.h>

#include <stdint.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include "kernel_density.hpp"
#include <map>
#include <vector>

//#include "dragen_exception.hpp"

void KernelDensity::CalcBandwidth()
{
  for (m_curr_var = 0; m_curr_var < m_data.size(); m_curr_var++) {
    if (m_bandwidth[m_curr_var] == -1.0) {
      DefaultBandwidth();
    }
  }
}

double KernelDensity::Pdf(double x)
{
  std::vector<double> tmp;
  tmp.push_back(x);
  return (Pdf(tmp));
}

double KernelDensity::Pdf(std::vector<double>& data)
{
  CalcBandwidth();
  double d = 0.0;
  for (unsigned int i = 0; i < m_data[0].size(); i++) {
    double a = 1.0;
    for (m_curr_var = 0; m_curr_var < m_data.size(); m_curr_var++) {
      a *= EpanechnikovPdf(data[m_curr_var], m_data[m_curr_var][i], m_bandwidth[m_curr_var]);
    }
    d += a;
  }

  return (d / m_count[0]);
}

void KernelDensity::DefaultBandwidth()
{
  if (!m_count[m_curr_var]) {
    assert(false);//ASSERT(false, "No data when attempting to configure default bandwidth");
  }

  double x     = m_sum_x[m_curr_var] / m_count[m_curr_var];
  double x2    = m_sum_x2[m_curr_var] / m_count[m_curr_var];
  double sigma = sqrt(x2 - (x * x));
  double bw    = sigma * (pow((3.0 * m_count[m_curr_var] / 4.0), (-1.0 / 5.0)));

  m_bandwidth[m_curr_var]         = bw;
  m_default_bandwidth[m_curr_var] = bw;
}

void KernelDensity::AddSample(double x)
{
  std::vector<double> tmp;
  tmp.push_back(x);
  AddSamples(tmp);
}

void KernelDensity::AddSamples(std::vector<double>& x)
{
  if (!m_data.size()) {
    for (size_t i = 0; i < x.size(); i++) {
      std::vector<double> tmp;
      tmp.push_back(x[i]);
      m_data.push_back(tmp);
      m_sum_x[i]     = x[i];
      m_sum_x2[i]    = x[i] * x[i];
      m_count[i]     = 1;
      m_min[i]       = x[i];
      m_max[i]       = x[i];
      m_bandwidth[i] = -1.0;
    }
  } else {
    for (size_t i = 0; i < x.size(); i++) {
      m_data[i].push_back(x[i]);
      m_sum_x[i] += x[i];
      m_sum_x2[i] += x[i] * x[i];
      m_count[i]++;
      m_min[i]       = x[i] < m_min[i] ? x[i] : m_min[i];
      m_max[i]       = x[i] > m_max[i] ? x[i] : m_max[i];
      m_bandwidth[i] = -1.0;
    }
  }
}

double KernelDensity::EpanechnikovPdf(const double& x, const double& mu, const double& sigma) const
{
  double z = (x - mu) / sigma;
  if (fabs(z) > 1.0) {
    return (0.0);
  }
  return 0.75 * (1.0 - (z * z)) / sigma;
}
