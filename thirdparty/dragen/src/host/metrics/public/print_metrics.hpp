// Copyright 2018 Edico Genome Corporation. All rights reserved.
// Copyright 2018-2020 Illumina, Inc. All rights reserved.

#ifndef __PRINT_M__
#define __PRINT_M__

#include "elapsed_timer.hpp"

#include <climits>
#include <iomanip>
#include <iostream>
#include <sstream>

#include <boost/format.hpp>
#include <boost/iostreams/device/null.hpp>

//--------------------------------------------------------------------------------
// PrintMetrics - a class for writing Dragen metrics to DRAGEN log and metrics CSVs
// we open and close the file streams outside of the class

class PrintMetrics {
private:
  std::ostream& m_stdOut;
  std::ostream& m_metricsCSV;

protected:
  std::string m_sample;
  std::string m_section;

private:
  uint64_t m_denominator;
  uint64_t m_scaleFactor;
  uint8_t  m_colWidth;
  bool     m_tumorNormalMode = false;

public:
  PrintMetrics(std::ostream& other_stdOut, std::ostream& other_metricsCSV, uint8_t other_colWidth)
    : m_stdOut(other_stdOut),
      m_metricsCSV(other_metricsCSV),
      m_sample(""),
      m_section(""),
      m_denominator(1),
      m_scaleFactor(1),
      m_colWidth(other_colWidth),
      m_tumorNormalMode(false)
  {
  }

  PrintMetrics(std::ostream& other_stdOut = std::cout)
    : m_stdOut(other_stdOut),
      m_metricsCSV(m_stdOut),
      m_sample(""),
      m_section(""),
      m_denominator(1),
      m_scaleFactor(1),
      m_colWidth(80),
      m_tumorNormalMode(false)
  {
  }

  virtual ~PrintMetrics() = default;

  void reset()
  {
    m_sample      = "";
    m_section     = "";
    m_denominator = 1;
    m_scaleFactor = 1;
  }

  void setScaleFactor(uint64_t other_scaleFactor) { m_scaleFactor = other_scaleFactor; }

  void setDenominator(uint64_t other_denumerator) { m_denominator = other_denumerator; }

  virtual void setSampleName(
      const std::string& other_sample, bool sj = false, const bool setPrintedSampleName = true)
  {
    (void)sj;  // Unused in base class
    if (not setPrintedSampleName) return;
    m_sample = other_sample;
  }

  virtual void setSection(const std::string& other_section) { m_section = other_section; }

  void setTumorNormalMode(const bool tnMode) { m_tumorNormalMode = tnMode; }

  // Needed to identify Tumor-Normal context in the SV caller's script-like conditions
  bool getTumorNormalMode() const { return m_tumorNormalMode; }

  enum { SECTION_W = 31, DESC1_W = 66, DESC2_W = 15, VALUE1_W = 15, VALUE2_W = 25 };

  // only print description and metric value
  virtual void printCount(const std::string& metricDescription, uint64_t count)
  {
    m_stdOut << std::left << std::setw(SECTION_W) << m_section.c_str() << std::setw(m_colWidth)
             << m_sample.c_str() << std::setw(DESC1_W) << metricDescription.c_str() << std::setw(VALUE1_W)
             << count << std::setw(VALUE2_W) << std::endl;

    m_metricsCSV << m_section.c_str() << "," << m_sample.c_str() << "," << metricDescription.c_str() << ","
                 << count << std::endl;
  }

  // --------------------------------------------------------------------------
  // output histogram row. If "ceiling" is INT_MAX, it is output as "inf"
  virtual void printHistFloat(
      const std::string& metricDescription,
      uint64_t           floor,
      uint64_t           ceiling,
      const float        ratio,
      int                precision = 2)
  {
    std::string myFloor = boost::str(boost::format("%3s") % std::to_string(floor));
    std::string myCeiling;
    if (ceiling == INT_MAX) {
      myCeiling = std::string(" inf)");
    } else {
      myCeiling = boost::str(boost::format("%3s") % ceiling) + std::string("x)");
    }

    m_stdOut << std::left << std::setw(SECTION_W) << m_section.c_str() << std::setw(m_colWidth)
             << m_sample.c_str() << std::setw(DESC1_W)
             << (metricDescription.c_str() + std::string("[") + myFloor + std::string("x:") + myCeiling)
             << std::setprecision(precision) << std::fixed << ratio << std::endl;

    m_metricsCSV << m_section.c_str() << "," << m_sample.c_str() << "," << metricDescription.c_str() << "["
                 << myFloor << "x:" << myCeiling << "," << std::setprecision(precision) << std::fixed << ratio
                 << std::endl;
  }

  // Sample Category Description Ratio
  // float ratio causes confusing roundings in printed results. Please keep it double.
  virtual void printFloat(const std::string& metricDescription, const double ratio, int precision = 2)
  {
    m_stdOut << std::left << std::setw(SECTION_W) << m_section.c_str() << std::setw(m_colWidth)
             << m_sample.c_str() << std::setw(DESC1_W) << metricDescription.c_str()
             << std::setprecision(precision) << std::fixed << ratio << std::endl;

    m_metricsCSV << m_section.c_str() << "," << m_sample.c_str() << "," << metricDescription.c_str() << ","
                 << std::setprecision(precision) << std::fixed << ratio << std::endl;
  }

  virtual void printNonzeroDuration(const std::string& msg, uint64_t duration_ms)
  {
    if (!duration_ms) return;
    printDuration(msg, duration_ms);
  }

  // Prints the time, given the total duration in milliseconds
  virtual void printDuration(const std::string& msg, uint64_t duration_ms)
  {
    const std::string ss_time    = timeToFormattedString(duration_ms);
    const std::string ss_seconds = timeToSeconds(duration_ms);

    m_stdOut << std::left << std::setw(SECTION_W) << m_section.c_str() << std::setw(m_colWidth)
             << m_sample.c_str() << std::setw(DESC1_W) << msg << std::setw(VALUE1_W) << ss_time << ss_seconds
             << std::endl;

    m_metricsCSV << m_section << "," << m_sample << "," << msg << "," << ss_time << "," << ss_seconds
                 << std::endl;
  }

  virtual void printString(const std::string& metricDescription, float ratio)
  {
    m_stdOut << std::left << std::setw(SECTION_W) << m_section << std::setw(m_colWidth) << m_sample
             << std::setw(DESC1_W) << metricDescription << std::setprecision(2) << std::fixed << ratio;

    m_metricsCSV << std::left << m_section << "," << m_sample << "," << metricDescription << ","
                 << std::setprecision(2) << std::fixed << ratio;
  }

  // Section, Sample, String1, String2, MetricValue
  virtual void printString(
      const std::string& metricDescription1,
      const std::string& metricDescription2,
      float              value,
      int                precision)
  {
    m_stdOut << std::left << std::setw(SECTION_W) << m_section << std::setw(m_colWidth) << m_sample
             << std::setw(DESC1_W) << metricDescription1 << std::setw(VALUE1_W) << metricDescription2
             << std::setprecision(precision) << std::fixed << value << std::endl;

    m_metricsCSV << m_section << "," << m_sample << "," << metricDescription1 << "," << metricDescription2
                 << "," << std::setprecision(precision) << std::fixed << value << std::endl;
  }

  // Section, Sample, string, value
  virtual void printString(const std::string& metricDescription, std::string value)
  {
    m_stdOut << std::left << std::setw(SECTION_W) << m_section << std::setw(m_colWidth) << m_sample
             << std::setw(DESC1_W) << metricDescription << value << std::endl;

    m_metricsCSV << std::left << m_section << "," << m_sample << "," << metricDescription << "," << value
                 << std::endl;
  }

  // SampleRead GroupDescription Count/ ratioPercentage
  virtual void printCountAndPercentage(
      const std::string& metricDescription, uint64_t count, uint64_t denominator, int precision)
  {
    float percent;
    if ((count > 0) && (denominator > 0)) {
      percent = 100.0 * static_cast<float>(count) / denominator;
    } else {
      percent = 0.0;
    }

    m_stdOut << std::left << std::setw(SECTION_W) << m_section.c_str() << std::setw(m_colWidth)
             << m_sample.c_str() << std::setw(DESC1_W) << metricDescription.c_str() << std::setw(VALUE1_W)
             << count << std::setw(VALUE2_W) << std::setprecision(precision) << std::fixed << percent
             << std::endl;

    m_metricsCSV << m_section.c_str() << "," << m_sample.c_str() << "," << metricDescription.c_str() << ","
                 << count << "," << std::setprecision(precision) << std::fixed << percent << std::endl;
  }

  inline virtual void printCountAndPercentage(
      const std::string& metricDescription, uint64_t count, uint64_t denominator)
  {
    printCountAndPercentage(metricDescription, count, denominator, 2);
  }

  // SampleRead GroupDescription Count/ ratioPercentage
  inline virtual void printCountAndPercentage(const std::string& metricDescription, uint64_t count)
  {
    printCountAndPercentage(metricDescription, count, m_denominator);
  }

  virtual void printPercentage(
      const std::string& metricDescription, uint64_t count, uint64_t denominator, int precision)
  {
    float percent;
    if ((count > 0) && (denominator > 0)) {
      percent = 100.0 * static_cast<float>(count) / denominator;
    } else {
      percent = 0.0;
    }

    m_stdOut << std::left << std::setw(SECTION_W) << m_section.c_str() << std::setw(m_colWidth)
             << m_sample.c_str() << std::setw(DESC1_W) << metricDescription.c_str() << std::setw(VALUE1_W)
             << std::setprecision(precision) << std::fixed << percent << std::endl;

    m_metricsCSV << m_section.c_str() << "," << m_sample.c_str() << "," << metricDescription.c_str() << ","
                 << std::setprecision(precision) << std::fixed << percent << std::endl;
  }

  // Essentially like printCountAndPercentage, but without the 100 percentage multiplier and with a
  // potentially floating point normalization factor
  virtual void printCountAndNormalized(
      const std::string& metricDescription, uint64_t value, double normFactor, int precision = 2)
  {
    double normValue = ((value > 0) && (normFactor > 0)) ? value / normFactor : 0.0;

    m_stdOut << std::left << std::setw(SECTION_W) << m_section.c_str() << std::setw(m_colWidth)
             << m_sample.c_str() << std::setw(DESC1_W) << metricDescription.c_str() << std::setw(VALUE1_W)
             << value << std::setw(VALUE2_W) << std::setprecision(precision) << std::fixed << normValue
             << std::endl;

    m_metricsCSV << m_section.c_str() << "," << m_sample.c_str() << "," << metricDescription.c_str() << ","
                 << value << "," << std::setprecision(precision) << std::fixed << normValue << std::endl;
  }

  // Essentially like printCountAndPercentage, but taking floating point values and not x 100
  virtual void printFloatAndNormalized(
      const std::string& metricDescription, double value, double normFactor, int precision)
  {
    double normValue = ((value > 0) && (normFactor > 0)) ? value / normFactor : 0.0;

    m_stdOut << std::left << std::setw(SECTION_W) << m_section.c_str() << std::setw(m_colWidth)
             << m_sample.c_str() << std::setw(DESC1_W) << metricDescription.c_str() << std::setw(VALUE1_W)
             << std::setprecision(precision) << std::fixed << value << std::setw(VALUE2_W)
             << std::setprecision(precision) << std::fixed << normValue << std::endl;

    m_metricsCSV << m_section.c_str() << "," << m_sample.c_str() << "," << metricDescription.c_str() << ","
                 << std::setprecision(precision) << std::fixed << value << "," << std::setprecision(precision)
                 << std::fixed << normValue << std::endl;
  }

  // SampleRead GroupDescription Count/ ratioPercentage
  virtual void printRatio(
      const std::string& metricDescription, float numerator, float denominator, int precision)
  {
    float percent;
    if (denominator == 0) {
      m_stdOut << std::left << std::setw(SECTION_W) << m_section.c_str() << std::setw(m_colWidth)
               << m_sample.c_str() << std::setw(DESC1_W) << metricDescription.c_str() << "inf" << std::endl;

      m_metricsCSV << m_section.c_str() << "," << m_sample.c_str() << "," << metricDescription.c_str() << ","
                   << "inf" << std::endl;
    } else {
      percent = static_cast<float>(numerator) / denominator;
      m_stdOut << std::left << std::setw(SECTION_W) << m_section.c_str() << std::setw(m_colWidth)
               << m_sample.c_str() << std::setw(DESC1_W) << metricDescription.c_str() << std::setw(VALUE2_W)
               << std::setprecision(precision) << std::fixed << percent << std::endl;

      m_metricsCSV << m_section.c_str() << "," << m_sample.c_str() << "," << metricDescription.c_str() << ","
                   << std::setprecision(precision) << std::fixed << percent << std::endl;
    }
  }

  // print NA when stats is not applicable
  virtual void printNA(const std::string& metricDescription)
  {
    m_stdOut << std::left << std::setw(SECTION_W) << m_section.c_str() << std::setw(m_colWidth)
             << m_sample.c_str() << std::setw(DESC1_W) << metricDescription.c_str() << "NA" << std::endl;

    m_metricsCSV << m_section.c_str() << "," << m_sample.c_str() << "," << metricDescription.c_str() << ","
                 << "NA" << std::endl;
  }

  void printCleanLine() const { m_stdOut << std::endl; }
};

#endif

