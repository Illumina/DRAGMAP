#pragma once

class InputDbamRecord;

//--------------------------------------------------------------------------------adamb
// Interface for classes that can pump InputDbamRecord's back through the
// mapper/aligner. Implement it in the various ReadSenders so they can
// use the StatsInitRemapper
//
class InputDbamRemapper {
public:
  virtual ~InputDbamRemapper(){};

  virtual void sendRead(InputDbamRecord& dbr, const bool lastRecord) = 0;
};
