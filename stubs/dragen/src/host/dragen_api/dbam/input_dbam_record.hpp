// Copyright 2016 Edico Genome Corporation. All rights reserved.
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

#ifndef __INPUT_DBAM_RECORD_HPP__
#define __INPUT_DBAM_RECORD_HPP__

#include "dragen_exception.hpp"

//--------------------------------------------------------------------------------adamb
// If the DBAM record is the first of a pair, it includes a sub-record describing
// the insert-size stats to use for rescue alignments and scoring:
//
struct DbamInsertStats {
  uint16_t peMeanInsert = 0;
  uint16_t peMinInsert = 0;
  uint16_t peMaxInsert = 0;
  uint16_t rescueMinInsert = 0;
  uint16_t rescueMaxInsert = 0;
  uint16_t insertSigmaFactor = 0;
  uint8_t  peOrientation : 2;
  uint8_t  unused0 : 6;
  uint8_t  unused1[3] = {0,0,0};

  DbamInsertStats() : peOrientation(0), unused0(0) {}
};


//--------------------------------------------------------------------------------adamb
// First thing in an input DBAM record is a header, always the same for all records:
//
struct InputDbamHeader {
  uint16_t recordQuadwords;
  uint16_t flag;
  uint8_t  sequenceLen[3];
  uint8_t  anchorDoublewords;
  uint32_t fragmentIndex;
  uint32_t qnameHash;

  InputDbamHeader() : recordQuadwords(0), flag(0), anchorDoublewords(0), fragmentIndex(0), qnameHash(0)
  {
    setSequenceLen(0);
  }
//
//  // Calculate the size of the record when the exact size of the blob, bloblen, is known.
//  // bloblen must be padded to a multiple of 8 bytes.
//  static uint32_t calculateRecordLen(
//      const uint16_t flag, const uint32_t seqlen, const uint32_t bloblen, const uint32_t numAnchorRecords)
//  {
//    const bool firstOfPair = ((flag & InputDbamUtils::pairedRead1Mask) == InputDbamUtils::pairedRead1Mask);
//    return sizeof(InputDbamHeader) + (firstOfPair ? sizeof(DbamInsertStats) : 0) + pad8(seqlen) + bloblen +
//           pad8(numAnchorRecords * sizeof(uint32_t));
//  }
//
//  static uint32_t calculateRecordLen(
//      const uint16_t     flag,
//      const uint32_t     seqlen,
//      const uint32_t     namelen,
//      const uint32_t     tagbloblen,
//      const uint32_t     numAnchorRecords,
//      const UMI::UMIType umiType)
//  {
//    // TODO - add anchors
//    return calculateRecordLen(
//        flag, seqlen, DbamBlobSection::getSize(namelen, tagbloblen, umiType), numAnchorRecords);
//  }
//
//  uint32_t setRecordSize(
//      const uint16_t     flag,
//      const uint32_t     seqlen,
//      const uint32_t     namelen,
//      const uint32_t     tagbloblen,
//      const uint32_t     numAnchorRecords,
//      const UMI::UMIType umiType)
//  {
//    const uint32_t reclen = calculateRecordLen(flag, seqlen, namelen, tagbloblen, numAnchorRecords, umiType);
//    recordQuadwords       = reclen / sizeof(uint64_t);
//    return reclen;
//  }
//
//  void setRecordLen(const uint32_t reclen) { recordQuadwords = reclen / sizeof(uint64_t); }
//
//  uint32_t getRecordLen() const { return recordQuadwords * sizeof(uint64_t); }
//
  void setSequenceLen(const uint32_t len) { memcpy(sequenceLen, &len, 3); }

  uint32_t getSequenceLen() const
  {
    uint32_t x = 0;
    memcpy(&x, sequenceLen, 3);
    return x;
  }

};


class InputDbamRecord
{
  InputDbamHeader dbh_;
  DbamInsertStats insertStats_;

public:
  uint8_t  peStatsInterval = 0;
  InputDbamRecord (std::size_t readLen = -1)
  {
    dbh_.setSequenceLen(readLen);
  }

  InputDbamRecord (uint8_t* dummy)
  {
    // ReadGroupInsertStats::saveForRemapping does this weird malloc, so...
    free(dummy);
  }
//
//  InputDbamRecord(InputDbamRecord&& that) :
//    alignment_(that.alignment_),
//    dbh_(that.dbh_),
//    insertSstats_(that.insertSstats_),
//    peStatsInterval(that.peStatsInterval)
//  {
//  }

  DbamInsertStats* getInsertStats()
  {
    if (isFirstOfPair())
      return &insertStats_;
    else
      return 0;
  }


  bool isFirstOfPair() const {
    return true;
//    assert(pair_);
//    return pair_->isFirstInTemplate();
  }

  InputDbamHeader* getHeader() { return &dbh_; }

  InputDbamRecord& getBlobs() { return *this; }

  InputDbamRecord* getFixedBlob() { return this; }

  uint32_t getRecordLen() const { return 0; }

  bool hasMate() const {
    return true;
//    assert(pair_);
//    return pair_->hasMultipleSegments();
  }

  uint8_t* getRawData() const { return 0; }
};

#endif
