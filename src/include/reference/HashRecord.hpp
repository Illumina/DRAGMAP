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

#ifndef REFERENCE_HASH_RECORD_HPP
#define REFERENCE_HASH_RECORD_HPP

#include <cstdint>

#include "common/Bits.hpp"

namespace dragenos {
namespace reference {
/**
 ** \brief Encapsulation of the interpretation of binary hashtable data
 **
 ** Only provides low level primitives to identify the type of record and to extract
 ** specific values. It is the responsibility of the user code to ensure that the
 ** extracted values are compatible with the actual type of the record.
 **
 ** Compatibility with hashtable v8.
 **/
class HashRecord {
public:
  /**
   ** \brief Default is HIT. When bits [31:28] are all set, then bits [27:24] identify the record type.
   **
   ** The structure for each type is as follows:
   ** HIT           : ThreadId[63:58] HashBits[57:35] EX[34] LF[33] RC[32] RefPos[31:0]
   ** EMPTY         : 0[63:32] 0xF[31:28] 0x0[27:24] 0[23:0]
   ** HIFREQ        : ThreadId[63:58] HashBits[57:35] EX[34] LF[33] RC[32] 0xF[31:28] 0x1[27:24] 0[23] AL[22] Frequency[21:0]
   ** EXTEND        : ThreadId[63:58] HashBits[57:35] EX[34] LF[33] RC[32] 0xF[31:28] 0x2[27:24] RF[23] AL[22] ExtensionLength[21:18] ExtensionId[17:0]
   ** CHAIN_BEG_MASK: FilterMask[63:32] 0xF[31:28] 0x4[27:24] 0[23:18] ChainPointer[17:0]
   ** CHAIN_BEG_LIST: FilterList4[63:56] FilterList3[55:48] FilterList2[47:40] FilterList1[39:32] 0xF[31:28] 0x5[27:24] 0[23:18] ChainPointer[17:0]
   ** CHAIN_CON_MASK: FilterMask[63:32] 0xF[31:28] 0x6[27:24] 0[23:18] ChainPointer[17:0]
   ** CHAIN_CON_LIST: FilterList4[63:56] FilterList3[55:48] FilterList2[47:40] FilterList1[39:32]  0xF[31:28] 0x7[27:24] 0[23:18] ChainPointer[17:0]
   ** INTERVAL_SL0  : ThreadId[63:58] HashBits[57:35] EX[34] LF[33] 0[32] 0xF[31:28] 0x8[27:24] Length[23:15] Start [14:0]
   ** INTERVAL_SL1  : ThreadId[63:58] HashBits[57:35] EX[34] LF[33] 1[32] 0xF[31:28] 0x8[27:24] Length[23:8] Start [7:0]
   ** INTERVAL_SLE  : ThreadId[63:58] HashBits[57:35] EX[34] LF[33] MSB[32] 0xF[31:28] 0x9[27:24] Exlifts[23:16] Length[15:8] Start [7:0]
   ** INTERVAL_S    : ThreadId[63:58] HashBits[57:35] EX[34] LF[33] MSB[32] 0xF[31:28] 0xA[27:24] Start[23:0]
   ** INTERVAL_L    : ThreadId[63:58] HashBits[57:35] EX[34] LF[33] MSB[32] 0xF[31:28] 0xB[27:24] Length[23:0]
   **
   ** Alt-aware specific records
   ** ALT HIT: Positions exceeding a threshold value
   ** DUMMY HIT: When there isn't ant Pri hit for the alt group, a HIT with POSITION=0 is inserted in the
   *table
   ** Each liftover group has 1+Alt hits followed by one Pri hit. If there isn't any Pri hit for the liftover
   *gropup, a dummy hit is inserted.
   **
   ** Valid combination of INTERVAL_* records in INTERVAL sets:
   ** - INTERVAL_SL0
   ** - INTERVAL_SL1, INTERVAL_S
   ** - INTERVAL_SLE (exlifts > 0, MSB =0)
   ** - INTERVAL_SLE (exlifts > 0, MSB =1), INTERVAL_S
   ** - INTERVAL_SLE (exlifts = 0, MSB =0), INTERVAL_L
   ** - INTERVAL_SLE (exlifts = 0, MSB =1), INTERVAL_S, INTERVAL_L
   ** - INTERVAL_S, INTERVAL_L
   **/
  enum RecordType {
    EMPTY          = 0,
    HIFREQ         = 1,
    EXTEND         = 2,
    REPAIR         = 3,  // obsolete
    CHAIN_BEG_MASK = 4,
    CHAIN_BEG_LIST = 5,
    CHAIN_CON_MASK = 6,
    CHAIN_CON_LIST = 7,
    INTERVAL_SL    = 8,
    INTERVAL_SLE   = 9,
    INTERVAL_S     = 0xA,
    INTERVAL_L     = 0xB,
    HIT            = 16  // strictly greater than any 4-bit op code
  };

  constexpr static unsigned THREAD_ID_START = 58;
  constexpr static unsigned THREAD_ID_BITS  = 6;
  constexpr static unsigned HASH_BITS_START = 35;
  constexpr static unsigned HASH_BITS_BITS  = 23;
  constexpr static unsigned EX_FLAG         = 34;  // primary vs extended seed
  constexpr static unsigned LF_FLAG         = 33;  // Last in Thread
  constexpr static unsigned RC_FLAG         = 32;  // Forward vs Reverse Complement
  constexpr static unsigned RS_FLAG         = 32;  // has random sample (HIFREQ and EXTEND records only
  constexpr static unsigned REFERENCE_POSITION_START = 0;
  constexpr static unsigned REFERENCE_POSITION_BITS  = 32;
  constexpr static unsigned NOT_HIT_START            = 28;
  constexpr static unsigned NOT_HIT_BITS             = 4;
  constexpr static unsigned OP_CODE_START            = 24;
  constexpr static unsigned OP_CODE_BITS             = 4;
  constexpr static unsigned AL_FLAG                  = 22;
  constexpr static unsigned FREQUENCY_START          = 0;
  constexpr static unsigned FREQUENCY_BITS           = 22;
  constexpr static unsigned RF_FLAG                  = 23;
  constexpr static unsigned EXTENSION_LENGTH_START   = 18;
  constexpr static unsigned EXTENSION_LENGTH_BITS    = 4;
  constexpr static unsigned EXTENSION_ID_START       = 0;
  constexpr static unsigned EXTENSION_ID_BITS        = 18;
  constexpr static unsigned CHAIN_POINTER_START      = 0;
  constexpr static unsigned CHAIN_POINTER_BITS       = 18;
  constexpr static unsigned FILTER_MASK_VALUE_BITS = 5;  // number of bits used to calculate the filter value
  constexpr static unsigned FILTER_MASK_START      = 32;
  constexpr static unsigned FILTER_MASK_BITS       = (1U << FILTER_MASK_VALUE_BITS);
  constexpr static unsigned FILTER_LIST_VALUE_BITS = 8;
  constexpr static unsigned FILTER_LIST_COUNT      = FILTER_MASK_BITS / FILTER_LIST_VALUE_BITS;
  constexpr static unsigned MATCH_BITS_START       = EX_FLAG;
  constexpr static unsigned MATCH_BITS_BITS        = THREAD_ID_BITS + HASH_BITS_BITS + 1;
  constexpr static unsigned EXLIFT_START           = 16;
  constexpr static unsigned EXLIFT_BITS            = 8;

  uint64_t getValue() const { return value_; }
  template <unsigned START, unsigned BITS>
  uint64_t getBits() const
  {
    return common::bits::getBits<START, BITS>(value_);
  }
  template <unsigned POSITION>
  bool getFlag() const
  {
    return common::bits::getFlag<POSITION>(value_);
  }
  bool     isHit() const { return 0xF != getBits<NOT_HIT_START, NOT_HIT_BITS>(); }
  bool     isDummyHit() const { return 0 == getPosition(); }
  uint64_t getOpCode() const { return getBits<OP_CODE_START, OP_CODE_BITS>(); }
  uint64_t getType() const { return isHit() ? HIT : getOpCode(); }
  uint8_t  getThreadId() const { return getBits<THREAD_ID_START, THREAD_ID_BITS>(); }
  uint32_t getMatchBits() const { return getBits<MATCH_BITS_START, MATCH_BITS_BITS>(); }
  bool     isLastInThread() const { return getFlag<LF_FLAG>(); }
  bool     isReverseComplement() const { return getFlag<RC_FLAG>(); }
  bool     isMsb() const { return isReverseComplement(); }     // specific meaning for INTERVAL_SLE
  bool     hasCarry() const { return isReverseComplement(); }  // specific meaning for INTERVAL_S
  bool     isExtendedSeed() const { return getFlag<EX_FLAG>(); }
  bool     isAltLiftover() const { return getFlag<AL_FLAG>(); }
  bool     hasRepairRecords() const { return getFlag<RF_FLAG>(); }
  bool     hasRandomSamples() const { return getFlag<RS_FLAG>(); }
  uint32_t getFrequency() const { return getBits<FREQUENCY_START, FREQUENCY_BITS>(); }
  bool     isChainBegin() const
  {
    const auto type = getType();
    return (type == CHAIN_BEG_MASK) || (type == CHAIN_BEG_LIST);
  }
  bool isChainCon() const
  {
    const auto type = getType();
    return (type == CHAIN_CON_MASK) || (type == CHAIN_CON_LIST);
  }
  bool     isChainRecord() const { return isChainBegin() || isChainCon(); }
  uint64_t getChainPointer() const { return getBits<CHAIN_POINTER_START, CHAIN_POINTER_BITS>(); }
  uint8_t  getExlift() const { return getBits<EXLIFT_START, EXLIFT_BITS>(); }

public:
  /// The reference position - if any - is just the 32 LSB of the value
  uint32_t getPosition() const { return static_cast<uint32_t>(value_); }
  /// Hash bits used to match queries to the correct hash records in a bucket
  uint32_t getHashBits() const { return getBits<HASH_BITS_START, HASH_BITS_BITS>(); }
  /// For CHAIN masks
  uint32_t getFilterMask() const { return getBits<FILTER_MASK_START, FILTER_MASK_BITS>(); }
  /// Total extension length (double the length of each wing)
  uint64_t getExtensionLength() const { return getBits<EXTENSION_LENGTH_START, EXTENSION_LENGTH_BITS>(); }
  /// Used to create the key to hash for the next seed extension
  uint64_t getExtensionId() const { return getBits<EXTENSION_ID_START, EXTENSION_ID_BITS>(); }

private:
  uint64_t value_;  // pointer to the record
};

}  // namespace reference
}  // namespace dragenos
#endif  // #ifndef REFERENCE_HASH_RECORD_HPP
