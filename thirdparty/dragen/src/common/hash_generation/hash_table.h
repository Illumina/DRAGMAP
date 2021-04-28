#ifndef _hash_table_h
#define _hash_table_h

#include <inttypes.h>
#include "gen_hash_table.h"
#include "liftover.h"

#define BASE_PAD 0
#define BASE_A 1
#define BASE_C 2
#define BASE_M_AC 3
#define BASE_G 4
#define BASE_R_AG 5
#define BASE_S_CG 6
#define BASE_V_ACG 7
#define BASE_T 8
#define BASE_W_AT 9
#define BASE_Y_CT 10
#define BASE_H_ACT 11
#define BASE_K_GT 12
#define BASE_D_AGT 13
#define BASE_B_CGT 14
#define BASE_N 15

#define SINGLE_BASE(x) ((x) == BASE_A || (x) == BASE_C || (x) == BASE_G || (x) == BASE_T)

// This array contains values 0-15 for letters defined above, or lowercase equivalents, and -1 for other
// characters. E.g. ENCODE_BASE['a'] = 1, ENCODE_BASE['R'] = 5, ENCODE_BASE[' '] = -1.
extern const char ENCODE_BASE[256];

// Convert 4-bit code to base value 0-3 for ACGT or 0 otherwise
extern const uint8_t baseCodeBase[16];
// Number of '1' bits in 4-bit codes 0-15
extern const uint8_t baseCodeNumBases[16];
// Number of '1' bits in 4-bit codes 0-15 (minimum of 1)
extern const uint8_t baseCodeNumBasesMin1[16];
// List of positions of '1' bits in 4-bit codes 0-15
extern const uint8_t baseCodeBaseList[16][4];

#define MAX_REF_SEQS ((1 << 24) - 1)
#define THINNING_MAX_PERIOD 16
#define MAX_SEED_HIT_FREQ 256
#define MAX_SEED_INDEXES 0xF0000000
#define TARGET_OCCUPANCY 0.75
#define THRESH_OCCUPANCY 0.80
#define MAPPER_REF_POS_BITS 36
#define KEY_ANCHOR_OFFSET 16
#define MIN_ANCHOR_BIN_BITS 8
#define REF_BASES_PER_SJ_RESERVE_BYTE 32
#define REF_SEQ_TRIM_GRAN 256
#define MAX_HASH_TABLE_CHUNKS 32

// Opcodes include the escape nibble
#define HASH_OPC_EMPTY 0xF0
// #define HASH_OPC_HIFREQ 0xF1  (no longer used)
#define HASH_OPC_EXTEND 0xF2
#define HASH_OPC_CHAIN_BEG_MASK 0xF4
#define HASH_OPC_CHAIN_BEG_LIST 0xF5
#define HASH_OPC_CHAIN_CON_MASK 0xF6
#define HASH_OPC_CHAIN_CON_LIST 0xF7
#define HASH_OPC_INTERVAL_SL 0xF8
#define HASH_OPC_INTERVAL_SLE 0xF9
#define HASH_OPC_INTERVAL_S 0xFA
#define HASH_OPC_INTERVAL_L 0xFB
#define HASH_OPC_CHAIN_LIST_FLAG 0x01

#define HASH_OPC_SPECIAL_HIT 0xFC  // Internal use

#define HASH_OPC_IS_INTVL(x) (((x)&0xFC) == HASH_OPC_INTERVAL_SL)
#define HASH_OPC_IS_CHAIN(x) (((x)&0xFC) == HASH_OPC_CHAIN_BEG_MASK)
#define HASH_OPC_IS_BEG(x) (((x)&0xFE) == HASH_OPC_CHAIN_BEG_MASK)
#define HASH_OPC_IS_CON(x) (((x)&0xFE) == HASH_OPC_CHAIN_CON_MASK)
#define HASH_OPC_IS_LIST(x) (((x)&0xFD) == HASH_OPC_CHAIN_BEG_LIST)
#define HASH_OPC_IS_MASK(x) (((x)&0xFD) == HASH_OPC_CHAIN_BEG_MASK)
#define HASH_OPC_IS_HIT(x) (((x) < HASH_OPC_EMPTY) | ((x) >= HASH_OPC_SPECIAL_HIT))
#define HASH_OPC_IS_NORM(x) ((x) < HASH_OPC_EMPTY)
#define HASH_OPC_IS_SPEC(x) ((x) >= HASH_OPC_SPECIAL_HIT)

#define HASH_BUCKET_BYTES_LOG2 6
#define HASH_BUCKET_BYTES (1 << HASH_BUCKET_BYTES_LOG2)
#define HASH_RECORD_BYTES_LOG2 3
#define HASH_RECORD_BYTES (1 << HASH_RECORD_BYTES_LOG2)
#define HASH_EXTRA_BITS 3
#define MAX_PROBES (1 << HASH_EXTRA_BITS)
#define HASH_RECORDS_PER_BUCKET_LOG2 (HASH_BUCKET_BYTES_LOG2 - HASH_RECORD_BYTES_LOG2)
#define HASH_RECORDS_PER_BUCKET (HASH_BUCKET_BYTES / HASH_RECORD_BYTES)
#define HASH_RECORD_EXT_ID_BITS 18
#define HASH_RECORD_HASH_BITS 23
#define EXT_ID_HASH_BITS_MIN 3
#define EXT_ID_HASH_BITS_MAX 24
#define HASH_THREADS_LOG2 HASH_RECORDS_PER_BUCKET_LOG2
#define HASH_THREADS HASH_RECORDS_PER_BUCKET
#define BUCKET_THREAD_MASK (HASH_THREADS - 1)
#define HASH_BYTE_ADDR_START 19
#define HASH_BUCKET_ADDR_START (HASH_BYTE_ADDR_START + HASH_BUCKET_BYTES_LOG2)
#define ADDR_THREAD_ID_START 3
#define THREAD_ID_BITS 6
#define INDEPENDENT_ADDR_BITS 30
#define INDEPENDENT_HASH_BITS (INDEPENDENT_ADDR_BITS + HASH_BYTE_ADDR_START)
#define MAX_SEC_CRC_BITS INDEPENDENT_HASH_BITS
#define MAX_PRI_SEED_LENGTH 30
#define MAX_NET_SEED_EXTENSION 128
#define MAX_EXTENDED_LENGTH (MAX_PRI_SEED_LENGTH + MAX_NET_SEED_EXTENSION)
#define MAX_SEED_EXTENSION_INCR 12
#define SEC_CRC_BITS_MINUS_EXT_ID_HASH_BITS (HASH_RECORD_EXT_ID_BITS + 2 * MAX_SEED_EXTENSION_INCR)

#define MAX_WRAP_BYTES_LOG2 15
#define MAX_WRAP_BYTES (1 << MAX_WRAP_BYTES_LOG2)
#define CHAIN_BLOCK_BUCKETS_LOG2 18
#define CHAIN_BLOCK_BUCKETS (1 << CHAIN_BLOCK_BUCKETS_LOG2)
#define CHAIN_LIST_HASH_BITS 8
#define CHAIN_MASK_HASH_BITS 5
#define CHAIN_LIST_HASH_MASK ((1 << CHAIN_LIST_HASH_BITS) - 1)
#define CHAIN_MASK_HASH_MASK ((1 << CHAIN_MASK_HASH_BITS) - 1)

typedef struct {
  uint32_t lo_fields : 24;
  uint32_t opcode : 8;
  uint32_t hi_fields : 32;
} hashrec_general_t;

typedef struct {
  uint32_t lo_fields : 24;
  uint32_t opcode : 8;
  uint32_t rc : 1;
  uint32_t lf : 1;
  uint32_t ex : 1;
  uint32_t hash_bits : 23;
  uint32_t thread_id : 6;
} hashrec_matchable_t;

typedef struct {
  uint32_t lo_fields : 24;
  uint32_t opcode : 8;
  uint32_t rc : 1;
  uint32_t lf : 1;
  uint32_t match_bits : 30;
} hashrec_match_bits_t;

typedef struct {
  uint32_t pos : 32;
  uint32_t rc : 1;
  uint32_t lf : 1;
  uint32_t ex : 1;
  uint32_t hash_bits : 23;
  uint32_t thread_id : 6;
} hashrec_hit_t;

typedef struct {
  uint32_t extend_id : 18;
  uint32_t extend_len : 4;
  uint32_t al : 1;
  uint32_t rf : 1;  // RF no longer used
  uint32_t opcode : 8;
  uint32_t rs : 1;  // RS no longer used
  uint32_t lf : 1;
  uint32_t ex : 1;
  uint32_t hash_bits : 23;
  uint32_t thread_id : 6;
} hashrec_extend_t;

typedef struct {
  uint32_t start : 15;
  uint32_t length : 9;
  uint32_t opcode : 8;
  uint32_t fmt : 1;
  uint32_t lf : 1;
  uint32_t ex : 1;
  uint32_t hash_bits : 23;
  uint32_t thread_id : 6;
} hashrec_interval_sl0_t;

typedef struct {
  uint32_t start : 8;
  uint32_t length : 16;
  uint32_t opcode : 8;
  uint32_t fmt : 1;
  uint32_t lf : 1;
  uint32_t ex : 1;
  uint32_t hash_bits : 23;
  uint32_t thread_id : 6;
} hashrec_interval_sl1_t;

typedef struct {
  uint32_t start : 8;
  uint32_t length : 8;
  uint32_t exlifts : 8;
  uint32_t opcode : 8;
  uint32_t msb : 1;
  uint32_t lf : 1;
  uint32_t ex : 1;
  uint32_t hash_bits : 23;
  uint32_t thread_id : 6;
} hashrec_interval_sle_t;

typedef struct {
  uint32_t start : 24;
  uint32_t opcode : 8;
  uint32_t carry : 1;
  uint32_t lf : 1;
  uint32_t ex : 1;
  uint32_t hash_bits : 23;
  uint32_t thread_id : 6;
} hashrec_interval_st_t;

typedef struct {
  uint32_t length : 24;
  uint32_t opcode : 8;
  uint32_t rsvd : 1;
  uint32_t lf : 1;
  uint32_t ex : 1;
  uint32_t hash_bits : 23;
  uint32_t thread_id : 6;
} hashrec_interval_ln_t;

typedef struct {
  uint32_t chain_ptr : 18;
  uint32_t chain_pad : 6;
  uint32_t opcode : 8;
  uint32_t filter : 32;
} hashrec_chain_t;

typedef struct {
  uint32_t chain_ptr : 18;
  uint32_t chain_pad : 6;
  uint32_t opcode : 8;
  uint8_t  filter_list[4];
} hashrec_chain_list_t;

typedef union {
  hashrec_general_t      general;
  hashrec_matchable_t    matchable;
  hashrec_match_bits_t   match_bits;
  hashrec_hit_t          hit;
  hashrec_extend_t       extend;
  hashrec_interval_sl0_t interval_sl0;
  hashrec_interval_sl1_t interval_sl1;
  hashrec_interval_sle_t interval_sle;
  hashrec_interval_st_t  interval_st;
  hashrec_interval_ln_t  interval_ln;
  hashrec_chain_t        chain;
  hashrec_chain_list_t   chain_list;
  uint64_t               qword;
} hashrec_t;

#define LIFT_CODE_NONE 0
#define LIFT_CODE_ALT 1
#define LIFT_CODE_PRI 2
#define LIFT_CODE_DIF_PRI 3

typedef struct {
  uint32_t pos : 32;
  uint32_t rc : 1;
  uint32_t lift_code : 2;
  uint32_t lift_group : 28;
  uint32_t literal : 1;
} extend_hit_t;

#define HASH_REC_EMPTY_QWORD 0x00000000F0000000ull
#define HASH_REC_CHAIN_TERM_QWORD 0x00000000F6000000ull
#define HASH_REC_POS_SPECIAL 0xFC000000
#define HASH_REC_SPECIAL_MASK 0x03FFFFFF

char* buildHashTable(
    hashTableConfig_t* config,
    uint8_t*           refSeq,
    uint64_t*          refCodeHist,
    uint8_t*           refMask,
    uint8_t*           refCode,
    uint8_t*           altMatches

);

// Copy N 2-bit bases from a possibly unaligned source to an aligned destination.
// srcOffset = 0-3 indicates the position of the first base in the first source byte.
// srcOffset may have higher bits set, which are ignored.
void copyBases(void* dst, void* src, int len, int srcOffset);

// Reverse-complement N 2-bit bases, either in-place or to a separate destination.
// The input and output are aligned.
void revComp(void* dst, void* src, int len);

// Functions to extract segments of hash bits
uint64_t qwordExtractBits(uint64_t hash, int start, int len);
uint64_t hashRecExtractHashBits(hashrec_t rec, int start, int len);

// Function returns a readable string for a byte length, e.g. "22.25 GB".
// No need to free this string; the function manages a pool of them.
char* bytesReadable(uint64_t bytes);

#endif
