#ifndef __FAST_NONVECTOR_CRC32C_H__
#define __FAST_NONVECTOR_CRC32C_H__

#ifdef __cplusplus
extern "C" {
#endif

uint32_t fastCrc32c1(uint32_t crc, const uint8_t* data);
uint32_t fastCrc32c1f(const uint8_t* data);
uint32_t fastCrc32c2(uint32_t crc, const uint8_t* data);
uint32_t fastCrc32c2f(const uint8_t* data);
uint32_t fastCrc32c4(uint32_t crc, const uint8_t* data);
uint32_t fastCrc32c4f(const uint8_t* data);
uint32_t fastCrc32c8(uint32_t crc, const uint8_t* data);
uint32_t fastCrc32c8f(const uint8_t* data);
uint32_t fastCrc32c16(uint32_t crc, const uint8_t* data);
uint32_t fastCrc32c16f(const uint8_t* data);

#ifdef __cplusplus
}
#endif

#endif  // __FAST_NONVECTOR_CRC32C_H__
