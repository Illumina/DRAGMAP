#ifndef __CRC32_HW_H__
#define __CRC32_HW_H__
#include <inttypes.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

uint32_t crc32c_hw(uint32_t crc, const void* buf, size_t len);
bool     machine_has_sse42();

#ifdef __cplusplus
}
#endif

#endif
