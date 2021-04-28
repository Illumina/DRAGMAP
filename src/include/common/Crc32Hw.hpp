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

#ifndef COMMON_CRC32_HW_HPP
#define COMMON_CRC32_HW_HPP

#include <cstddef>
#include <cstdint>

namespace dragenos {
namespace common {

uint32_t crc32c_hw(uint32_t crc, const void* buf, std::size_t len);
bool     machine_has_sse42();

}  // namespace common
}  // namespace dragenos

#endif  //#ifndef COMMON_CRC32_HW_HPP
