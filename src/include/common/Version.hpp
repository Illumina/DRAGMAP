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
#ifndef COMMON_VERSION_HPP
#define COMMON_VERSION_HPP

#include <sstream>
#include <string>

#define STRINGIFY(s) XSTRINGIFY(s)
#define XSTRINGIFY(s) #s

namespace dragenos {
namespace common {

struct Version {
  static std::string string() { return std::string(STRINGIFY(VERSION_STRING)); }
};

}  // namespace common
}  // namespace dragenos

#endif  // #ifndef COMMON_VERSION_HPP
