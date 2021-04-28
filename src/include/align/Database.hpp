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

#ifndef ALIGN_DATABASE_HPP
#define ALIGN_DATABASE_HPP

#include <vector>

namespace dragenos {
namespace align {

/**
 ** \brief Proxy for the reference for the Smith-waterman algorithm
 **/
class Database : public std::vector<unsigned char> {
  typedef std::vector<unsigned char> BaseT;

public:
  // conveniently forward all valid costructors from std::vector<char>
  template <typename... Args>
  Database(Args... args) : std::vector<unsigned char>(std::forward<Args>(args)...)
  {
  }
  //
private:
  friend std::ostream& operator<<(std::ostream& os, const Database& db)
  {
    for (const auto& d : db) {
      os << sequences::Read::decodeBase(d);
    }
    return os;
  }
};

}  // namespace align
}  // namespace dragenos

#endif  // #ifndef ALIGN_DATABASE_HPP
