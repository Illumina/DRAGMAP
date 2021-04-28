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

#ifndef SEQUENCES_READ_PAIR_HPP
#define SEQUENCES_READ_PAIR_HPP

#include <array>

#include "sequences/Read.hpp"

namespace dragenos {
namespace sequences {

class ReadPair : public std::array<Read, 2> {
public:
  // ReadPair(Read &&r0, Read &&r1) : std::array<Read, 2>({r0, r1})
  //{
  //}
  // ReadPair(const Read &r0, const Read &r1) : std::array<Read, 2>({r0, r1})
  //{
  //}

  double getLength() const { return (double(at(0).getLength()) + double(at(1).getLength())) / 2; }

  friend std::ostream& operator<<(std::ostream& os, const ReadPair& pair)
  {
    return os << "ReadPair(" << pair[0] << "," << pair[1] << ")" << std::endl;
  }
};

}  // namespace sequences
}  // namespace dragenos

#endif  // #ifndef SEQUENCES_READ_PAIR_HPP
