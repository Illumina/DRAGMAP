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

#include "reference/ReferenceSequence.hpp"

#include <boost/format.hpp>

namespace dragenos {
namespace reference {

char ReferenceSequence::decodeBase(unsigned char base)
{
  const static std::string translate{"PACMGRSVTWYHKDBN"};
  return translate[base & 0xF];
}

unsigned char ReferenceSequence::translateTo2bpb(unsigned char base4bpb)
{
  const static std::array<unsigned char, 16> translate{0, 0, 1, 0, 2, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0};
  return translate[base4bpb & 0xF];
}

unsigned char ReferenceSequence::translateToR2bpb(unsigned char base4bpb)
{
  const static std::array<unsigned char, 16> translate{0, 3, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  return translate[base4bpb & 0xF];
}

}  // namespace reference
}  // namespace dragenos
