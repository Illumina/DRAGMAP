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

#include <emmintrin.h>

#include <cassert>
#include <sstream>
#include <vector>

#include "common/Exceptions.hpp"
#include "sequences/Read.hpp"

namespace dragenos {
namespace sequences {

Read& Read::operator=(Read&& that)
{
  id_       = that.id_;
  position_ = that.position_;
  //  name_ = std::move(that.name_);
  //  bases_ = std::move(that.bases_);
  //  qualities_ = std::move(that.qualities_);
  name_.swap(that.name_);
  bases_.swap(that.bases_);
  rcBases_.swap(that.rcBases_);
  qualities_.swap(that.qualities_);
  return *this;
}

void reverseComplement4bpb(const Read::Bases& bases, Read::Bases& rcBases)
{
  const char A = 1;
  const char C = 2;
  const char G = 4;
  const char T = 8;
  const char M = A | C;
  const char R = A | G;
  const char S = C | G;
  const char V = A | C | G;
  const char W = A | T;
  const char Y = C | T;
  const char H = A | C | T;
  const char K = G | T;
  const char D = A | G | T;
  const char B = C | G | T;
  const char N = A | C | G | T;
  //static const char bases[] = {'N', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'};
  //static const char bases[] = {'N', 'T', 'G', 'K', 'C', 'Y', 'S', 'B', 'A', 'W', 'R', 'D', 'M', 'H', 'V', 'N'};
  static const char rc[] = {N, T, G, K, C, Y, S, B, A, W, R, D, M, H, V, N};
  rcBases.resize(bases.size());
  std::transform(bases.rbegin(), bases.rend(), rcBases.begin(), [](char c) { return rc[c & 0xf]; });
}

void Read::init(Name&& name, Bases&& bases, Qualities&& qualities, uint64_t id, unsigned position)
{
  id_       = id;
  position_ = position;
  //  name_ = std::move(name);
  name_.swap(name);

  //  bases_ = std::move(bases);
  bases_.swap(bases);
  //  qualities_ = std::move(qualities);
  qualities_.swap(qualities);
  reverseComplement4bpb(bases_, rcBases_);
}

std::ostream& operator<<(std::ostream& os, const __m128i& i128)
{
  for (size_t i = 0; sizeof(__m128i) != i; ++i) {
    os << std::setfill('0') << std::hex << std::setw(2) << (unsigned int)(((unsigned char*)&i128)[i]) << " ";
  }

  return os;
}

}  // namespace sequences
}  // namespace dragenos
