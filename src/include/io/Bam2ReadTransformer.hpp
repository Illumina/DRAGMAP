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

#pragma once

#include <boost/range.hpp>

#include "bam/Bam.hpp"
#include "sequences/Read.hpp"

namespace dragenos {
namespace io {

class BamToReadTransformer {
public:
  sequences::Read::Name      tmpName_;
  sequences::Read::Bases     tmpBases_;
  sequences::Read::Qualities tmpQscores_;
  static const char          DEFAULT_QNAME_DELIM = ' ';
  char                       qnameSuffixDelim_   = DEFAULT_QNAME_DELIM;

  // testability hook
  template <typename DumpT>
  friend DumpT& dump(DumpT& dump, const BamToReadTransformer& transformer);

public:
  BamToReadTransformer(const char qnameSuffixDelim = DEFAULT_QNAME_DELIM)
    : qnameSuffixDelim_(qnameSuffixDelim)
  {
  }

  void operator()(const bam::BamRecordAccessor& bra, unsigned pos, uint64_t fragmentId, sequences::Read& read)
  {
    const auto name    = bra.getName(qnameSuffixDelim_);
    const auto bases   = bra.getBases();
    const auto qscores = bra.getQscores();
    tmpName_.assign(name.first, name.second);
    tmpBases_.clear();
    if (bra.reverse()) {
      for (auto it = bases.first; bases.second != it; ++it) {
        unsigned char b = *it;
        b               = (b & 0xCC) >> 2 | (b & 0x33) << 2;
        b               = (b & 0xAA) >> 1 | (b & 0x55) << 1;

        tmpBases_.push_back(b >> 4);
        tmpBases_.push_back(b & 0x0f);
      }
    } else {
      for (auto it = bases.first; bases.second != it; ++it) {
        tmpBases_.push_back(*it >> 4);
        tmpBases_.push_back(*it & 0x0f);
      }
    }
    tmpBases_.resize(bra.readLength());
    tmpQscores_.assign(qscores.first, qscores.second);
    if (bra.reverse()) {
      std::reverse(tmpBases_.begin(), tmpBases_.end());
      std::reverse(tmpQscores_.begin(), tmpQscores_.end());
    }
    assert(tmpBases_.size() == tmpQscores_.size());

    read.init(std::move(tmpName_), std::move(tmpBases_), std::move(tmpQscores_), fragmentId, pos);
  }
};

}  // namespace io
}  // namespace dragenos
