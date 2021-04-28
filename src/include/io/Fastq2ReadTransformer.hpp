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

#include "fastq/Tokenizer.hpp"
#include "sequences/Read.hpp"

namespace dragenos {
namespace io {

#define VECTOR_REGISTER_WIDTH 16
//#define VECTOR_REGISTER_WIDTH 0

class FastqToReadTransformer {
public:
  static const char DEFAULT_Q0          = 33;
  static const char DEFAULT_QNAME_DELIM = ' ';
  char              qnameSuffixDelim_   = DEFAULT_QNAME_DELIM;

  sequences::Read::Name      tmpName_;
  sequences::Read::Bases     tmpBases_;
  sequences::Read::Qualities tmpQscores_;

  sequences::Read::Qscore q0_[VECTOR_REGISTER_WIDTH ? VECTOR_REGISTER_WIDTH : 1];
  sequences::Read::Qscore q2_[VECTOR_REGISTER_WIDTH ? VECTOR_REGISTER_WIDTH : 1];

  // testability hook
  template <typename DumpT>
  friend DumpT& dump(DumpT& dump, const FastqToReadTransformer& transformer);

public:
  FastqToReadTransformer(const char qnameSuffixDelim = DEFAULT_QNAME_DELIM, const char q0 = DEFAULT_Q0)
    : qnameSuffixDelim_(qnameSuffixDelim)
  {
    std::fill(q0_, q0_ + (VECTOR_REGISTER_WIDTH ? VECTOR_REGISTER_WIDTH : 1), q0);
    std::fill(q2_, q2_ + (VECTOR_REGISTER_WIDTH ? VECTOR_REGISTER_WIDTH : 1), 2);
  }

  void operator()(
      const fastq::Tokenizer::Token& fastqToken, unsigned pos, uint64_t fragmentId, sequences::Read& read)
  {
    const auto name    = fastqToken.getName(qnameSuffixDelim_);
    const auto bases   = fastqToken.getBases();
    const auto qscores = fastqToken.getQscores();
    tmpName_.assign(name.first, name.second);
    tmpBases_.assign(bases.first, bases.second);
    tmpQscores_.assign(qscores.first, qscores.second);

#if VECTOR_REGISTER_WIDTH
    const std::size_t minCapacity =
        ((tmpBases_.size() + VECTOR_REGISTER_WIDTH - 1) / VECTOR_REGISTER_WIDTH) * VECTOR_REGISTER_WIDTH;
    const std::size_t size = tmpBases_.size();
    tmpBases_.resize(minCapacity);
    tmpBases_.resize(size);
    tmpQscores_.resize(minCapacity);
    tmpQscores_.resize(size);
    convertBasesAndQualities();
#else
    convertBasesAndQualities(q0_[0]);
#endif  // VECTOR_REGISTER_WIDTH

    read.init(std::move(tmpName_), std::move(tmpBases_), std::move(tmpQscores_), fragmentId, pos);
  }

private:
  void convertBasesAndQualities();
  void convertBasesAndQualities(const char q0);

  void convertBases(sequences::Read::Base* b, const std::size_t size);
  void convertBases2(sequences::Read::Base* b, const std::size_t size);
  void convertQualities(const sequences::Read::Base* b, sequences::Read::Qscore* q, const std::size_t size);
};

}  // namespace io
}  // namespace dragenos
