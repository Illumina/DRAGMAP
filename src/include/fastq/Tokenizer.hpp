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
 ** \file fastq/Token.hh
 **
 ** Component to read FASTQ files.
 **
 ** \author Roman Petrovski
 **/

#pragma once

#include <algorithm>
#include <exception>
#include <iostream>
#include <string>
#include <vector>

#include "common/Debug.hpp"
#include "fastq/Token.hpp"

namespace dragenos {
namespace fastq {

class Tokenizer {
  typedef std::vector<char> BufferType;

public:
  typedef BasicToken<BufferType::iterator> Token;

private:
  const bool               mixedNewline_ = false;
  std::istream&            input_;
  static const std::size_t DEFAULT_BUFFER_SIZE_ = 1024 * 1024;
  BufferType               buffer_;
  BufferType::iterator     bufferIterator_;
  Token                    currentToken_;

public:
  Tokenizer(
      std::istream&     input,
      const std::size_t bufferSize   = DEFAULT_BUFFER_SIZE_,
      const bool        mixedNewline = false)
    : mixedNewline_(mixedNewline), input_(input)
  {
    buffer_.reserve(bufferSize);
    bufferIterator_ = buffer_.begin();
  }
  const Token& token() const { return currentToken_; }

  bool next();
};

}  // namespace fastq
}  // namespace dragenos
