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
#include <string>

namespace dragenos {
namespace fastq {

/**
 ** \brief Exception thrown when fastq violation is encountered.
 **
 **/
class FastqInvalidFormat : public std::logic_error {
public:
  FastqInvalidFormat(const std::string& message) : std::logic_error("Corrupt fastq: " + message) {}
};

template <typename IT>
class BasicToken {
private:
  bool valid_;
  IT   headerBegin_;
  IT   headerEnd_;
  IT   baseCallsBegin_;
  IT   baseCallsEnd_;
  IT   qScoresBegin_;
  IT   end_;

public:
  BasicToken() : valid_(false) {}

  bool        valid() const { return valid_; }
  IT          end() const { return end_; }
  bool        empty() const { return end_ == headerBegin_; }
  std::size_t readLength() const { return std::distance(baseCallsBegin_, baseCallsEnd_); }

  // @return true if complete record has been read
  bool reset(IT begin, IT end)
  {
    valid_ = false;

    headerBegin_    = findNotNewLine(begin, end);
    headerEnd_      = findNewLine(headerBegin_, end);
    baseCallsBegin_ = findNotNewLine(headerEnd_, end);

    // special case for zero-length reads
    const bool zeroLengthRead = end != baseCallsBegin_ && '+' == *baseCallsBegin_;

    if (zeroLengthRead) {
      baseCallsEnd_ = baseCallsBegin_;
      qScoresBegin_ = findNotNewLine(baseCallsEnd_ + 1, end);
      end_          = qScoresBegin_;
      valid_        = headerBegin_ != headerEnd_;
      return valid_ || end != end_;
    } else {
      baseCallsEnd_ = findNewLine(baseCallsBegin_, end);
      qScoresBegin_ = findNotNewLine(baseCallsEnd_, end);
      if (end != qScoresBegin_) {
        if ('+' != *qScoresBegin_) {
          throw FastqInvalidFormat("'+' is missing in fastq record");
        }
        qScoresBegin_ = findNotNewLine(findNewLine(qScoresBegin_ + 1, end), end);
        end_          = findNewLine(qScoresBegin_, end);
        valid_        = std::distance(qScoresBegin_, end_) == std::distance(baseCallsBegin_, baseCallsEnd_);
        return valid_ || end != end_;
      } else {
        // else ended too soon, stay invalid
        end_ = end;
      }
    }
    return false;
  }

  // skip @ at the start of name
  std::pair<IT, IT> getName(const char qnameSuffixDelim) const
  {
    return make_pair(headerBegin_ + 1, std::find(headerBegin_ + 1, headerEnd_, qnameSuffixDelim));
  }
  std::pair<IT, IT> getBases() const { return std::make_pair(baseCallsBegin_, baseCallsEnd_); }
  std::pair<IT, IT> getQscores() const { return std::make_pair(qScoresBegin_, end_); }

private:
  template <typename IteratorT>
  static IteratorT findNotNewLine(IteratorT itBegin, IteratorT itEnd)
  {
    // this usually ends after first comparison or so
    while (itEnd != itBegin && '\n' == *itBegin) {
      ++itBegin;
    }
    return itBegin;
  }

  template <typename IteratorT>
  static IteratorT findNewLine(IteratorT itBegin, IteratorT itEnd)
  {
    return std::find(itBegin, itEnd, '\n');
  }

  friend std::ostream& operator<<(std::ostream& os, const BasicToken& token)
  {
    return os << "BasicToken(" << std::string(token.headerBegin_, token.headerEnd_) << " "
              << std::string(token.baseCallsBegin_, token.baseCallsEnd_) << " "
              << std::string(token.qScoresBegin_, token.end_)
              << "):" << (token.empty() ? "empty" : token.valid() ? "valid" : "invalid");
  }
};

}  // namespace fastq
}  // namespace dragenos
