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

#include <istream>
#include <vector>

#include <boost/iostreams/operations.hpp>

#include "common/Debug.hpp"

namespace dragenos {
namespace fastq {

class FastqNRecordReader {
  std::istream& stream_;

public:
  FastqNRecordReader(std::istream& stream) : stream_(stream) {}

  /**
   * \brief         read up to n fastq records from the underlying stream
   * \return        number of records read
   * \postcondition Subsequent call will result in buffered record returned
   */
  template <typename InsertIt>
  std::size_t read(InsertIt it, std::size_t n)
  {
    static const int FASTQ_LINES_PER_RECORD = 4;
    int              lines                  = n * FASTQ_LINES_PER_RECORD;

    while (lines && readLine(it)) {
      --lines;
    }

    return n - lines / FASTQ_LINES_PER_RECORD;
  }

  bool eof() const { return stream_.eof(); }

private:
  template <typename InsertIt>
  bool readLine(InsertIt it)
  {
    int c = 0;
    while (EOF != (c = boost::iostreams::get(stream_))) {
      assert(boost::iostreams::WOULD_BLOCK != c);
      if ('\r' == c) {
        continue;
      }

      *(it++) = c;
      if ('\n' == c) {
        break;
      }
    }

    if (!stream_ && !stream_.eof()) {
      throw std::logic_error(
          std::string("Error '") << strerror(errno) << "' reading input stream at " << stream_.tellg());
    }

    return !stream_.eof();
  }
};

}  // namespace fastq
}  // namespace dragenos
