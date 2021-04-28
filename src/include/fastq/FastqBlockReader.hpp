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

#include <istream>
#include <vector>

#include <boost/iostreams/stream_buffer.hpp>

namespace dragenos {
namespace fastq {

/**
 * \brief Reads a block of data from fastq input stream and then buffers the last incomplete record so that
 *        parsers don't have to deal with it
 */
class FastqBlockReader {
  std::istream&     stream_;
  std::vector<char> buffer_;

public:
  typedef std::istream::char_type char_type;

  FastqBlockReader(std::istream& stream) : stream_(stream) {}

  /**
   * \brief         read up to n bytes from underlying stream
   * \return        number of bytes in s ending at the end of last complete record
   * \postcondition Subsequent call will result in buffered record returned
   */
  std::size_t read(char_type* s, std::size_t n)
  {
    const std::size_t toCopy = std::min<std::size_t>(buffer_.size(), n);
    std::copy(buffer_.begin(), buffer_.begin() + toCopy, s);
    buffer_.erase(buffer_.begin(), buffer_.begin() + toCopy);
    const std::size_t extra = n - toCopy;
    std::size_t       read  = 0;
    if (extra) {
      if (!stream_.read(s + toCopy, extra) && !stream_.eof()) {
        throw std::logic_error(std::string("Error reading input stream at ") << stream_.tellg());
      }
      read = stream_.gcount();
    }

    std::size_t ret = stream_.eof() ? toCopy + read : findFirstIncomplete(s, toCopy + read);
    buffer_.insert(buffer_.end(), s + ret, s + toCopy + read);
    return ret;
  }

  bool eof() const
  {
    if (!buffer_.empty()) {
      return false;
    }

    return stream_.eof();
  }

private:
  const char_type* findPrevName(const char_type* s, std::size_t n) const
  {
    for (const char_type* p = s + n - 1; n; --n, --p) {
      if ('\n' == *p) {
        --p;
        --n;
        if (4 > n) {
          // not enough chars to skip delimeter, bases and read name
          throw std::logic_error(
              std::string("Unable to parse fastq around position ")
              << stream_.tellg() << " block starts with:\n"
              << std::string(s, s + std::min<int>(n, 300)));
        }
        // handle option Microsoft newline garbage
        if ('\r' == *p) {
          --p;
          --n;
        }
        if ('+' == *p) {
          --p;
          --n;
          if ('\n' == *p) {
            const auto start   = std::reverse_iterator<const char_type*>(s);
            const auto nameEnd = std::find(std::reverse_iterator<const char_type*>(p - 1), start, '\n');
            if (start == nameEnd) {
              throw std::logic_error(
                  std::string("Unable to find read name end in fastq around position ")
                  << stream_.tellg() << " block starts with:\n"
                  << std::string(s, s + std::min<int>(n, 300)));
            }
            const auto nameStart = std::find(nameEnd + 1, start, '\n');
            if (start == nameStart) {
              if ('@' == *start) {
                return s;
              }
            } else if ('@' == *(nameStart - 1)) {
              return nameStart.base();
            }
            throw std::logic_error(
                std::string("Unable to find read name start in fastq around position ")
                << stream_.tellg() << " block starts with:\n"
                << std::string(s, s + std::min<int>(n, 300)));
          }
        }
      }
    }

    throw std::logic_error(
        std::string("Unable to find read name start in fastq around position ")
        << stream_.tellg() << " block starts with:\n"
        << std::string(s, s + std::min<int>(n, 300)));
  }

  bool nameMatches(const char_type* n1, const char_type* n2) const
  {
    static const char delimiters[] = {'#', ' ', '\r', '\n'};

    // 1000 is arbitrary, the caller is responsible to feed the pointers containing at least one full line of
    // the fastq record
    const char_type* n1End = std::find_first_of(n1, n1 + 1000, delimiters, delimiters + sizeof(delimiters));
    const char_type* n2End = std::find_first_of(n2, n2 + 1000, delimiters, delimiters + sizeof(delimiters));

    if (std::distance(n1, n1End) != std::distance(n2, n2End)) {
      return false;
    }

    if (n1End != std::mismatch(n1, n1End, n2).first) {
      return false;
    }

    return true;
  }

  /**
   * \brief find first incomplete record
   * \return number of bytes in s ending at the end of last complete record
   */
  std::size_t findFirstIncomplete(const char_type* s, std::size_t n) const
  {
    const char_type* p = findPrevName(s, n);

    const char_type* pp = findPrevName(s, std::distance(s, p));

    if (nameMatches(p, pp)) {
      return std::distance(s, pp);
    }

    return std::distance(s, p);
  }
};

}  // namespace fastq
}  // namespace dragenos
