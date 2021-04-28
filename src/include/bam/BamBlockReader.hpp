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

#include "bam/Bam.hpp"
#include "common/Debug.hpp"
#include "common/Exceptions.hpp"

namespace dragenos {
namespace bam {

/**
 * \brief Same as logic_error
 */
class BamException : public std::logic_error {
public:
  BamException(const std::string& message) : std::logic_error(message) {}
};

/**
 * \brief Reads a block of data from bam input stream and then buffers the last incomplete record so that
 *        parsers don't have to deal with it
 */
class BamBlockReader {
  std::istream&     stream_;
  std::vector<char> buffer_;
  char              inputQnameSuffixDelim_;

public:
  typedef std::istream::char_type char_type;

  BamBlockReader(std::istream& stream, const char inputQnameSuffixDelim) : stream_(stream)
  {
    skipToFirstRecord();
  }

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
      try {
        if (!stream_.read(s + toCopy, extra) && !stream_.eof()) {
          throw std::logic_error(std::string("Error reading input stream at ") << stream_.tellg());
        }
      } catch (boost::iostreams::gzip_error& e) {
        // it is ok to have CRC error as bam in particular sets it to 0
        if (boost::iostreams::gzip::bad_crc != e.error()) {
          throw;
        }
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
  /**
   * \brief find first incomplete record
   * \return number of bytes in s ending at the end of last complete record
   */
  std::size_t findFirstIncomplete(const char_type* s, std::size_t n) const
  {
    const char_type*  prevPrev = 0;
    const char_type*  prev     = 0;
    BamRecordAccessor bra(s);
    while (sizeof(BamRecordHeader) < n && bra.size() <= n) {
      n -= bra.size();
      if (!bra.secondary() && !bra.suplementary()) {
        prevPrev = prev;
        prev     = bra.current();
      }
      //      std::cerr << bra << std::endl;
      bra.envelop(bra.next());
    }

    if (!prev) {
      BOOST_THROW_EXCEPTION(BamException("Insufficient buffer to read a single bam record"));
    }

    BamRecordAccessor prevBra(prev);
    if (prevBra.paired()) {
      if (!prevPrev || BamRecordAccessor(prevPrev).getName(inputQnameSuffixDelim_) ==
                           prevBra.getName(inputQnameSuffixDelim_)) {
        return std::distance(s, (const char*)prevBra.next());
      }
    }
    return std::distance(s, (const char*)prev);
  }

  void readMagic()
  {
    char magic[4];
    if (!stream_.read(magic, sizeof(magic))) {
      BOOST_THROW_EXCEPTION(BamException("Unable to read magic from bam stream"));
    }

    static const char expected[sizeof(magic)] = {'B', 'A', 'M', 1};
    if (!std::equal(magic, magic + sizeof(magic), expected)) {
      BOOST_THROW_EXCEPTION(BamException("Incorrect magic read from bam stream"));
    }
  }

  void readName()
  {
    uint32_t l_text = 0;
    if (!stream_.read((char*)&l_text, sizeof(l_text))) {
      BOOST_THROW_EXCEPTION(BamException("Unable to read l_text from bam stream"));
    }

    char text[l_text];
    if (!stream_.read(text, l_text)) {
      BOOST_THROW_EXCEPTION(
          BamException(std::string("Unable to read text from bam stream. l_text") << l_text));
    }
  }

  void readRef()
  {
    uint32_t l_name = 0;
    if (!stream_.read((char*)&l_name, sizeof(l_name))) {
      BOOST_THROW_EXCEPTION(BamException("Unable to read l_text from bam stream"));
    }

    char name[l_name];
    if (!stream_.read(name, l_name)) {
      BOOST_THROW_EXCEPTION(
          BamException(std::string("Unable to read text from bam stream. l_name") << l_name));
    }

    uint32_t l_ref = 0;
    if (!stream_.read((char*)&l_ref, sizeof(l_ref))) {
      BOOST_THROW_EXCEPTION(BamException("Unable to read l_ref from bam stream"));
    }
  }

  void readReferences()
  {
    uint32_t n_ref = 0;
    if (!stream_.read((char*)&n_ref, sizeof(n_ref))) {
      BOOST_THROW_EXCEPTION(BamException("Unable to read n_ref from bam stream"));
    }

    while (n_ref--) {
      readRef();
    }
  }

  void skipToFirstRecord()
  {
    readMagic();
    readName();
    readReferences();
  }
};

}  // namespace bam
}  // namespace dragenos
