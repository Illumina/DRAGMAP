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

#include <assert.h>
#include <cstdint>

namespace dragenos {
namespace bam {

enum Flag : uint16_t {
  MULTIPLE_SEGMENTS               = 0x1,
  ALL_PROPERLY_ALIGNED            = 0x2,
  UNMAPPED                        = 0x4,
  UNMAPPD_NEXT_SEGMENT            = 0x8,
  REVERSE_COMPLEMENT              = 0x10,
  REVERSE_COMPLEMENT_NEXT_SEGMENT = 0x20,
  FIRST_IN_TEMPLATE               = 0x40,
  LAST_IN_TEMPLATE                = 0x80,
  SECONDARY_ALIGNMENT             = 0x100,
  FAILED_FILTERS                  = 0x200,
  DUPLICATE                       = 0x400,
  SUPPLEMENTARY_ALIGNMENT         = 0x800
};

/**
 * \brief Reads a block of data from bam input stream and then buffers the last incomplete record so that
 *        parsers don't have to deal with it
 */
struct BamRecordHeader {
  uint32_t block_size;
  int32_t  refID;
  int32_t  pos;
  uint8_t  l_read_name;
  uint8_t  mapq;
  uint16_t bin;
  uint16_t n_cigar_op;
  uint16_t flag;
  uint32_t l_seq;
  int32_t  next_refID;
  int32_t  next_pos;
  int32_t  tlen;
  char     read_name[];
};

class BamRecordAccessor {
  struct Name : std::pair<const char*, const char*> {
    typedef std::pair<const char*, const char*> Base;
    Name(const char* begin, const char* end) : Base(begin, end) {}

    bool operator==(const Name& that) const
    {
      return std::distance(first, second) == std::distance(that.first, that.second) &&
             std::equal(first, second, that.first);
    }

    friend std::ostream& operator<<(std::ostream& os, const Name& name)
    {
      return os << std::string_view(name.first, std::distance(name.first, name.second));
    }
  };

  const BamRecordHeader* h_;

public:
  BamRecordAccessor(const char* bh = 0) : h_((const BamRecordHeader*)bh) {}

  bool empty() const { return !h_; }
  void envelop(const char* bh) { h_ = (const BamRecordHeader*)bh; }

  std::size_t size() const { return sizeof(h_->block_size) + h_->block_size; }

  const char* current() const { return (const char*)h_; }
  const char* next() const { return (const char*)h_ + size(); }

  bool              secondary() const { return h_->flag & Flag::SECONDARY_ALIGNMENT; }
  bool              suplementary() const { return h_->flag & Flag::SUPPLEMENTARY_ALIGNMENT; }
  bool              reverse() const { return h_->flag & Flag::REVERSE_COMPLEMENT; }
  bool              paired() const { return h_->flag & Flag::MULTIPLE_SEGMENTS; }
  bool              first() const { return h_->flag & Flag::FIRST_IN_TEMPLATE; }
  bool              last() const { return h_->flag & Flag::LAST_IN_TEMPLATE; }
  const std::size_t readLength() const { return h_->l_seq; }

  const char* nameBegin() const { return h_->read_name; }
  const char* nameEnd() const { return nameBegin() + h_->l_read_name; }
  const Name  getName(const char qnameSuffixDelim) const
  {
    assert(nameBegin() != nameEnd());  // unexpected empty name
    const char* end = nameEnd() - 1;
    assert(!*end);  // name includes 0 terminator according to specs
    Name ret(nameBegin(), std::find(nameBegin(), end, qnameSuffixDelim));
    return ret;
  }

  const char* cigarBegin() const { return nameEnd(); }
  const char* cigarEnd() const { return cigarBegin() + h_->n_cigar_op * 4; }

  const unsigned char* basesBegin() const { return (const unsigned char*)cigarEnd(); }
  const unsigned char* basesEnd() const { return basesBegin() + (h_->l_seq + 1) / 2; }
  const std::pair<const unsigned char*, const unsigned char*> getBases() const
  {
    return std::make_pair(basesBegin(), basesEnd());
  }

  const unsigned char* qscoresBegin() const { return basesEnd(); }
  const unsigned char* qscoresEnd() const { return qscoresBegin() + h_->l_seq; }
  const std::pair<const unsigned char*, const unsigned char*> getQscores() const
  {
    return std::make_pair(qscoresBegin(), qscoresEnd());
  }

  friend std::ostream& operator<<(std::ostream& os, const BamRecordAccessor& bra)
  {
    return os << "BamRecordAccessor(" << bra.size() << "s " << std::string(bra.nameBegin(), bra.nameEnd())
              << ")";
  }
};

}  // namespace bam
}  // namespace dragenos
