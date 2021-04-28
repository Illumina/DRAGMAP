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

#ifndef ALIGN_QUERY_HPP
#define ALIGN_QUERY_HPP

#include "sequences/Read.hpp"

namespace dragenos {
namespace align {

/**
 ** \brief Proxy for a read for Smith-Waterman algorithm
 **/
class Query : public sequences::Read::Bases  // public std::vector<unsigned char>
{
public:
  typedef sequences::Read::Bases Bases;
  Query() : beginOffset_(0) {}
  Query(Bases::const_iterator begin, Bases::const_iterator end)
    : sequences::Read::Bases(begin, end), beginOffset_(0)
  {
  }
  template <class I>
  Query(I begin, const I end) : beginOffset_(0)
  {
    while (end != begin) {
      push_back(static_cast<Bases::value_type>(*begin));
      ++begin;
    }
  }
  // conveniently forward all valid costructors from std::vector<unsigned char>
  //template<typename ... Args>
  //Query(Args ... args) : std::vector<char>(std::forward<Args>(args)...), beginOffset_(0) {}
  Query& operator=(const Bases bases)
  {
    beginOffset_   = 0;
    Bases::operator=(bases);
    return *this;
  }
  void setBeginOffset(const size_t offset) { beginOffset_ = offset; }
  void incrementBeginOffset(const size_t offset) { ++beginOffset_; }
  // reimplement the begin method to support an offset
  //auto begin() const -> decltype(begin()) {return begin() + beginOffset_;}
private:
  size_t beginOffset_;

  friend std::ostream& operator<<(std::ostream& os, const Query& query)
  {
    for (const auto& q : query) {
      os << sequences::Read::decodeBase(q);
    }
    return os;
  }
};

}  // namespace align
}  // namespace dragenos

#endif  // #ifndef ALIGN_QUERY_HPP
