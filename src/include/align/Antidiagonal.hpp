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

#ifndef ALIGN_ANTIDIAGONAL_HPP
#define ALIGN_ANTIDIAGONAL_HPP

#include <array>
#include <cstddef>
#include <iterator>
#include <limits>
#include <ostream>

namespace dragenos {
namespace align {

/**
 ** \brief encapsulates implementation detail for wavefront data
 **
 ** The important details include:
 ** - the type of data (assuming that short will do for short reads)
 ** - the number of cells (48 as in the existing HW implementation)
 ** - the alignment of the data (16 for integer SIMD instructions)
 **
 ** Note that this could have been a simple typedef except that alignas can't
 ** be specified in a typedef and there isn't any easy way to have an array
 ** with each of the elements properly aligned, which is needed because the
 ** wavefront needs to keep track of the short term history (previous E and F
 ** component and two previous H components).
 **/
template <typename T, int WIDTH, int ALIGN>
struct AntidiagonalT {
  static const T                             IMPOSSIBLE_SCORE_ = std::numeric_limits<T>::min();
  static constexpr size_t                    width             = WIDTH;
  typedef T                                  value_type;
  typedef std::array<value_type, width>      data_type;
  typedef typename data_type::const_iterator const_iterator;
  AntidiagonalT() : data() {}
  template <typename... E>
  AntidiagonalT(E&&... e) : data{{static_cast<value_type>(e)...}}
  {
  }  // explicit narrowing of int
  //AntidiagonalT(E&&...e) : data{{std::forward<E>(e)...}} {} // this would generate lots of narrowing warnings
  template <typename I>
  AntidiagonalT(const I begin, const I end)
    : data([begin, end]() {
        data_type tmp;
        for (size_t i = 0; (tmp.size() > i) && (begin + i != end); ++i) {
          tmp[i] = *(begin + i);
        };
        return tmp;
      }())
  {
  }
  alignas(ALIGN) data_type data;
  // expose relevant methods from the data for convenience of use
  auto size() const -> decltype(data.size()) { return data.size(); }
  auto operator[](size_t pos) -> decltype(data[pos]) { return data[pos]; }
  auto operator[](size_t pos) const -> decltype(data[pos]) { return data[pos]; }
  auto begin() -> decltype(data.begin()) { return data.begin(); }
  auto begin() const -> decltype(data.begin()) { return data.begin(); }
  auto cbegin() const -> decltype(data.cbegin()) { return data.cbegin(); }
  auto rbegin() const -> decltype(data.rbegin()) { return data.rbegin(); }
  auto end() -> decltype(data.end()) { return data.end(); }
  auto end() const -> decltype(data.end()) { return data.end(); }
  auto cend() const -> decltype(data.cend()) { return data.cend(); }
  auto front() -> decltype(data.front()) { return data.front(); }
  auto front() const -> decltype(data.front()) { return data.front(); }
  auto back() -> decltype(data.back()) { return data.back(); }
  auto back() const -> decltype(data.back()) { return data.back(); }
  T    at(int pos) const { return pos < 0 || int(data.size()) <= pos ? 0 : data.at(pos); }

  friend std::ostream& operator<<(std::ostream& os, const AntidiagonalT& ad)
  {
    // print in reverse to match the order of sim map_sw_scores.log
    for (auto it = ad.data.rbegin(); ad.data.rend() != it; ++it) {
      os << *it << ",";
    }
    return os;
  }
};

// for testing
typedef AntidiagonalT<short, 48, 16> Antidiagonal48;

}  // namespace align
}  // namespace dragenos

#endif  // #ifndef ALIGN_ANTIDIAGONAL_HPP
