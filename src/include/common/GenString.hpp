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

#ifndef COMMON_GEN_STRING_HPP
#define COMMON_GEN_STRING_HPP

#include <sstream>

////////////////////////////////////////////////////////////////////////////////
// Template metaprogram to generate gen_string_() function for each of 1 to 30
// variable-type arguments.  The preprocessor output looks like this:
//
//  template < typename T0 >
//  std::string gen_string_( const T0 &t0) {
//    std::ostringstream oss;
//    oss << t0;
//    return oss.str();
//}
//
//  template < typename T0 , typename T1 >
//  std::string gen_string_( const T0 &t0 , const T1 &t1) {
//    std::ostringstream oss;
//    oss << t0 << t1;
//    return oss.str();
//  }
//
// and so on, up to thirty classes.
#include <boost/preprocessor/repetition/enum_binary_params.hpp>
#include <boost/preprocessor/repetition/enum_params.hpp>
#include <boost/preprocessor/repetition/repeat_from_to.hpp>
#define GEN_STRING_LSHIFTER(z, ITER, data) data##ITER

#define GEN_STRING(z, ITER, data)                                         \
  template <BOOST_PP_ENUM_PARAMS(ITER, typename T)>                       \
  std::string gen_string_(BOOST_PP_ENUM_BINARY_PARAMS(ITER, const T, &t)) \
  {                                                                       \
    std::ostringstream oss;                                               \
    oss                BOOST_PP_REPEAT(ITER, GEN_STRING_LSHIFTER, << t);  \
    return oss.str();                                                     \
  }
BOOST_PP_REPEAT_FROM_TO(1, 30, GEN_STRING, ~)
//
// ... end of boost gen_string_ metaprogram.
////////////////////////////////////////////////////////////////////////////////

#endif
