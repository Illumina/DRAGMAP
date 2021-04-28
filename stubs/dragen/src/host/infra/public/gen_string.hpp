// Copyright 2013 Edico Genome Corporation. All rights reserved.
//
// This file contains confidential and proprietary information of the Edico
// Genome Corporation and is protected under the U.S. and international
// copyright and other intellectual property laws.
//
// $Id$
// $Author$
// $Change$
// $DateTime$
//
#ifndef __GEN_STRING_HPP__
#define __GEN_STRING_HPP__

#include <sstream>

/// \brief Convert variable arguments into a concatenated string, based on each
/// argument's operator<< function
///
template <typename... Args> std::string gen_string_(Args &&... args) {
  // The dummy braced-init-list type for variadic parameter pack expansion below
  // uses an empty struct to reduce generated code size
  struct empty {};
  std::ostringstream oss;
  (void)(empty[]){(oss << args, empty{})...};
  return oss.str();
}

#endif
