// Copyright 2016 Edico Genome Corporation. All rights reserved.
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

#pragma once

#include <assert.h>
#include <string.h>

#include <iostream>

#include "gen_string.hpp"

template <typename... Args>
void __attribute__((noinline))
dragen_assert_fail_(const char *file, const int line, const char *pred,
                    Args &&... args) {
  struct empty {};
  (void)(empty[]){(std::cerr << "Assertion failed in " << file << " line "
                             << line << " -- " << pred << " -- " << args,
                   empty{})...};
  std::cerr << std::endl;
  assert(false);
}

#define ASSERT(pred, ...)                                                      \
  if (__builtin_expect(!(bool)(pred), false)) {                                \
    dragen_assert_fail_(__FILE__, __LINE__, #pred, ##__VA_ARGS__);             \
  }

/////////////////////////////////////////////////////////////////////////////////////
// Base class for non-fatal exceptions, which are not handled in the same way as
// fatal exceptions.  For example, they are logged differently, do not generate
// pstacks, and do not call dragen_reset
//
class NonFatalException : public std::exception {
public:
  NonFatalException(std::string message);
  virtual ~NonFatalException() throw();

  const char *what() const throw();

protected:
  std::string m_className;
  std::string m_message;
};

template <class T> void dragen_throw(const std::string &msg) { throw T(msg); }

/////////////////////////////////////////////////////////////////////////////
// Use this macro with non-fatal exceptions, which are caught and logged
// differently than fatal exceptions.
#define THROW_NONFATAL(classname, ...)                                         \
  {                                                                            \
    std::string msg = gen_string_(__VA_ARGS__);                                \
    dragen_throw<classname>(msg);                                              \
  }

/////////////////////////////////////////////////////////////////////////////
// Define additional exception classes with this macro.  You specify the new
// exception classname, and the name of the parent class
//
#define DEFINE_EXCEPTION(exception, parent)                                    \
  class exception : public parent {                                            \
  public:                                                                      \
    exception(std::string msg) : parent(msg) { m_className = #exception; }     \
  };
