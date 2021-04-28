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

#ifndef COMMON_EXCEPTIONS_HPP
#define COMMON_EXCEPTIONS_HPP

#include <boost/cerrno.hpp>
#include <boost/exception/all.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/throw_exception.hpp>
#include <ios>
#include <stdexcept>
#include <string>

namespace dragenos {
namespace common {

/**
 ** \brief Virtual base class to all the exception classes in iSAAC.
 **
 ** Use BOOST_THROW_EXCEPTION to get the context info (file, function, line)
 ** at the throw site.
 **/
class ExceptionData : public boost::exception {
public:
  ExceptionData(int errorNumber = 0, const std::string& message = "");
  ExceptionData(const ExceptionData& e)
    : boost::exception(e), errorNumber_(e.errorNumber_), message_(e.message_)
  {
  }
  virtual ~ExceptionData() throw() {}
  int                getErrorNumber() const { return errorNumber_; }
  const std::string& getMessage() const { return message_; }
  std::string        getContext() const;

private:
  const int         errorNumber_;
  const std::string message_;
  ExceptionData&    operator=(const ExceptionData&);
};

class DragenOsException : public std::exception, public ExceptionData {
public:
  DragenOsException(int errorNumber, const std::string& message) : ExceptionData(errorNumber, message) {}
  DragenOsException(const std::string& message) : ExceptionData(0, message) {}
  DragenOsException(const DragenOsException& e) : std::exception(e), ExceptionData(e) {}
  virtual const char* what() const throw() { return getMessage().c_str(); }

private:
  DragenOsException& operator=(const DragenOsException&);
};

/**
 * \brief Exception thrown when there are problems with the IO operations
 */
class IoException : public std::ios_base::failure, public ExceptionData {
public:
  IoException(int errorNumber, const std::string& message);
};

/**
 * \brief Exception thrown when there is insufficient resources to perform an operation. For example
 *        if the adjusting the soft ulimit fails due to a set hard limit
 */
class ResourceException : public std::exception, public ExceptionData {
public:
  ResourceException(int errorNumber, const std::string& message);
};

/**
 * \brief Same as bad_alloc but with a message
 */
class MemoryException : public std::bad_alloc, public ExceptionData {
public:
  MemoryException(const std::string& message);
};

/**
 ** \brief Exception thrown when the client supplied and unsupported version number.
 **
 ** Particularly relevant to data format and software versions
 ** (Pipeline, IPAR, Phoenix, etc.). It should not be used in
 ** situations where the client didn't have the possibility to check
 ** the version (for instance when reading the version of a data
 ** format from the header of a file).
 **
 **/
class UnsupportedVersionException : public std::logic_error, public ExceptionData {
public:
  UnsupportedVersionException(const std::string& message);
};

/**
 ** \brief Thrown when the requested functionality is not available.
 **
 **/
class FeatureNotAvailable : public std::logic_error, public ExceptionData {
public:
  FeatureNotAvailable(const std::string& message);
};

/**
 ** \brief Exception thrown when the client supplied an invalid parameter.
 **
 **/
class InvalidParameterException : public std::logic_error, public ExceptionData {
public:
  InvalidParameterException(const std::string& message);
};

/**
 ** \brief Exception thrown when an invalid command line option was detected.
 **
 **/
class InvalidOptionException : public std::logic_error, public ExceptionData {
public:
  InvalidOptionException(const std::string& message);
};

/**
 ** \brief Exception thrown when a method invocation violates the pre-conditions.
 **
 **/
class PreConditionException : public std::logic_error, public ExceptionData {
public:
  PreConditionException(const std::string& message);
};

/**
 ** \brief Exception thrown when a method invocation violates the post-conditions.
 **
 **/
class PostConditionException : public std::logic_error, public ExceptionData {
public:
  PostConditionException(const std::string& message);
};

/**
 ** \brief Exception thrown when a libxslt method invocation fails.
 **
 **/
class LibXsltException : public DragenOsException {
public:
  LibXsltException();
};

}  // namespace common
}  // namespace dragenos

#endif  // #ifndef COMMON_EXCEPTIONS_HPP
