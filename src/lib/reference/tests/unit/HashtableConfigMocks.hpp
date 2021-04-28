#include "common/Exceptions.hpp"

namespace dragenos {
namespace common {
ExceptionData::ExceptionData(int errorNumber, const std::string& message)
  : boost::exception(), errorNumber_(errorNumber), message_(message)
{
}

InvalidParameterException::InvalidParameterException(const std::string& message)
  : std::logic_error(message), ExceptionData(EINVAL, message)
{
}

}  // namespace common
}  // namespace dragenos
