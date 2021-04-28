#define COMMON_EXCEPTIONS_HPP
#include <boost/throw_exception.hpp>
#include <stdexcept>
#include <string>

namespace dragenos {
namespace common {
class InvalidParameterException : public std::logic_error {
public:
  InvalidParameterException(const std::string& message) : std::logic_error(message) {}
};

}  // namespace common
}  // namespace dragenos
