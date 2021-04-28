#define SEQUENCES_READ_HPP
#define SEQUENCES_SEED_HPP

#include <cstddef>
#include <cstdint>
#include <ostream>

namespace dragenos {
namespace sequences {

class Read {
public:
  size_t getLength() const { return 151; }
};

class Seed {
public:
  Seed(
      const Read* read,
      unsigned    readPosition,
      unsigned    primaryLength /* , unsigned extensionLength = 0, Direction direction = FORWARD */)
    : read_(read), readPosition_(readPosition), primaryLength_(primaryLength)
  {
  }
  unsigned getPrimaryLength() const { return primaryLength_; }
  unsigned getReadPosition() const { return readPosition_; }
  unsigned getFirstBaseReadPosition(const unsigned halfExtensionLength) const
  {
    return getReadPosition() - halfExtensionLength;
  }
  unsigned getLastBaseReadPosition(const unsigned halfExtensionLength) const
  {
    return getReadPosition() + getPrimaryLength() + halfExtensionLength - 1;
  }
  const Read* getRead() const { return read_; }

  bool operator==(const Seed& that) const
  {
    return read_ == that.read_ && readPosition_ == that.readPosition_ &&
           primaryLength_ == that.primaryLength_;
  }

private:
  const Read* const read_;
  const unsigned    readPosition_;
  const unsigned    primaryLength_;

  friend std::ostream& operator<<(std::ostream& os, const Seed& s)
  {
    return os << "Seed(" << s.readPosition_ << "," << s.primaryLength_ << ")";
  }
};

}  // namespace sequences
}  // namespace dragenos
