#define SEQUENCES_READ_HPP

#include <cstddef>
#include <cstdint>

namespace dragenos {
namespace sequences {

class Read {
public:
  typedef unsigned char Base;
  Base                  getBase4bpb(size_t i) const { return 1 << (getBase2bpb(i)); }
  char                  getBase2bpb(size_t i) const { return ((i * (i + 1)) / 2) % 4; }
  uint64_t              getId() const { return 17; }
  size_t                getLength() const { return 151; }
};

}  // namespace sequences
}  // namespace dragenos
