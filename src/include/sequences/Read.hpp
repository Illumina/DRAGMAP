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

#ifndef SEQUENCES_READ_HPP
#define SEQUENCES_READ_HPP

#include <iostream>
#include <string>
#include <vector>

namespace dragenos {
namespace sequences {

class Read {
public:
  typedef std::vector<char>   Name;
  typedef unsigned char       Base;
  typedef std::vector<Base>   Bases;
  typedef unsigned char       Qscore;
  typedef std::vector<Qscore> Qualities;

  explicit Read() : id_(0), position_(0) {}

  // The expected usage pattern is to create a bunch of persistent Read objects
  // and init them for each new input sequence from the corresponding parser.
  // No temporaries, no copying
  Read(Read&& that)      = delete;
  Read(const Read& that) = delete;
  // move is ok
  Read& operator=(Read&& that);

  void init(Name&& name, Bases&& bases, Qualities&& qualities, uint64_t id, unsigned position);

  uint64_t          getId() const { return id_; }
  unsigned          getPosition() const { return position_; }
  const Name&       getName() const { return name_; }
  const std::string getNameAsString() const { return std::string(name_.begin(), name_.end()); }
  Base              getBase4bpb(size_t i) const { return bases_[i]; }
  Base              getBase2bpb(size_t i) const
  {
    constexpr char    A           = 0;
    constexpr char    C           = 1;
    constexpr char    G           = 2;
    constexpr char    T           = 3;
    constexpr char    N           = G;
    static const char bases2bpb[] = {N, A, C, N, G, N, N, N, T, N, N, N, N, N, N, N};
    return bases2bpb[getBase4bpb(i)];
  }
  const Bases&     getBases() const { return bases_; }
  const Bases&     getRcBases() const { return rcBases_; }
  const Qualities& getQualities() const { return qualities_; }
  static char      decodeBase(unsigned int b)
  {
    static const char bases[] = {
        'N', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'};
    return b >= sizeof(bases) ? 'N' : bases[b];
  }
  static char decodeRcBase(unsigned int b)
  {
    static const char bases[] = {
        'N', 'T', 'G', 'K', 'C', 'Y', 'S', 'B', 'A', 'W', 'R', 'D', 'M', 'H', 'V', 'N'};
    return b >= sizeof(bases) ? 'N' : bases[b];
  }
  size_t               getLength() const { return bases_.size(); }
  friend std::ostream& operator<<(std::ostream& os, const Read& read)
  {
    return os << "Read(" << read.id_ << "," << read.position_ << ","
              << std::string(read.name_.begin(), read.name_.end()) << "," << read.getLength() << ","
              << (read.bases_.empty() ? std::string("empty") : decodeBase(read.bases_[0]) + std::string(".."))
              << ","
              << (read.qualities_.empty() ? std::string("empty") : std::to_string(read.qualities_[0]) + "..")
              << ")";
  }

private:
  uint64_t  id_;
  unsigned  position_;
  Name      name_;
  Bases     bases_;
  Bases     rcBases_;
  Qualities qualities_;
};

class SerializedRead {
  short nameLen_;
  short readLen_;
  char  bytes_[];

public:
  void operator<<(const Read& read)
  {
    nameLen_        = read.getName().size();
    readLen_        = read.getLength();
    char* bases     = std::copy(read.getName().begin(), read.getName().end(), bytes_);
    char* qualities = std::copy(read.getBases().begin(), read.getBases().end(), bases);
    std::copy(read.getQualities().begin(), read.getQualities().end(), qualities);
  }

  std::size_t getByteSize() const { return sizeof(SerializedRead) + nameLen_ + readLen_ * 2; }

  static std::size_t getByteSize(const Read& read)
  {
    return sizeof(SerializedRead) + read.getName().size() + read.getLength() * 2;
  }

  class SerializationString {
    const char* begin_;
    const char* end_;

  public:
    SerializationString(const char* begin, const char* end) : begin_(begin), end_(end) {}

    const char* begin() const { return begin_; }
    const char* end() const { return end_; }

    friend std::ostream& operator<<(std::ostream& os, const SerializationString& str)
    {
      return os << std::string(str.begin_, str.end_);
    }
  };

  SerializationString getName() const { return SerializationString(bytes_, bytes_ + nameLen_); }
  std::pair<const Read::Base*, const Read::Base*> getBases() const
  {
    return std::make_pair(
        reinterpret_cast<const Read::Base*>(bytes_) + nameLen_,
        reinterpret_cast<const Read::Base*>(bytes_) + nameLen_ + readLen_);
  }
  std::pair<const Read::Qscore*, const Read::Qscore*> getQualities() const
  {
    return std::make_pair(
        reinterpret_cast<const Read::Qscore*>(bytes_) + nameLen_ + readLen_,
        reinterpret_cast<const Read::Qscore*>(bytes_) + nameLen_ + readLen_ * 2);
  }

  static char decodeBase(unsigned int b) { return Read::decodeBase(b); }
  static char decodeRcBase(unsigned int b) { return Read::decodeRcBase(b); }
};

}  // namespace sequences
}  // namespace dragenos

#endif  // #ifndef SEQUENCES_READ_HPP
