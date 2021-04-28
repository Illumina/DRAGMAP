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

#ifndef REFERENCE_REFERENCE_DIR_HPP
#define REFERENCE_REFERENCE_DIR_HPP

#include <boost/filesystem.hpp>
#include <functional>
#include <memory>
#include <vector>

#include "reference/HashtableConfig.hpp"
#include "reference/ReferenceSequence.hpp"

namespace dragenos {
namespace reference {

class ReferenceDir {
public:
  virtual ~ReferenceDir();
  virtual const reference::HashtableConfig& getHashtableConfig() const   = 0;
  virtual const uint64_t*                   getHashtableData() const     = 0;
  virtual const uint64_t*                   getExtendTableData() const   = 0;
  virtual const ReferenceSequence&          getReferenceSequence() const = 0;
};

class ReferenceDir7 : public ReferenceDir {
  const bool mmap_;
  const bool load_;

public:
  ReferenceDir7(const boost::filesystem::path& path, bool mmap, bool load);
  ~ReferenceDir7();
  virtual const reference::HashtableConfig& getHashtableConfig() const { return hashtableConfig_; };
  virtual const uint64_t*                   getHashtableData() const { return hashtableData_.get(); }
  virtual const uint64_t*                   getExtendTableData() const { return extendTableData_.get(); }
  virtual const ReferenceSequence&          getReferenceSequence() const { return *referenceSequencePtr_; }
  size_t                                    getHashtableConfigSize() const;
  size_t                                    getHashtableDataSize() const;

protected:
  static constexpr auto         hashtableConfigBin = "hash_table.cfg.bin";
  static constexpr auto         hashtableBin       = "hash_table.bin";
  static constexpr auto         extendTableBin     = "extend_table.bin";
  static constexpr auto         referenceBin       = "reference.bin";
  static constexpr auto         hashTableCmp       = "hash_table.cmp";
  const boost::filesystem::path path_;
  // raw binary content of the hashtable config file
  const std::vector<char> hashtableConfigData_;
  // TODO: replace with a placement new
  const reference::HashtableConfig hashtableConfig_;
  /**
   ** \brief memory mapped hashtable data
   **
   ** Note: the custom destructor would be munmap, but munmap needs to know
   ** the size of the memory segment (hashtableConfig_.getHashtableBytes())
   ** which requires using a lambda as a constructor, which in turn requires
   ** declaring the type of the destructor as "std::function<void(void*)>>"
   ** instead of using the type of the function pointer "void(*)(uint64_t*)"
   **
   ** Note: the type can't be const because of munmap signature.
   **/
  typedef std::unique_ptr<uint64_t, std::function<void(uint64_t*)>> Uint64Ptr;
  Uint64Ptr                                                         hashtableData_;
  /**
   ** \brief memory mapped hashtable data
   **
   ** Note: using mmap and munmap like hashtableData_
   **/
  Uint64Ptr         extendTableData_;
  std::vector<char> getHashtableConfigData() const;
  template <typename T>
  std::unique_ptr<T, std::function<void(T*)>> mmapData(
      std::string binFile, size_t expectedBinFileBytes) const;
  template <typename T>
  std::unique_ptr<T, std::function<void(T*)>> readData(
      const std::string binFile, const size_t expectedBinFileBytes) const;
  typedef std::unique_ptr<unsigned char, std::function<void(unsigned char*)>> UcharPtr;
  UcharPtr                                                                    referenceData_;
  std::unique_ptr<ReferenceSequence>                                          referenceSequencePtr_;

  UcharPtr ReadFileIntoBuffer(const boost::filesystem::path& directory, std::streamsize& size);
};

}  // namespace reference
}  // namespace dragenos

#endif  // #ifndef REFERENCE_REFERENCE_DIR_HPP
