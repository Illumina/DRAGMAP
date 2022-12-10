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

#include <fcntl.h>
#include <sys/mman.h>
#include <fstream>
#include <iostream>
#include <thread>

#include <boost/format.hpp>

#include "common/hash_generation/hash_table_compress.h"
#include "reference/ReferenceDir.hpp"

namespace dragenos {
namespace reference {

ReferenceDir::~ReferenceDir() {}

//-------------------------------------------------------------------------------swhitmore
// ReadFileIntoBuffer - Read givin .bin file into buffer.  Memory is allocated in this
// method.  Returns a pointer to the buffer and #size# of the buffer.  Throws an
// exception on error.
//
typename ReferenceDir7::UcharPtr ReferenceDir7::ReadFileIntoBuffer(
    const boost::filesystem::path& filePath, std::streamsize& size)
{
  std::string path = filePath.string();
  if (!boost::filesystem::exists(path)) {
    //    THROW(DragenException, "Could not load reference - ", path, " does not exist");
    throw std::logic_error("Could not load reference - " + path + " does not exist");
  }

  std::ifstream file(path, std::ios::binary | std::ios::ate);

  if (!file) {
    //    THROW(DragenException, "Could not load reference - could not open ", path);
    throw std::logic_error("Could not load reference - could not open" + path);
  }

  size = file.tellg();
  file.seekg(0, file.beg);

  UcharPtr bufPtr(new uint8_t[size], [](uint8_t* p) -> void { delete [](p); });
  file.read(reinterpret_cast<char*>(bufPtr.get()), size);
  if (!file) {
    //    THROW(DragenException, "Could not load reference - could not read ", path);
    throw std::logic_error("Could not load reference - could not read" + path);
  }
  file.close();
  return bufPtr;
}

ReferenceDir7::ReferenceDir7(const boost::filesystem::path& path, bool mmap, bool load)
  : path_(path),
    hashtableConfigData_(getHashtableConfigData()),
    hashtableConfig_(hashtableConfigData_.data(), hashtableConfigData_.size())
{
  if (mmap) {
    hashtableData_   = mmapData<uint64_t>(hashtableBin, hashtableConfig_.getHashtableBytes());
    extendTableData_ = (exists(path_ / extendTableBin))
                           ? mmapData<uint64_t>(extendTableBin, hashtableConfig_.getExtendTableBytes())
                           : nullptr;
    referenceData_ = mmapData<unsigned char>(referenceBin, hashtableConfig_.getReferenceSequenceLength() / 2);
  } else if (load) {
    hashtableData_   = readData<uint64_t>(hashtableBin, hashtableConfig_.getHashtableBytes());
    extendTableData_ = (exists(path_ / extendTableBin))
                           ? readData<uint64_t>(extendTableBin, hashtableConfig_.getExtendTableBytes())
                           : nullptr;
    referenceData_ = readData<unsigned char>(referenceBin, hashtableConfig_.getReferenceSequenceLength() / 2);
  } else  // uncompress
  {
    std::streamsize hashcmpsize = 0, refsize = 0;
    referenceData_         = ReadFileIntoBuffer(path_ / referenceBin, refsize);
    UcharPtr hashcmpbufPtr = ReadFileIntoBuffer(path_ / hashTableCmp, hashcmpsize);

    uint64_t hashsize        = hashtableConfig_.getHashtableBytes();
    uint64_t extendTableSize = hashtableConfig_.getExtendTableBytes();

    // I don't know why decompHashTable needs pointer to pointer to hashbuf and extendTableBuf RP.
    uint8_t* hashbuf        = (uint8_t*)malloc(hashsize);
    uint8_t* extendTableBuf = (uint8_t*)malloc(extendTableSize);

    const int numThreads = std::thread::hardware_concurrency();

    // dragen likes to log to stdout. Redirect to stderr. Affects both c and c++ code
    int stdoutori = dup(1);
    dup2(2, 1);

    char* err = decompHashTable(
        numThreads,
        hashcmpbufPtr.get(),
        hashcmpsize,
        referenceData_.get(),
        refsize,
        &hashbuf,
        &hashsize,
        &extendTableBuf,
        &extendTableSize,
        nullptr,
        nullptr);

    if (err) {
      BOOST_THROW_EXCEPTION(std::logic_error(err));
    }

    // restore stdout
    dup2(stdoutori, 1);

    hashtableData_ =
        Uint64Ptr(reinterpret_cast<uint64_t*>(hashbuf), [](uint64_t* p) -> void { free(p); });
    extendTableData_ =
        Uint64Ptr(reinterpret_cast<uint64_t*>(extendTableBuf), [](uint64_t* p) -> void { free(p); });
  }

  referenceSequencePtr_ = std::unique_ptr<ReferenceSequence>(new ReferenceSequence(
      hashtableConfig_.getTrimmedRegions(),
      referenceData_.get(),
      hashtableConfig_.getReferenceSequenceLength() / 2));
}

ReferenceDir7::~ReferenceDir7() {}

static void checkDirectoryAndFile(const boost::filesystem::path& dir, const boost::filesystem::path& file)
{
  using namespace dragenos::common;
  if (!exists(dir))
    BOOST_THROW_EXCEPTION(
        IoException(ENOENT, std::string("ERROR: directory ") + dir.string() + " doesn't exist"));
  const auto filePath = dir / file;
  if (!exists(filePath))
    BOOST_THROW_EXCEPTION(IoException(ENOENT, std::string("ERROR: file not found: ") + filePath.string()));
  if (!is_regular_file(filePath))
    BOOST_THROW_EXCEPTION(IoException(ENOENT, std::string("ERROR: not a file: ") + filePath.string()));
}

static uintmax_t getFileSize(const boost::filesystem::path& filePath)
{
  using namespace dragenos::common;
  boost::system::error_code ec;
  const auto                fileSize = file_size(filePath, ec);
  if (ec)
    BOOST_THROW_EXCEPTION(
        IoException(ec.value(), std::string("ERROR: failed to get stats for file: ") + filePath.string()));
  return fileSize;
}

std::vector<char> ReferenceDir7::getHashtableConfigData() const
{
  using namespace dragenos::common;
  checkDirectoryAndFile(path_, hashtableConfigBin);
  const auto        hashtableConfigFile = path_ / hashtableConfigBin;
  const auto        fileSize            = getFileSize(hashtableConfigFile);
  std::vector<char> data(fileSize);
  std::ifstream     is(hashtableConfigFile.string());
  if ((!is.read(data.data(), fileSize)) || (static_cast<long>(fileSize) != is.gcount())) {
    boost::format message =
        boost::format("ERROR: failed to read %i bytes from binary hashtable config: %i bytes read from %s") %
        fileSize % is.gcount() % hashtableConfigFile.string();
    BOOST_THROW_EXCEPTION(IoException(errno, message.str()));
  }
  return data;
}

template <typename T>
std::unique_ptr<T, std::function<void(T*)>> ReferenceDir7::mmapData(
    const std::string binFile, const size_t expectedBinFileBytes) const
{
  using namespace dragenos::common;
  checkDirectoryAndFile(path_, binFile);
  const auto dataFile = path_ / binFile;
  const auto fileSize = getFileSize(dataFile);
  if (fileSize != expectedBinFileBytes) {
    boost::format message =
        boost::format(
            "ERROR: file size different from size in config file: %s: expected %i bytes: actual %i bytes") %
        dataFile % expectedBinFileBytes % fileSize;
    BOOST_THROW_EXCEPTION(IoException(errno, message.str()));
  }
  if (0 == expectedBinFileBytes) {
    return nullptr;
  }
  const auto hashtableFd = open(dataFile.c_str(), O_RDONLY, 0);
  if (-1 == hashtableFd) {
    BOOST_THROW_EXCEPTION(
        IoException(errno, std::string("ERROR: failed to open data file ") + dataFile.string()));
  }
  const int prot   = PROT_READ;
  const int flags  = MAP_PRIVATE | MAP_NORESERVE;
  const int offset = 0;
  auto      table  = mmap(NULL, fileSize, prot, flags, hashtableFd, offset);
  if (MAP_FAILED == table) {
    BOOST_THROW_EXCEPTION(
        IoException(errno, std::string("ERROR: failed to map hashtable data file ") + dataFile.string()));
  }
  return std::unique_ptr<T, std::function<void(T*)>>(
      reinterpret_cast<T*>(table), [this](T* p) -> void { munmap(p, hashtableConfig_.getHashtableBytes()); });
}

template <typename T>
std::unique_ptr<T, std::function<void(T*)>> ReferenceDir7::readData(
    const std::string binFile, const size_t expectedBinFileBytes) const
{
  using namespace dragenos::common;
  checkDirectoryAndFile(path_, binFile);
  const auto dataFile = path_ / binFile;
  const auto fileSize = getFileSize(dataFile);
  if (fileSize != expectedBinFileBytes) {
    boost::format message =
        boost::format(
            "ERROR: file size different from size in config file: %s: expected %i bytes: actual %i bytes") %
        dataFile % expectedBinFileBytes % fileSize;
    BOOST_THROW_EXCEPTION(IoException(errno, message.str()));
  }
  if (0 == expectedBinFileBytes) {
    return nullptr;
  }
  const auto hashtableFd = open(dataFile.c_str(), O_RDONLY, 0);
  if (-1 == hashtableFd) {
    BOOST_THROW_EXCEPTION(
        IoException(errno, std::string("ERROR: failed to open data file ") + dataFile.string()));
  }
  char* table = (char*)malloc(fileSize);
  if (!table) {
    BOOST_THROW_EXCEPTION(std::bad_alloc());
  }
  auto    toRead    = fileSize;
  ssize_t bytesRead = 0;
  do {
    bytesRead = read(hashtableFd, table + fileSize - toRead, toRead);
    toRead -= bytesRead;
  } while (bytesRead);
  close(hashtableFd);
  if (toRead) {
    free(table);
    BOOST_THROW_EXCEPTION(IoException(
        errno,
        std::string("ERROR: failed to read ") + std::to_string(fileSize) + " bytes from " +
            dataFile.string() + " read: " + std::to_string(fileSize - toRead) +
            " error: " + std::strerror(errno)));
  }
  return std::unique_ptr<T, std::function<void(T*)>>(
      reinterpret_cast<T*>(table), [](T* p) -> void { free(p); });
}

size_t ReferenceDir7::getHashtableConfigSize() const
{
  return hashtableConfigData_.size();
}
size_t ReferenceDir7::getHashtableDataSize() const
{
  return hashtableConfig_.getHashtableBytes();
}

}  // namespace reference
}  // namespace dragenos
