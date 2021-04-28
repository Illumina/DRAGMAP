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

#ifndef COMMON_BUFFER_POOL_HPP
#define COMMON_BUFFER_POOL_HPP

#include <condition_variable>
#include <deque>
#include <mutex>
#include <queue>
#include "inttypes.h"

namespace dragenos {
namespace common {

class BufferPool {
private:
  enum { BUFFER_ALIGNMENT = 1 << 12 };

private:
  typedef uint8_t*           Buffer;
  typedef std::deque<Buffer> Buffer_c;
  typedef Buffer_c::iterator Buffer_it;

  typedef std::queue<Buffer> Buffer_q;

public:
  BufferPool(uint32_t numBuffers, size_t bufferSize);

  virtual ~BufferPool();

  // Borrow a buffer from the pool
  Buffer borrow(const bool wait = true, const bool clear = true);  // borrow a buffer from the pool

  // Put a buffer back on the freelist
  void replace(void* buf);  // the buffer to replace in the freelist

  // trigger our condition variable to allow users to wake up and exit
  void exit()
  {
    std::unique_lock<std::mutex> lock(m_mutex);
    m_exit = true;
    m_cond.notify_one();
  }

  void unexit()
  {
    std::unique_lock<std::mutex> lock(m_mutex);
    m_exit = false;
  }

  // Return the number of free buffers in the buffer pool
  size_t getNumFree()
  {
    std::lock_guard<std::mutex> lock(m_mutex);
    return (m_totalSize - m_totalBorrowed) / m_bufferSize;
  }
  // Return the number of free buffers in the buffer pool (debug only)
  // Caller shall accept stale values
  size_t getNumFreeDebugOnly() { return (m_totalSize - m_totalBorrowed) / m_bufferSize; }

  // Return the total number of buffers in the buffer pool
  uint32_t getNumBuffers() { return m_totalSize / m_bufferSize; }

  // Add #num_buffers# worth of capacity to the pool
  void addBufferCapacity(size_t num_buffers)
  {
    m_totalSize += (m_bufferSize * num_buffers);
    m_cond.notify_one();
  }

  // Returns when a free buffer is avaiable.
  void waitForFreeBuffer()
  {
    std::unique_lock<std::mutex> lock(m_mutex);
    while (m_totalBorrowed + m_bufferSize >= m_totalSize && !m_exit) {  //m_free.empty() ) {
      m_cond.wait(lock);
    }
  }

  // Return size of each buffer
  size_t getBufferSize() const { return m_bufferSize; }

  // Return number of bytes currently on loan (i.e. borrowed)
  size_t getTotalBorrowed() const { return m_totalBorrowed; }

  // Get the high water mark for the buffer pool
  size_t getHighWatermark() const { return m_highWatermark; }

  // Get number of bytes still on loan, but about to be replace()'d due to being in an
  // i/o queue or some other such latency.
  size_t getPendingReplaces() const { return m_pendingReplaces; }

  // Take note that we should expect #n# bytes to be coming back via replace() soon.
  void incPendingReplaces(const size_t n)
  {
    std::unique_lock<std::mutex> lock(m_mutex);
    m_pendingReplaces += n;
  }

  // Return total amount of memory (both borrowed and available) in the pool
  size_t getTotalSize() const { return m_totalSize; }

  bool SetPreallocation(bool preallocateFlag);

private:
  size_t                  m_bufferSize;       // Size of each buffer
  size_t                  m_totalSize;        // Total number of bytes used for our buffers
  size_t                  m_totalBorrowed;    // How many bytes have been borrow()'ed
  size_t                  m_highWatermark;    // The high water mark for the Buffer pool
  size_t                  m_pendingReplaces;  // How many bytes are about to be returned
  std::mutex              m_mutex;            // Protect changes to the free/allocated lists
  std::condition_variable m_cond;             // Signals freelist non-emptiness
  bool                    m_exit;

  bool     m_preallocated;
  Buffer_q m_preallocatedBuffers;

  void ActivatePreallocation();
  void DeactivatePreallocation();
};

}  // namespace common
}  // namespace dragenos
#endif  // #ifndef COMMON_BUFFER_POOL_HPP
