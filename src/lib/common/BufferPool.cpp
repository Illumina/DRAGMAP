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

#include "common/BufferPool.hpp"

#include <cerrno>
#include <cstring>

#include <boost/assert.hpp>
#include <boost/format.hpp>

#include <cstring>
#include "common/DragenConstants.hpp"
#include "common/Exceptions.hpp"
//#include "dragen_exception.hpp"

/* ============ end of #include directives ============ */
//#include "edico_memdebug.h"

namespace dragenos {
namespace common {

//--------------------------------------------------------------------------------adam
// constructor - sets up the pool freelist of the specified size
BufferPool::BufferPool(
    uint32_t numBuffers,  // how many buffers to put in the pool
    size_t   bufferSize)    // size of each buffer
  : m_bufferSize(bufferSize),
    m_totalSize(numBuffers * bufferSize),
    m_totalBorrowed(0),
    m_highWatermark(0),
    m_pendingReplaces(0),
    m_exit(false),
    m_preallocated(false)
{
  BOOST_ASSERT(m_totalSize > 0);
}

//--------------------------------------------------------------------------------adam
// destructor - free the entire pool.
BufferPool::~BufferPool()
{
  if (m_preallocated) {
    DeactivatePreallocation();
  }
}

//------------------------------------------------------------------------ alain
bool BufferPool::SetPreallocation(bool preallocateFlag)
{
  std::unique_lock<std::mutex> lock(m_mutex);
  bool                         prevMode = m_preallocated;

  // ensure no allocation was made
  BOOST_ASSERT(m_totalBorrowed == 0);

  if (m_preallocated != preallocateFlag) {
    // there is a change of mode
    if (m_preallocated) {
      // was preallocated, now non-preallocated
      DeactivatePreallocation();

    } else {
      // was non-preallocated, now preallocated
      ActivatePreallocation();
    }

  } else {
    // no mode change, nothing to do
  }

  return prevMode;
}

//------------------------------------------------------------------------ alain
// preallocate buffers
void BufferPool::ActivatePreallocation()
{
  Buffer   buf;
  uint32_t numBuffers = m_totalSize / m_bufferSize;

  BOOST_ASSERT(m_preallocatedBuffers.size() == 0);

  while (m_preallocatedBuffers.size() != numBuffers) {
    buf = NULL;
    if (!posix_memalign(reinterpret_cast<void**>(&buf), DRAGEN_PAGE_SIZE, m_bufferSize)) {
      BOOST_THROW_EXCEPTION(
          MemoryException(std::string("Failed to preallocate buffers: ") + strerror(errno)));
    }
    m_preallocatedBuffers.push(buf);
  }

  BOOST_ASSERT(m_preallocatedBuffers.size() == numBuffers);

  m_preallocated = true;
}

//------------------------------------------------------------------------ alain
// free up preallocated buffers
void BufferPool::DeactivatePreallocation()
{
  while (m_preallocatedBuffers.size() > 0) {
    Buffer buf = reinterpret_cast<Buffer>(m_preallocatedBuffers.front());
    m_preallocatedBuffers.pop();
    free(buf);
  }

  BOOST_ASSERT(m_preallocatedBuffers.size() == 0);

  m_preallocated = false;
}

//--------------------------------------------------------------------------------adam
// borrow -- wait until the pool is non-empty, and return a buffer from the freelist
//
BufferPool::Buffer BufferPool::borrow(const bool wait, const bool clear)
{
  Buffer b = NULL;

  {
    // If necessary, wait until there is a buffer available
    std::unique_lock<std::mutex> lock(m_mutex);
    while (m_totalBorrowed >= m_totalSize) {
      if (!wait) {
        // return immediately with invalid value
        return NULL;
      }
      m_cond.wait(lock);
    }
    m_totalBorrowed += m_bufferSize;

    // Update the max borrowed if needed
    if (m_highWatermark < m_totalBorrowed) m_highWatermark = m_totalBorrowed;

    if (m_preallocated) {
      // get from the pool
      b = reinterpret_cast<Buffer>(m_preallocatedBuffers.front());
      m_preallocatedBuffers.pop();
      if (clear) {
        memset(b, 0, m_bufferSize);
      }
      return b;
    }
  }

  if (m_preallocated) {
    // b already is set in the critical section above
  } else {
    // dynamically allocate
    if (!posix_memalign(reinterpret_cast<void**>(&b), DRAGEN_PAGE_SIZE, m_bufferSize)) {
      BOOST_THROW_EXCEPTION(MemoryException(std::string("Failed to allocate buffer: ") + strerror(errno)));
    }
  }

  // alaing: add "and b" to the if condition to make SCA happy.
  // In the report, clang SCA assumes m_preallocated could have changed value
  // between the test in locked section and the test outside the locked section
  // which we know won't happen because mode is not changed during the lifetime
  // of the object
  if (clear and b) {
    // JimmyB make this an option for e.g. if used in the VC where the user fills the
    // buffer this seems to be unnecesarry
    //
    // You may wonder whether this memset is necessary, and is too expensive.  Well,
    // I tried commenting it out, and turns out it *is* necessary. I think because
    // we are relying on it to initialize the dbam header scructures.
    memset(b, 0, m_bufferSize);
  }
  return b;
}

//--------------------------------------------------------------------------------adam
// replace - put #buf# back into the freelist
//
void BufferPool::replace(void* abuf)
{
  if (abuf ==
      0) {  // Can occur when a buffer must be given to give 'last buffer' signal, but no data is available
    m_cond.notify_one();  // Wake up other threads in case anyone is still waiting
    return;
  }

  Buffer buf = reinterpret_cast<Buffer>(abuf);

  if (m_preallocated) {
    // return to the pool in the critical section below
  } else {
    // free dynamically-allocated buffer
    free(buf);
  }

  {
    // Return the buffer to the tail of the freelist
    std::lock_guard<std::mutex> lock(m_mutex);
    BOOST_ASSERT_MSG(
        m_totalBorrowed >= m_bufferSize,
        (boost::format("total borrowed: %i, buffer size: %i") % m_totalBorrowed % m_bufferSize)
            .str()
            .c_str());
    if (m_pendingReplaces >= m_bufferSize) {
      m_pendingReplaces -= m_bufferSize;
    }

    if (m_preallocated) {
      // back to the pool
      m_preallocatedBuffers.push(buf);
    }
  }

  // Let anyone waiting know that there's a buffer available now.
  m_cond.notify_one();
}
}  // namespace common
}  // namespace dragenos
