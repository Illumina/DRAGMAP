// Copyright 2019 Edico Genome Corporation. All rights reserved.
//
// This file contains confidential and proprietary information of the Edico Genome
// Corporation and is protected under the U.S. and international copyright and other
// intellectual property laws.
//
// $Id$
// $Author$
// $Change$
// $DateTime$
//

#pragma once

#include <assert.h>

#include <algorithm>
#include <condition_variable>
#include <functional>
#include <iostream>
#include <mutex>
#include <thread>
#include <vector>

namespace dragenos {
namespace common {

namespace detail {
/**
 * \brief Base non-template class for ScopeEndCall
 */
class ScopeEndCallBase {
public:
  virtual ~ScopeEndCallBase() = default;

  explicit operator bool() const { return false; }
};

/**
 * \brief Holds a functor which is called at the destruction time
 */
template <typename FuncT>
class ScopeEndCall : public ScopeEndCallBase {
  FuncT f_;

public:
  explicit ScopeEndCall(FuncT f) : f_(f) {}

  ~ScopeEndCall() override { f_(std::uncaught_exceptions()); }
};

/**
 * \brief Helper to create ScopeEndCallHolder
 */
template <typename FuncT>
const detail::ScopeEndCall<FuncT> makeScopeEndCallHolder(FuncT f)
{
  return detail::ScopeEndCall<FuncT>(f);
}

}  // namespace detail

/**
 * \brief ensures f is called during the stack unwind of the scope following the macro
 */
#define ASYNC_BLOCK_WITH_CLEANUP(f)                              \
  if (const dragenos::common::detail::ScopeEndCallBase& b =      \
          dragenos::common::detail::makeScopeEndCallHolder(f)) { \
    (void)b;                                                     \
  } else

//-------------------------------------------------------------------------------rpetrovski
//   Inversion of the std::unique_lock and such
//
template <typename Lock>
class unlock_guard {
private:
  Lock& l;

public:
  unlock_guard(unlock_guard&) = delete;
  unlock_guard& operator=(unlock_guard&) = delete;
  explicit unlock_guard(Lock& m_) : l(m_) { l.unlock(); }

  ~unlock_guard() { l.lock(); }
};

struct ThreadPoolException : public std::logic_error {
  using std::logic_error::logic_error;
};

template <bool crashOnExceptions>
class BasicThreadPool {
public:
  typedef std::unique_lock<std::mutex> lock_type;

private:
  struct Executor {
    Executor(const std::size_t maxThreads) : maxThreads_(maxThreads) {}
    virtual void execute(lock_type& lock) = 0;
    virtual ~Executor()                   = default;

    const std::size_t    maxThreads_;
    Executor*            next_      = 0;
    Executor*            prev_      = 0;
    std::size_t          threadsIn_ = 0;
    bool                 complete_  = false;
    friend std::ostream& operator<<(std::ostream& os, const Executor& e)
    {
      return os << "Executor(" << e.maxThreads_ << "mt," << e.threadsIn_ << "ti," << e.complete_ << "c)";
    }
  } * head_, *tail_;
  static std::mutex       mutex_;
  std::condition_variable stateChangedCondition_;

  // true when the whole thing goes down
  bool terminateRequested_;

  // number of threads that have entered ready state
  std::size_t threadsReady_ = 0;

  std::exception_ptr firstThreadException_;
  bool               aThreadFailed_ = false;

  typedef std::vector<std::thread> ThreadVector;
  typedef ThreadVector::size_type  size_type;
  ThreadVector                     threads_{};

public:
  /**
   * \return number of threads in th pool + 1. This is because the thread calling execute()
   *         is counted as a worker.
   */
  std::size_t size() const { return threads_.size() + 1; }

  /**
   * \brief clears and repopulates ThreadPool. Use at your own risk. The only reasonable situation you'd might
   *        want to do this is in the sequence of unit tests that require different parallelization.
   *        also resets terminateRequested_ state.
   */
  void reset(std::size_t newSize)
  {
    clear();
    lock_type lock(mutex_);
    assert(0 < newSize);  //, "Inadequate pool size";
    // thread calling the execute will be one of the workers
    while (--newSize) {
      threads_.push_back(std::thread(&BasicThreadPool::threadFunc, this));
    }
  }

  /**
   * \brief Constructs a vector of size threads. All memory allocations that are required happen
   *        at this point.
   *
   *  @param size         number of threads to produce
   */
  explicit BasicThreadPool(size_type size = std::thread::hardware_concurrency())
    : head_(0), tail_(0), terminateRequested_(false)
  {
    reset(size);
  }

  /**
   * \brief Tells all threads to terminate and releases them.
   */
  ~BasicThreadPool() { clear(); }

  template <typename F>
  void execute(F func, const unsigned threads)
  {
    lock_type lock(mutex_);
    execute(lock, func, threads);
  }

  /**
   * \brief Executes func on requested number of threads.
   *
   * \threads number of threads to use. This must be less or equal than size()
   * \param func Unary function to execute.
   */
  template <typename F>
  void execute(lock_type& lock, F func, const unsigned threads)
  {
    if (aThreadFailed_) {
      throw ThreadPoolException("Attempt to execute a failed threadpool");
    }
    if (firstThreadException_) {
      //           std::cerr << "[" << std::this_thread::get_id() << "]\tWARNING: execute called when an exception is pending " << std::endl;
      std::rethrow_exception(firstThreadException_);
    }

    //       assert(threads <= size()); //, "Request must not exceed the amount of threads available");
    struct FuncExecutor : public Executor {
      F& func_;
      FuncExecutor(F& func, const unsigned threads) : Executor(threads), func_(func) {}
      virtual void execute(lock_type& lock) { func_(lock); }
    } executor(func, threads);

    //       std::cerr << "created " << executor << std::endl;

    enque(&executor);
    //       std::cerr << "enqued " << executor << std::endl;

    stateChangedCondition_.notify_all();

    if (crashOnExceptions) {
      unsafeExecute(lock, &executor);
    } else {
      safeExecute(lock, &executor);
    }

    // thread does not leave an executor unless there is nothing more to be done
    assert(executor.complete_);

    // wait for helper threads to exit
    while (!terminateRequested_ && (executor.threadsIn_)) {
      stateChangedCondition_.wait(lock);
    }

    unque(&executor);
    //       std::cerr << "unqued " << executor << std::endl;

    if (firstThreadException_) {
      //           std::cerr << "[" << std::this_thread::get_id() << "]\tWARNING: rethrowing a thread exception " << std::endl;
      std::rethrow_exception(firstThreadException_);
    }
  }

  /*
   * \return always false or throws an exception
   */
  bool checkThreadFailed() const
  {
    if (aThreadFailed_) {
      throw ThreadPoolException("Another thread failed in the threadpool");
    }

    return false;
  }

  void waitForChange(lock_type& lock)
  {
    if (firstThreadException_) {
      //       std::cerr << "[" << std::this_thread::get_id() << "]\tWARNING: waitForChange called when an exception is pending " << std::endl;
      std::rethrow_exception(firstThreadException_);
    }

    if (aThreadFailed_) {
      throw ThreadPoolException("Attempt to wait in a failed threadpool");
    }

    if (!terminateRequested_) {
      stateChangedCondition_.wait(lock);
    }

    if (firstThreadException_) {
      //       std::cerr << "[" << std::this_thread::get_id() << "]\tWARNING: rethrowing a thread exception " << std::endl;
      std::rethrow_exception(firstThreadException_);
    }

    if (aThreadFailed_) {
      throw ThreadPoolException("Another thread failed while waiting for state change");
    }
  }

  void notify_all() { stateChangedCondition_.notify_all(); }

  /**
   * \brief Executes func on requested number of size() threads.
   **/
  template <typename F>
  void execute(F func)
  {
    execute(func, static_cast<const unsigned int>(size()));
  }

  lock_type lock() { return lock_type(mutex_); }

private:
  void clear()
  {
    {
      lock_type lock(mutex_);
      // if clear is called immediately after construction, some threads might not be ready to be joined
      waitAllReady(lock);

      terminateRequested_ = true;
      stateChangedCondition_.notify_all();
    }
    std::for_each(threads_.begin(), threads_.end(), [](std::thread& t) { t.join(); });
    threads_.clear();
    terminateRequested_ = false;
    threadsReady_       = 0;
  }

  void unque(Executor* executor)
  {
    if (executor->next_) {
      executor->next_->prev_ = executor->prev_;
    } else {
      tail_ = executor->prev_;
    }
    if (executor->prev_) {
      executor->prev_->next_ = executor->next_;
    } else {
      head_ = executor->next_;
    }
  }

  void enque(Executor* executor)
  {
    executor->next_ = head_;
    if (head_) {
      head_->prev_ = executor;
    } else {
      tail_ = executor;
    }
    head_ = executor;
  }

  Executor* findNext()
  {
    Executor* ret = 0;
    for (Executor* e = tail_; 0 != e; e = e->prev_) {
      if (!e->complete_ && e->threadsIn_ < e->maxThreads_) {
        ret = e;
        break;
      }
    }
    return ret;
  }
  /**
   * \brief executes and allows the exception to escape
   */
  void unsafeExecute(lock_type& lock, Executor* executor)
  {
    assert(!executor->complete_);
    assert(executor->threadsIn_ < executor->maxThreads_);

    ++executor->threadsIn_;
    try {
      executor->execute(lock);
    } catch (...) {
      aThreadFailed_      = true;
      executor->complete_ = true;
      --executor->threadsIn_;
      stateChangedCondition_.notify_all();
      throw;
    }
    executor->complete_ = true;
    --executor->threadsIn_;
    stateChangedCondition_.notify_all();
  }

  /**
   * \brief executes and stores the exception information so that it can be rethrown on the main thread
   */
  void safeExecute(lock_type& lock, Executor* executor)
  {
    try {
      unsafeExecute(lock, executor);
    } catch (...) {
      if (!firstThreadException_) {
        firstThreadException_ = std::current_exception();
        std::cerr << "[" << std::this_thread::get_id() << "]\tERROR: This thread caught an exception first"
                  << std::endl;
      } else {
        //             std::cerr << "[" << std::this_thread::get_id() << "]\tERROR: This thread also caught an exception" << std::endl;
      }
    }
  }

  void threadFunc()
  {
    lock_type lock(mutex_);
    ++threadsReady_;
    stateChangedCondition_.notify_all();
    while (!terminateRequested_) {
      executeAllPending(lock);
      if (!terminateRequested_) {
        stateChangedCondition_.wait(lock);
      }
    }
  }

  void waitAllReady(lock_type& lock)
  {
    while (size() != (threadsReady_ + 1)) {
      stateChangedCondition_.wait(lock);
    }
  }

  bool executeAllPending(lock_type& lock)
  {
    bool ret = false;
    for (Executor* executor = findNext(); !terminateRequested_ && 0 != executor; executor = findNext()) {
      if (crashOnExceptions) {
        unsafeExecute(lock, executor);
      } else {
        safeExecute(lock, executor);
      }
      ret = true;
    }

    return ret;
  }
};
template <bool crashOnExceptions>
std::mutex BasicThreadPool<crashOnExceptions>::mutex_;

typedef BasicThreadPool<false> SafeThreadPool;
typedef BasicThreadPool<true>  UnsafeThreadPool;
typedef SafeThreadPool         ThreadPool;

/**
 * \param threadsMax  maximum number of threads. Used in the singleton initialization.
 *                    will produce failure if changes on subsequent calls.
 *                    special value 0 will result in pool allocation with std::thread::hardware_concurrency(),
 *                    also will not lead to a failure on subsequent calls regardless of the size of allocated
 * pool \return reference to a thread pool singleton
 */
ThreadPool& CPU_THREADS(std::size_t threadsMax = 0);

}  // namespace common
}  // namespace dragenos
