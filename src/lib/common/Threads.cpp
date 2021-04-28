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

#include "common/Threads.hpp"

namespace dragenos {

namespace common {

ThreadPool& CPU_THREADS(std::size_t threadsMax /* = 0*/)
{
  static ThreadPool pool(threadsMax ? threadsMax : std::thread::hardware_concurrency());

  // if does not match, means that static was constructed with a different parameter
  assert(!threadsMax || pool.size() == threadsMax);

  return pool;
}

}  // namespace common
}  // namespace dragenos
