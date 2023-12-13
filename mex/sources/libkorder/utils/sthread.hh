/*
 * Copyright © 2004 Ondra Kamenik
 * Copyright © 2019-2023 Dynare Team
 *
 * This file is part of Dynare.
 *
 * Dynare is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Dynare is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Dynare.  If not, see <https://www.gnu.org/licenses/>.
 */

// Simple threads.

/* This file defines types making a simple interface to multi-threading.

   The file provides the following interfaces:

   — detach_thread is a pure virtual class, which must be inherited and a
     method operator()() be implemented as the running code of the thread.

   — detach_thread_group allows insertion of detach_thread’s and running all of
     them simultaneously. The threads are not joined, they are synchronized by
     means of a counter counting running threads. A change of the counter is
     checked by waiting on an associated condition. The number of maximum
     parallel threads can be controlled. See below. The group also provides a
     mutex to be shared between the workers for their own synchronization
     purposes.

   The number of maximum parallel threads is controlled via a static member of
   the detach_thread_group class. */

#ifndef STHREAD_HH
#define STHREAD_HH

#include <condition_variable>
#include <map>
#include <memory>
#include <mutex>
#include <thread>
#include <utility>
#include <vector>

namespace sthread
{
class detach_thread
{
public:
  virtual ~detach_thread() = default;
  virtual void operator()(std::mutex& mut) = 0;
};

class detach_thread_group
{
  std::vector<std::unique_ptr<detach_thread>> tlist;
  std::mutex mut_cv; // For the condition variable and the counter
  std::condition_variable cv;
  int counter {0};
  std::mutex mut_threads; // Passed to the workers and shared between them
public:
  static int max_parallel_threads;

  void
  insert(std::unique_ptr<detach_thread> c)
  {
    tlist.push_back(std::move(c));
  }

  ~detach_thread_group() = default;

  void run();
};

int default_threads_number();
};

#endif
