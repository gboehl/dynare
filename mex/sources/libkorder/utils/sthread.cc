/*
 * Copyright © 2004 Ondra Kamenik
 * Copyright © 2019 Dynare Team
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

#include <algorithm>

#include "sthread.hh"

namespace sthread
{
  /* We set the default value for ‘max_parallel_threads’ to half the number of
     logical CPUs */
  int
  default_threads_number()
  {
    return std::max(1, static_cast<int>(std::thread::hardware_concurrency()) / 2);
  }

  int detach_thread_group::max_parallel_threads = default_threads_number();

  /* We cycle through all threads in the group, and in each cycle we wait for
     the change in the ‘counter’. If the counter indicates less than maximum
     parallel threads running, then a new thread is run, and the iterator in
     the list is moved.

     At the end we have to wait for all thread to finish. */
  void
  detach_thread_group::run()
  {
    std::unique_lock<std::mutex> lk{mut_cv};
    auto it = tlist.begin();
    while (it != tlist.end())
      {
        counter++;
        std::thread th{[&, it] {
                         // The ‘it’ variable is captured by value, because otherwise the iterator may move
                         (*it)->operator()(mut_threads);
                         std::unique_lock<std::mutex> lk2{mut_cv};
                         counter--;
                         /* First notify the thread waiting on the condition variable, then
                            unlock the mutex. We must do these two operations in that order,
                            otherwise there is a possibility of having the main process
                            destroying the condition variable before the thread tries to
                            notify it (if all other threads terminate at the same time and
                            bring the counter down to zero).
                            For that reason, we cannot use std::notify_all_at_thread_exit() */
                         cv.notify_one();
                         lk2.unlock();
                       }};
        th.detach();
        ++it;
        cv.wait(lk, [&] { return counter < max_parallel_threads; });
      }
    cv.wait(lk, [&] { return counter == 0; });
  }
}
