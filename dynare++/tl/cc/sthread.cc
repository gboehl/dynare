// Copyright 2004, Ondra Kamenik


#include "sthread.hh"

namespace sthread
{
  /* We set the default value for |max_parallel_threads| to 2, i.e.
     uniprocessor machine with hyper-threading */
  int detach_thread_group::max_parallel_threads = 2;

  /* We cycle through all threads in the group, and in each cycle we wait
     for the change in the |counter|. If the counter indicates less than
     maximum parallel threads running, then a new thread is run, and the
     iterator in the list is moved.

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
            // The "it" variable is captured by value, because otherwise the iterator may move
            (*it)->operator()(mut_threads);
            std::unique_lock<std::mutex> lk2{mut_cv};
            counter--;
            std::notify_all_at_thread_exit(cv, std::move(lk2));
          }};
        th.detach();
        ++it;
        cv.wait(lk, [&] { return counter < max_parallel_threads; });
      }
    cv.wait(lk, [&] { return counter == 0; });
  }
}
