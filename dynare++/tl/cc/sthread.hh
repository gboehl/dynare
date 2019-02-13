// Copyright 2004, Ondra Kamenik

// Simple threads.

/* This file defines types making a simple interface to
   multi-threading.

   The file provides the following interfaces:
   \unorderedlist
   \li |detach_thread| is a pure virtual class, which must be inherited and a
   method |operator()()| be implemented as the running code of the
   thread.
   \li |detach_thread_group| allows insertion of |detach_thread|s and running
   all of them simultaneously. The threads
   are not joined, they are synchronized by means of a counter counting
   running threads. A change of the counter is checked by waiting on an
   associated condition. The number of maximum parallel
   threads can be controlled. See below. The group also provides a mutex to be
   shared between the workers for their own synchronization purposes.
   \endunorderedlist

   The number of maximum parallel threads is controlled via a static
   member of the |detach_thread_group| class. */

#ifndef STHREAD_H
#define STHREAD_H

#include <vector>
#include <map>
#include <memory>
#include <utility>
#include <thread>
#include <mutex>
#include <condition_variable>

namespace sthread
{
  /* The detached thread is the same as joinable |thread|. We only
     re-implement |run| method to call |thread_traits::detach_run|, and add
     a method which installs a counter. The counter is increased and
     decreased on the body of the new thread. */

  class detach_thread
  {
  public:
    virtual ~detach_thread() = default;
    virtual void operator()(std::mutex &mut) = 0;
  };

  /* The detach thread group is (by interface) the same as
     |thread_group|. The extra thing we have here is the |counter|. The
     implementation of |insert| and |run| is different. */

  class detach_thread_group
  {
    std::vector<std::unique_ptr<detach_thread>> tlist;
    std::mutex mut_cv; // For the condition variable and the counter
    std::condition_variable cv;
    int counter{0};
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
