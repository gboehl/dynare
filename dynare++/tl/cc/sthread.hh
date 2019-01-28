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
   threads can be controlled. See below.
   \li |synchro| object locks a piece of code to be executed only serially
   for a given data and specified entry-point. It locks the code until it
   is destructed. So, the typical use is to create the |synchro| object
   on the stack of a function which is to be synchronized. The
   synchronization can be subjected to specific data (then a pointer can
   be passed to |synchro|'s constructor), and can be subjected to
   specific entry-point (then |std::string| is passed to the
   constructor).
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
  using mmkey = std::pair<const void *, std::string>;

  /* Here we define a map of mutexes keyed by a pair of address, and a
     string. A purpose of the map of mutexes is that, if synchronizing, we
     need to publish mutexes locking some piece of codes (characterized by
     the string) accessing the data (characterized by the pointer). So, if
     any thread needs to pass a |synchro| object, it creates its own with
     the same address and string, and must look to some public storage to
     unlock the mutex. If the |synchro| object is created for the first
     time, the mutex is created and inserted to the map. We count the
     references to the mutex (number of waiting threads) to know, when it
     is save to remove the mutex from the map. This is the only purpose of
     counting the references. Recall, that the mutex is keyed by an address
     of the data, and without removing, the number of mutexes would only
     grow.

     The map itself needs its own mutex to avoid concurrent insertions and
     deletions. */

  struct ltmmkey
  {
    bool
    operator()(const mmkey &k1, const mmkey &k2) const
    {
      return k1.first < k2.first
                        || (k1.first == k2.first && k1.second < k2.second);
    }
  };

  using mutex_int_map = std::map<mmkey, std::pair<std::mutex, int>, ltmmkey>;

  class mutex_map : public mutex_int_map
  {
    using mmval = std::pair<std::mutex, int>;
    std::mutex m;
  public:
    mutex_map() = default;
    void
    insert(const void *c, std::string id)
    {
      // We cannot use emplace(), because std::mutex is neither copyable nor moveable
      operator[](mmkey{c, std::move(id)}).second = 0;
    }
    bool
    check(const void *c, std::string id) const
    {
      return find(mmkey{c, std::move(id)}) != end();
    }
    /* This returns the pair of mutex and count reference number. */
    mmval &
    get(const void *c, std::string id)
    {
      return operator[](mmkey{c, std::move(id)});
    }

    /* This removes unconditionally the mutex from the map regardless its
       number of references. The only user of this class should be |synchro|
       class, it implementation must not remove referenced mutex. */
    void
    remove(const void *c, std::string id)
    {
      auto it = find(mmkey{c, std::string{id}});
      if (it != end())
        erase(it);
    }
    void
    lock_map()
    {
      m.lock();
    }
    void
    unlock_map()
    {
      m.unlock();
    }

  };


  // The global map used by the synchro class
  static mutex_map mutmap;

  /* This is the |synchro| class. The constructor of this class tries to
     lock a mutex for a particular address (identification of data) and
     string (identification of entry-point). If the mutex is already
     locked, it waits until it is unlocked and then returns. The destructor
     releases the lock. The typical use is to construct the object on the
     stacked of the code being synchronized. */
  class synchro
  {
  private:
    const void *caller;
    const std::string iden;
  public:
    synchro(const void *c, std::string id);
    ~synchro();
  };

  /* The detached thread is the same as joinable |thread|. We only
     re-implement |run| method to call |thread_traits::detach_run|, and add
     a method which installs a counter. The counter is increased and
     decreased on the body of the new thread. */

  class detach_thread
  {
  public:
    virtual ~detach_thread() = default;
    virtual void operator()() = 0;
  };

  /* The detach thread group is (by interface) the same as
     |thread_group|. The extra thing we have here is the |counter|. The
     implementation of |insert| and |run| is different. */

  class detach_thread_group
  {
    std::vector<std::unique_ptr<detach_thread>> tlist;
    std::mutex m; // For the condition variable and the counter
    std::condition_variable cv;
    int counter{0};
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
};

#endif
