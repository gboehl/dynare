// Copyright 2004, Ondra Kamenik


#include "sthread.hh"

namespace sthread
{
  /* We set the default value for |max_parallel_threads| to 2, i.e.
     uniprocessor machine with hyper-threading */
  int detach_thread_group::max_parallel_threads = 2;

  /* The constructor acquires the mutex in the map. First it tries to
     get an exclusive access to the map. Then it increases a number of
     references of the mutex (if it does not exists, it inserts it). Then
     unlocks the map, and finally tries to lock the mutex of the map. */
  synchro::synchro(const void *c, std::string id)
    : caller{c}, iden{std::move(id)}
  {
    mutmap.lock_map();
    if (!mutmap.check(caller, iden))
      mutmap.insert(caller, iden);
    mutmap.get(caller, iden).second++;
    mutmap.unlock_map();
    mutmap.get(caller, iden).first.lock();
  }

  /* The destructor first locks the map. Then releases the lock,
     and decreases a number of references. If it is zero, it removes the
     mutex. */
  synchro::~synchro()
  {
    mutmap.lock_map();
    if (mutmap.check(caller, iden))
      {
        mutmap.get(caller, iden).first.unlock();
        mutmap.get(caller, iden).second--;
        if (mutmap.get(caller, iden).second == 0)
          mutmap.remove(caller, iden);
      }
    mutmap.unlock_map();
  }

  /* We cycle through all threads in the group, and in each cycle we wait
     for the change in the |counter|. If the counter indicates less than
     maximum parallel threads running, then a new thread is run, and the
     iterator in the list is moved.

     At the end we have to wait for all thread to finish. */
  void
  detach_thread_group::run()
  {
    std::unique_lock<std::mutex> lk{m};
    auto it = tlist.begin();
    while (it != tlist.end())
      {
        counter++;
        std::thread th{[&, it] {
            // The "it" variable is captured by value, because otherwise the iterator may move
            (*it)->operator()();
            std::unique_lock<std::mutex> lk2{m};
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
