#include "seed_generator.hh"

#include <limits>
#include <mutex>

namespace seed_generator
{
  std::mutex mut;

  std::mt19937 rng;

  std::uniform_int_distribution<std::mt19937::result_type> seed_generator(std::numeric_limits<std::mt19937::result_type>::min(),
                                                                          std::numeric_limits<std::mt19937::result_type>::max());

  std::mt19937::result_type
  get_new_seed()
  {
    std::lock_guard<std::mutex> lk{mut};
    return seed_generator(rng);
  }

  void
  set_meta_seed(std::mt19937::result_type s)
  {
    std::lock_guard<std::mutex> lk{mut};
    rng.seed(s);
  }
};
