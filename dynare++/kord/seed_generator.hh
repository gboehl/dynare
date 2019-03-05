#ifndef RANDOM_H
#define RANDOM_H

#include <random>

namespace seed_generator
{
  // Produces seeds that can be used with Mersenne-Twister generators (thread-safe)
  std::mt19937::result_type get_new_seed();

  // Sets the seed for the seed generator (!)
  void set_meta_seed(std::mt19937::result_type s);
};

#endif
