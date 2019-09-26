/*
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
 * along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
 */

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
