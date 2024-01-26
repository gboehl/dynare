/* Interface to quasi Monte Carlo sequences (à la Sobol) routines.
 *
 * Copyright © 2010-2024 Dynare Team
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

#ifndef SOBOL_HH
#define SOBOL_HH

#include <concepts>
#include <cstdint>
#include <dynblas.h> // For the FORTRAN_WRAPPER macro

constexpr int DIM_MAX = 1111;

// Interface to Fortran code
extern "C"
{
#define i8_sobol FORTRAN_WRAPPER(i8_sobol)
  void i8_sobol(const int64_t* dim_num, int64_t* seed, double* quasi);
}

inline void
next_sobol(int dim_num, int64_t* seed, double* quasi)
{
  int64_t dim_num2 {dim_num};
  i8_sobol(&dim_num2, seed, quasi);
}

inline int64_t
sobol_block(int dimension, int block_size, int64_t seed, double* block)
{
  for (int iter = 0; iter < block_size; iter++)
    next_sobol(dimension, &seed, &block[iter * dimension]);
  return seed;
}

template<floating_point T>
void
expand_unit_hypercube(int dimension, int block_size, T* block, const T* lower_bound,
                      const T* upper_bound)
{
  T* hypercube_length = new T[dimension];
  for (int dim = 0; dim < dimension; dim++)
    hypercube_length[dim] = upper_bound[dim] - lower_bound[dim];

  int base = 0;
  for (int sim = 0; sim < block_size; sim++)
    {
      for (int dim = 0; dim < dimension; dim++)
        block[base + dim] = lower_bound[dim] + hypercube_length[dim] * block[base + dim];
      base += dimension;
    }
  delete[] hypercube_length;
}

#endif
