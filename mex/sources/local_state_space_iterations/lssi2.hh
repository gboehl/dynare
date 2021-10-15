// Prototypes for state space iteration at order 2 routines written in assembly

/*
 * Copyright © 2022 Dynare Team
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

/*
  WARNING: s must be a multiple of 4 for the AVX2 routines.

  Sizes of the vector/matrix arguments:
  – y: m×s
  – yhat: n×s
  – epsilon: q×s
  – ghx: m×n
  – constant: m
  – ghxx: m×n²
  – ghuu: m×q²
  – ghxu: m×nq
*/

#if defined(__x86_64__) && defined(__LP64__)

extern "C" void lssi2_avx2(double *y, const double *yhat, const double *epsilon,
                           const double *ghx, const double *ghu, const double *constant,
                           const double *ghxx, const double *ghuu, const double *ghxu,
                           int m, int n, int q, int s);

extern "C" void lssi2_avx512(double *y, const double *yhat, const double *epsilon,
                             const double *ghx, const double *ghu, const double *constant,
                             const double *ghxx, const double *ghuu, const double *ghxu,
                             int m, int n, int q, int s);

#endif
