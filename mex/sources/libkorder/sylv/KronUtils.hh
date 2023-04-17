/*
 * Copyright © 2004-2011 Ondra Kamenik
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

#ifndef KRON_UTILS_H
#define KRON_UTILS_H

#include "KronVector.hh"
#include "QuasiTriangular.hh"

class KronUtils
{
public:
  /* Computes x = (Iₘ⊗…⊗Iₘ⊗T⊗Iₘ⊗…⊗Iₘ⊗Iₙ)·x, where x is n×mᵈ.
     T must be m×m, number of ⊗ is d, level is the number of Iₘ’s
     between T and Iₙ plus 1. If level=0, then we multiply by Iₘ⊗…⊗Iₘ⊗T,
     T must be n×n. */
  static void multAtLevel(int level, const QuasiTriangular &t,
                          KronVector &x);
  static void multAtLevelTrans(int level, const QuasiTriangular &t,
                               KronVector &x);

  // Computes x=(Fᵀ⊗Fᵀ⊗…⊗K)·x
  static void multKron(const QuasiTriangular &f, const QuasiTriangular &k,
                       KronVector &x);
};

#endif /* KRON_UTILS_H */
