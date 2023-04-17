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

#ifndef SYM_SCHUR_DECOMP_H
#define SYM_SCHUR_DECOMP_H

#include "SylvMatrix.hh"

class SymSchurDecomp
{
protected:
  Vector lambda;
  SqSylvMatrix q;
public:
  /* Computes the factorization A = Q·Λ·Qᵀ, where A is assummed to be
     symmetric and Λ real diagonal, hence a vector. */
  SymSchurDecomp(const ConstGeneralMatrix &a);
  SymSchurDecomp(const SymSchurDecomp &ssd) = default;
  virtual ~SymSchurDecomp() = default;
  const Vector &
  getLambda() const
  {
    return lambda;
  }
  const SqSylvMatrix &
  getQ() const
  {
    return q;
  }
  /* Return factor F·Fᵀ = A, raises and exception if A is not
     positive semidefinite, F must be square. */
  void getFactor(GeneralMatrix &f) const;
  // Returns true if A is positive semidefinite.
  bool isPositiveSemidefinite() const;
  /* Correct definitness. This sets all eigenvalues between minus
     tolerance and zero to zero. */
  void correctDefinitness(double tol);
};

#endif
