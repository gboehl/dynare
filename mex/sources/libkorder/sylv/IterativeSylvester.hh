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

#ifndef ITERATIVE_SYLVESTER_H
#define ITERATIVE_SYLVESTER_H

#include "SylvesterSolver.hh"
#include "KronVector.hh"
#include "QuasiTriangular.hh"
#include "SimilarityDecomp.hh"

class IterativeSylvester : public SylvesterSolver
{
public:
  IterativeSylvester(const QuasiTriangular &k, const QuasiTriangular &f)
    : SylvesterSolver(k, f)
  {
  }
  IterativeSylvester(const SchurDecompZero &kdecomp, const SchurDecomp &fdecomp)
    : SylvesterSolver(kdecomp, fdecomp)
  {
  }
  IterativeSylvester(const SchurDecompZero &kdecomp, const SimilarityDecomp &fdecomp)
    : SylvesterSolver(kdecomp, fdecomp)
  {
  }
  void solve(SylvParams &pars, KronVector &x) const override;
private:
  double performFirstStep(KronVector &x) const;
  static double performStep(const QuasiTriangular &k, const QuasiTriangular &f,
                            KronVector &x);
};

#endif /* ITERATIVE_SYLVESTER_H */
