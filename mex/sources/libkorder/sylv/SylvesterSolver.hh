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

#ifndef SYLVESTER_SOLVER_H
#define SYLVESTER_SOLVER_H

#include "KronVector.hh"
#include "QuasiTriangular.hh"
#include "QuasiTriangularZero.hh"
#include "SimilarityDecomp.hh"
#include "SylvParams.hh"
#include "SchurDecomp.hh"

#include <memory>

class SylvesterSolver
{
protected:
  const std::unique_ptr<const QuasiTriangular> matrixK;
  const std::unique_ptr<const QuasiTriangular> matrixF;
private:
  /* Return true when it is more efficient to use QuasiTriangular
     than QuasiTriangularZero */
  static bool
  zeroPad(const SchurDecompZero &kdecomp)
  {
    return ((kdecomp.getZeroCols()*3 < kdecomp.getDim()*2)
            || (kdecomp.getZeroCols() < 10));
  }
public:
  SylvesterSolver(const QuasiTriangular &k, const QuasiTriangular &f)
    : matrixK(std::make_unique<QuasiTriangular>(k)),
      matrixF(std::make_unique<QuasiTriangular>(f))
  {
  }
  SylvesterSolver(const SchurDecompZero &kdecomp, const SchurDecomp &fdecomp)
    : matrixK((zeroPad(kdecomp)) ?
              std::make_unique<QuasiTriangular>(kdecomp) :
              std::make_unique<QuasiTriangularZero>(kdecomp)),
      matrixF(std::make_unique<QuasiTriangular>(fdecomp))
  {
  }
  SylvesterSolver(const SchurDecompZero &kdecomp, const SimilarityDecomp &fdecomp)
    : matrixK((zeroPad(kdecomp)) ?
              std::make_unique<QuasiTriangular>(kdecomp) :
              std::make_unique<QuasiTriangularZero>(kdecomp)),
      matrixF(std::make_unique<BlockDiagonal>(fdecomp.getB()))
  {
  }
  virtual ~SylvesterSolver() = default;
  virtual void solve(SylvParams &pars, KronVector &x) const = 0;
};

#endif /* SYLVESTER_SOLVER_H */
