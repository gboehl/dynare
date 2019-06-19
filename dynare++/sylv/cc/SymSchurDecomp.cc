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
 * along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "SymSchurDecomp.hh"
#include "SylvException.hh"

#include <dynlapack.h>

#include <algorithm>
#include <cmath>
#include <memory>

SymSchurDecomp::SymSchurDecomp(const ConstGeneralMatrix &mata)
  : lambda(mata.nrows()), q(mata.nrows())
{
  // check mata is square
  if (mata.nrows() != mata.ncols())
    throw SYLV_MES_EXCEPTION("Matrix is not square in SymSchurDecomp constructor");

  // prepare for dsyevr
  lapack_int n = mata.nrows();
  GeneralMatrix tmpa(mata);
  double *a = tmpa.base();
  lapack_int lda = tmpa.getLD();
  double dum;
  double *vl = &dum;
  double *vu = &dum;
  lapack_int idum;
  lapack_int *il = &idum;
  lapack_int *iu = &idum;
  double abstol = 0.0;
  lapack_int m = n;
  double *w = lambda.base();
  double *z = q.base();
  lapack_int ldz = q.getLD();
  auto isuppz = std::make_unique<lapack_int[]>(2*std::max(1, static_cast<int>(m)));
  double tmpwork;
  lapack_int lwork = -1;
  lapack_int tmpiwork;
  lapack_int liwork = -1;
  lapack_int info;

  // query for lwork and liwork
  dsyevr("V", "A", "U", &n, a, &lda, vl, vu, il, iu, &abstol,
         &m, w, z, &ldz, isuppz.get(), &tmpwork, &lwork, &tmpiwork, &liwork, &info);
  lwork = static_cast<lapack_int>(tmpwork);
  liwork = tmpiwork;
  // allocate work arrays
  auto work = std::make_unique<double[]>(lwork);
  auto iwork = std::make_unique<lapack_int[]>(liwork);

  // do the calculation
  dsyevr("V", "A", "U", &n, a, &lda, vl, vu, il, iu, &abstol,
         &m, w, z, &ldz, isuppz.get(), work.get(), &lwork, iwork.get(), &liwork, &info);

  if (info < 0)
    throw SYLV_MES_EXCEPTION("Internal error in SymSchurDecomp constructor");
  if (info > 0)
    throw SYLV_MES_EXCEPTION("Internal LAPACK error in DSYEVR");
}

void
SymSchurDecomp::getFactor(GeneralMatrix &f) const
{
  if (f.nrows() != q.nrows())
    throw SYLV_MES_EXCEPTION("Wrong dimension of factor matrix in SymSchurDecomp::getFactor");
  if (f.nrows() != f.ncols())
    throw SYLV_MES_EXCEPTION("Factor matrix is not square in SymSchurDecomp::getFactor");
  if (!isPositiveSemidefinite())
    throw SYLV_MES_EXCEPTION("Symmetric decomposition not positive semidefinite in SymSchurDecomp::getFactor");

  f = q;
  for (int i = 0; i < f.ncols(); i++)
    {
      Vector fi{f.getCol(i)};
      fi.mult(std::sqrt(lambda[i]));
    }
}

/* LAPACK says that eigenvalues are ordered in ascending order, but we
   do not rely on it */
bool
SymSchurDecomp::isPositiveSemidefinite() const
{
  for (int i = 0; i < lambda.length(); i++)
    if (lambda[i] < 0)
      return false;
  return true;
}

void
SymSchurDecomp::correctDefinitness(double tol)
{
  for (int i = 0; i < lambda.length(); i++)
    if (lambda[i] < 0 && lambda[i] > -tol)
      lambda[i] = 0.0;
}
