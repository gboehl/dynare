/*
 * Copyright © 2005 Ondra Kamenik
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

#include "vector_function.hh"

#include <dynlapack.h>

#include <cmath>
#include <algorithm>

/* Just an easy constructor of sequence of booleans defaulting to change
   everywhere. */

ParameterSignal::ParameterSignal(int n)
  : data(n, true)
{
}

/* This sets ‘false’ (no change) before a given parameter, and ‘true’ (change)
   after the given parameter (including). */

void
ParameterSignal::signalAfter(int l)
{
  for (size_t i = 0; i < std::min(static_cast<size_t>(l), data.size()); i++)
    data[i] = false;
  for (size_t i = l; i < data.size(); i++)
    data[i] = true;
}

/* This constructs a function set hardcopying also the first. */
VectorFunctionSet::VectorFunctionSet(const VectorFunction &f, int n)
  : funcs(n)
{
  for (int i = 0; i < n; i++)
    {
      func_copies.push_back(f.clone());
      funcs[i] = func_copies.back().get();
    }
}

/* This constructs a function set with shallow copy in the first and hard
   copies in others. */
VectorFunctionSet::VectorFunctionSet(VectorFunction &f, int n)
  : funcs(n)
{
  if (n > 0)
    funcs[0] = &f;
  for (int i = 1; i < n; i++)
    {
      func_copies.push_back(f.clone());
      funcs[i] = func_copies.back().get();
    }
}


/* Here we construct the object from the given function f and given
   variance-covariance matrix Σ=vcov. The matrix A is calculated as lower
   triangular and yields Σ=AAᵀ. */

GaussConverterFunction::GaussConverterFunction(VectorFunction &f, const GeneralMatrix &vcov)
  : VectorFunction(f), func(&f), A(vcov.nrows(), vcov.nrows()),
    multiplier(calcMultiplier())
{
  // TODO: raise if A.nrows() ≠ indim()
  calcCholeskyFactor(vcov);
}

GaussConverterFunction::GaussConverterFunction(std::unique_ptr<VectorFunction> f, const GeneralMatrix &vcov)
  : VectorFunction(*f), func_storage{move(f)}, func{func_storage.get()}, A(vcov.nrows(), vcov.nrows()),
    multiplier(calcMultiplier())
{
  // TODO: raise if A.nrows() ≠ indim()
  calcCholeskyFactor(vcov);
}

GaussConverterFunction::GaussConverterFunction(const GaussConverterFunction &f)
  : VectorFunction(f), func_storage{f.func->clone()}, func{func_storage.get()}, A(f.A),
    multiplier(f.multiplier)
{
}

/* Here we evaluate the function

    g(y) = 1/√(πⁿ) f(√2·Ay).

   Since the matrix A is lower triangular, the change signal for the function f
   will look like (0,…,0,1,…,1) where the first 1 is in the same position as
   the first change in the given signal ‘sig’ of the input y=point. */

void
GaussConverterFunction::eval(const Vector &point, const ParameterSignal &sig, Vector &out)
{
  ParameterSignal s(sig);
  int i = 0;
  while (i < indim() && !sig[i])
    i++;
  s.signalAfter(i);

  Vector x(indim());
  x.zeros();
  A.multaVec(x, point);
  x.mult(sqrt(2.0));

  func->eval(x, s, out);

  out.mult(multiplier);
}

/* This returns 1/√(πⁿ). */

double
GaussConverterFunction::calcMultiplier() const
{
  return sqrt(pow(M_PI, -1*indim()));
}

void
GaussConverterFunction::calcCholeskyFactor(const GeneralMatrix &vcov)
{
  A = vcov;

  lapack_int rows = A.nrows(), lda = A.getLD();
  for (int i = 0; i < rows; i++)
    for (int j = i+1; j < rows; j++)
      A.get(i, j) = 0.0;

  lapack_int info;
  dpotrf("L", &rows, A.base(), &lda, &info);
  // TODO: raise if info≠1
}
