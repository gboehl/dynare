/*
 * Copyright © 2005-2011 Ondra Kamenik
 * Copyright © 2019-2022 Dynare Team
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

#include "dynmex.h"
#include "mex.h"

#include "GeneralSylvester.hh"
#include "SylvException.hh"
#include "int_power.hh"

void
gen_sylv_solve(int order, int n, int m, int zero_cols,
               const ConstVector &A, const ConstVector &B,
               const ConstVector &C, Vector &X)
{
  GeneralSylvester s(order, n, m, zero_cols, A, B, C, X, false);
  s.solve();
}

mxArray *
gen_sylv_solve_and_check(int order, int n, int m, int zero_cols,
                         const ConstVector &A, const ConstVector &B,
                         const ConstVector &C, const ConstVector &D,
                         Vector &X)
{
  GeneralSylvester s(order, n, m, zero_cols, A, B, C, X, true);
  s.solve();
  s.check(D);
  return s.getParams().createStructArray();
}

extern "C" {
  void
  mexFunction(int nlhs, mxArray *plhs[],
              int nhrs, const mxArray *prhs[])
  {
    if (nhrs != 5 || nlhs > 2 || nlhs < 1)
      mexErrMsgTxt("Gensylv: Must have exactly 5 input args and either 1 or 2 output args.");

    if (!(mxIsScalar(prhs[0]) && mxIsNumeric(prhs[0])))
      mexErrMsgTxt("First argument should be a numeric scalar");

    auto order = static_cast<int>(mxGetScalar(prhs[0]));
    const mxArray *const A = prhs[1];
    const mxArray *const B = prhs[2];
    const mxArray *const C = prhs[3];
    const mxArray *const D = prhs[4];

    if (!mxIsDouble(A) || mxIsComplex(A) || mxIsSparse(A))
      mexErrMsgTxt("Matrix A must be a real dense matrix.");
    if (!mxIsDouble(B) || mxIsComplex(B) || mxIsSparse(B))
      mexErrMsgTxt("Matrix B must be a real dense matrix.");
    if (!mxIsDouble(C) || mxIsComplex(C) || mxIsSparse(C))
      mexErrMsgTxt("Matrix C must be a real dense matrix.");
    if (!mxIsDouble(D) || mxIsComplex(D) || mxIsSparse(D))
      mexErrMsgTxt("Matrix D must be a real dense matrix.");

    const mwSize *const Adims = mxGetDimensions(A);
    const mwSize *const Bdims = mxGetDimensions(B);
    const mwSize *const Cdims = mxGetDimensions(C);
    const mwSize *const Ddims = mxGetDimensions(D);

    if (Adims[0] != Adims[1])
      mexErrMsgTxt("Matrix A must be a square matrix.");
    if (Adims[0] != Bdims[0])
      mexErrMsgTxt("Matrix A and matrix B must have the same number of rows.");
    if (Adims[0] != Ddims[0])
      mexErrMsgTxt("Matrix A and matrix B must have the same number of rows.");
    if (Cdims[0] != Cdims[1])
      mexErrMsgTxt("Matrix C must be square.");
    if (Bdims[0] < Bdims[1])
      mexErrMsgTxt("Matrix B must not have more columns than rows.");
    if (Ddims[1] != static_cast<mwSize>(power(Cdims[0], order)))
      mexErrMsgTxt("Matrix D has wrong number of columns.");

    auto n = static_cast<int>(Adims[0]);
    auto m = static_cast<int>(Cdims[0]);
    auto zero_cols = static_cast<int>(Bdims[0] - Bdims[1]);
    mxArray *X = mxCreateDoubleMatrix(Ddims[0], Ddims[1], mxREAL);
    // copy D to X
    ConstVector Avec{A}, Bvec{B}, Cvec{C}, Dvec{D};
    Vector Xvec{X};
    Xvec = Dvec;
    // solve (or solve and check)
    try
      {
        if (nlhs == 1)
          gen_sylv_solve(order, n, m, zero_cols, Avec, Bvec, Cvec, Xvec);
        else if (nlhs == 2)
          plhs[1] = gen_sylv_solve_and_check(order, n, m, zero_cols, Avec, Bvec, Cvec, Dvec, Xvec);
      }
    catch (const SylvException &e)
      {
        mexErrMsgTxt(e.getMessage().c_str());
      }
    plhs[0] = X;
  }
};
