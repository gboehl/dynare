// Copyright (C) 2005-2011, Ondra Kamenik

#include "dynmex.h"
#include "mex.h"

#include "GeneralSylvester.hh"
#include "SylvException.hh"

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
    if (nhrs != 5 || nlhs > 3 || nlhs < 2)
      DYN_MEX_FUNC_ERR_MSG_TXT("Gensylv: Must have exactly 5 input args and either 2 or 3 output args.");

    auto order = static_cast<int>(mxGetScalar(prhs[0]));
    const mxArray *const A = prhs[1];
    const mxArray *const B = prhs[2];
    const mxArray *const C = prhs[3];
    const mxArray *const D = prhs[4];
    const mwSize *const Adims = mxGetDimensions(A);
    const mwSize *const Bdims = mxGetDimensions(B);
    const mwSize *const Cdims = mxGetDimensions(C);
    const mwSize *const Ddims = mxGetDimensions(D);

    if (Adims[0] != Adims[1])
      DYN_MEX_FUNC_ERR_MSG_TXT("Matrix A must be a square matrix.");
    if (Adims[0] != Bdims[0])
      DYN_MEX_FUNC_ERR_MSG_TXT("Matrix A and matrix B must have the same number of rows.");
    if (Adims[0] != Ddims[0])
      DYN_MEX_FUNC_ERR_MSG_TXT("Matrix A and matrix B must have the same number of rows.");
    if (Cdims[0] != Cdims[1])
      DYN_MEX_FUNC_ERR_MSG_TXT("Matrix C must be square.");
    if (Bdims[0] < Bdims[1])
      DYN_MEX_FUNC_ERR_MSG_TXT("Matrix B must not have more columns than rows.");
    if (Ddims[1] != static_cast<mwSize>(power(Cdims[0], order)))
      DYN_MEX_FUNC_ERR_MSG_TXT("Matrix D has wrong number of columns.");

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
        if (nlhs == 2)
          gen_sylv_solve(order, n, m, zero_cols, Avec, Bvec, Cvec, Xvec);
        else if (nlhs == 3)
          plhs[2] = gen_sylv_solve_and_check(order, n, m, zero_cols, Avec, Bvec, Cvec, Dvec, Xvec);
      }
    catch (const SylvException &e)
      {
        DYN_MEX_FUNC_ERR_MSG_TXT(e.getMessage().c_str());
      }
    plhs[1] = X;
    plhs[0] = mxCreateDoubleScalar(0);
  }
};
