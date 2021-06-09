/*
 * Copyright © 2009-2020 Dynare Team.
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
 *
 * This mex file calls fortran routines from the Slicot library.
 */

/*
  ++    INPUTS
  ++    ======
  ++
  ++
  ++      [0]  T       (double)   n×n transition matrix.
  ++
  ++      [1]  QQ      (double)   n×n matrix (=R·Q·Rᵀ, where Q is the covariance matrix of the structural innovations).
  ++
  ++      [2]  Z       (double)   n×p selection matrix.
  ++
  ++      [3]  H       (double)   p×p covariance matrix of the measurement errors.
  ++
  ++
  ++
  ++
  ++    OUTPUTS
  ++    =======
  ++
  ++
  ++      [0]  P       (double)   n×n covariance matrix of the state vector.
  ++
  ++
  ++    NOTES
  ++    =====
  ++
  ++    [1] T = transpose(dynare transition matrix) and Z = transpose(dynare selection matrix).
*/

#include <algorithm>
#include <memory>
#include <dynmex.h>
#include <dynlapack.h>

#define sb02od FORTRAN_WRAPPER(sb02od)

extern "C"
{
  /* Note: matrices q, r and l may be modified internally (though they are
     restored on exit), hence their pointers are not declared as const */
  int sb02od(const char *dico, const char *jobb, const char *fact, const char *uplo,
             const char *jobl, const char *sort, const lapack_int *n, const lapack_int *m,
             const lapack_int *p, const double *a, const lapack_int *lda, const double *b,
             const lapack_int *ldb, double *q, const lapack_int *ldq, double *r, const lapack_int *ldr,
             double *l, const lapack_int *ldl, double *rcond, double *x, const lapack_int *ldx,
             double *alfar, double *alfai, double *beta, double *s, const lapack_int *lds, double *t,
             const lapack_int *ldt, double *u, const lapack_int *ldu, const double *tol,
             lapack_int *iwork, double *dwork, const lapack_int *ldwork, lapack_int *bwork,
             lapack_int *info);
}

void
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  // Check the number of arguments and set some flags.
  bool measurement_error_flag = true;
  if (nrhs < 3 || 4 < nrhs)
    mexErrMsgTxt("kalman_steady_state accepts either 3 or 4 input arguments!");

  if (nlhs != 1)
    mexErrMsgTxt("kalman_steady_state accepts exactly one output argument!");

  if (nrhs == 3)
    measurement_error_flag = false;

  // Check the type of the input arguments and get the size of the matrices.
  lapack_int n = mxGetM(prhs[0]);
  if (static_cast<size_t>(n) != mxGetN(prhs[0]))
    mexErrMsgTxt("kalman_steady_state: The first input argument (T) must be a square matrix!");
  if (mxIsNumeric(prhs[0]) == 0 || mxIsComplex(prhs[0]) == 1)
    mexErrMsgTxt("kalman_steady_state: The first input argument (T) must be a real matrix!");

  lapack_int q = mxGetM(prhs[1]);
  if (static_cast<size_t>(q) != mxGetN(prhs[1]))
    mexErrMsgTxt("kalman_steady_state: The second input argument (QQ) must be a square matrix!");
  if (mxIsNumeric(prhs[1]) == 0 || mxIsComplex(prhs[1]) == 1)
    mexErrMsgTxt("kalman_steady_state: The second input argument (QQ) must be a real matrix!");
  if (q != n)
    mexErrMsgTxt("kalman_steady_state: The size of the second input argument (QQ) must match the size of the first argument (T)!");
  lapack_int p = mxGetN(prhs[2]);
  if (mxGetM(prhs[2]) != static_cast<size_t>(n))
    mexErrMsgTxt("kalman_steady_state: The number of rows of the third argument (Z) must match the number of rows of the first argument (T)!");
  if (mxIsNumeric(prhs[2]) == 0 || mxIsComplex(prhs[2]) == 1)
    mexErrMsgTxt("kalman_steady_state: The third input argument (Z) must be a real matrix!");
  if (measurement_error_flag)
    {
      if (mxGetM(prhs[3]) != mxGetN(prhs[3]))
        mexErrMsgTxt("kalman_steady_state: The fourth input argument (H) must be a square matrix!");
      if (mxGetM(prhs[3]) != static_cast<size_t>(p))
        mexErrMsgTxt("kalman_steady_state: The number of rows of the fourth input argument (H) must match the number of rows of the third input argument!");
      if (mxIsNumeric(prhs[3]) == 0 || mxIsComplex(prhs[3]) == 1)
        mexErrMsgTxt("kalman_steady_state: The fifth input argument (H) must be a real matrix!");
    }
  // Get input matrices.
  const double *T = mxGetPr(prhs[0]);
  auto QQ = std::make_unique<double[]>(n*n);
  std::copy_n(mxGetPr(prhs[1]), n*n, QQ.get());
  const double *Z = mxGetPr(prhs[2]);
  auto H = std::make_unique<double[]>(p*p);
  if (measurement_error_flag)
    std::copy_n(mxGetPr(prhs[3]), p*p, H.get());
  // L will not be used.
  auto L = std::make_unique<double[]>(n*p);
  lapack_int nn = 2*n;
  lapack_int LDA = std::max(static_cast<lapack_int>(1), n);
  lapack_int LDQ = LDA;
  lapack_int LDU = std::max(static_cast<lapack_int>(1), nn);
  lapack_int LDS = std::max(static_cast<lapack_int>(1), nn+p);
  lapack_int LIWORK = std::max(static_cast<lapack_int>(1), std::max(p, nn));
  lapack_int LDR = std::max(static_cast<lapack_int>(1), p);
  lapack_int LDB = LDA, LDL = LDA, LDT = LDS, LDX = LDA;
  lapack_int LDWORK = std::max(static_cast<lapack_int>(7)*(static_cast<lapack_int>(2)*n + static_cast<lapack_int>(1)) + static_cast<lapack_int>(16), static_cast<lapack_int>(16)*n);
  LDWORK = std::max(LDWORK, std::max(static_cast<lapack_int>(2)*n + p, static_cast<lapack_int>(3)*p));
  double tolerance = 1e-16;
  lapack_int INFO;
  // Outputs of subroutine sb02OD
  double rcond;
  auto WR = std::make_unique<double[]>(nn);
  auto WI = std::make_unique<double[]>(nn);
  auto BETA = std::make_unique<double[]>(nn);
  auto S = std::make_unique<double[]>(LDS*(nn+p));
  auto TT = std::make_unique<double[]>(LDT*nn);
  auto UU = std::make_unique<double[]>(LDU*nn);
  // Working arrays
  auto IWORK = std::make_unique<lapack_int[]>(LIWORK);
  auto DWORK = std::make_unique<double[]>(LDWORK);
  auto BWORK = std::make_unique<lapack_int[]>(nn);
  // Initialize the output of the mex file
  plhs[0] = mxCreateDoubleMatrix(n, n, mxREAL);
  double *P = mxGetPr(plhs[0]);
  // Call the slicot routine
  sb02od("D", // We want to solve a discrete Riccati equation.
         "B", // Matrices Z and H are given.
         "N", // Given matrices H and QQ are not factored.
         "U", // Upper triangle of matrix H is stored.
         "Z", // L matrix is zero.
         "S", // Stable eigenvalues come first.
         &n, &p, &p, T, &LDA, Z, &LDB, QQ.get(), &LDQ, H.get(), &LDR, L.get(), &LDL,
         &rcond, P, &LDX, WR.get(), WI.get(), BETA.get(), S.get(), &LDS, TT.get(), &LDT, UU.get(),
         &LDU, &tolerance, IWORK.get(), DWORK.get(), &LDWORK, BWORK.get(), &INFO);
  switch (INFO)
    {
    case 0:
      break;
    case 1:
      mexErrMsgTxt("The computed extended matrix pencil is singular, possibly due to rounding errors.");
      break;
    case 2:
      mexErrMsgTxt("The QZ (or QR) algorithm failed!");
      break;
    case 3:
      mexErrMsgTxt("The reordering of the (generalized) eigenvalues failed!");
      break;
    case 4:
      mexErrMsgTxt("After reordering, roundoff changed values of some complex eigenvalues so that leading eigenvalues\n in the (generalized) Schur form no longer satisfy the stability condition; this could also be caused due to scaling.");
      break;
    case 5:
      mexErrMsgTxt("The computed dimension of the solution does not equal n!");
      break;
    case 6:
      mexErrMsgTxt("A singular matrix was encountered during the computation of the solution matrix P!");
      break;
    default:
      mexErrMsgTxt("Unknown problem!");
    }
}
