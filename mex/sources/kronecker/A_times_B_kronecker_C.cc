/*
 * Copyright © 2007-2020 Dynare Team
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

/*
 * This mex file computes A·(B⊗C) or A·(B⊗B) without explicitly building B⊗C or B⊗B, so that
 * one can consider large matrices B and/or C.
 */

#include <dynmex.h>
#include <dynblas.h>

void
full_A_times_kronecker_B_C(const double *A, const double *B, const double *C, double *D,
                           blas_int mA, blas_int nA, blas_int mB, blas_int nB, blas_int mC, blas_int nC)
{
  const blas_int shiftA = mA*mC;
  const blas_int shiftD = mA*nC;
  blas_int kd = 0, ka = 0;
  double one = 1.0;
  for (blas_int col = 0; col < nB; col++)
    {
      ka = 0;
      for (blas_int row = 0; row < mB; row++)
        {
          dgemm("N", "N", &mA, &nC, &mC, &B[mB*col+row], &A[ka], &mA, C, &mC, &one, &D[kd], &mA);
          ka += shiftA;
        }
      kd += shiftD;
    }
}

void
full_A_times_kronecker_B_B(const double *A, const double *B, double *D, blas_int mA, blas_int nA, blas_int mB, blas_int nB)
{
  const blas_int shiftA = mA*mB;
  const blas_int shiftD = mA*nB;
  blas_int kd = 0, ka = 0;
  double one = 1.0;
  for (blas_int col = 0; col < nB; col++)
    {
      ka = 0;
      for (blas_int row = 0; row < mB; row++)
        {
          dgemm("N", "N", &mA, &nB, &mB, &B[mB*col+row], &A[ka], &mA, B, &mB, &one, &D[kd], &mA);
          ka += shiftA;
        }
      kd += shiftD;
    }
}

void
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  // Check input and output:
  if (nrhs > 3 || nrhs < 2 || nlhs != 1)
    {
      mexErrMsgTxt("A_times_B_kronecker_C takes 2 or 3 input arguments and provides 1 output argument.");
      return; // Needed to shut up some GCC warnings
    }

  // Get & Check dimensions (columns and rows):
  size_t mA = mxGetM(prhs[0]);
  size_t nA = mxGetN(prhs[0]);
  size_t mB = mxGetM(prhs[1]);
  size_t nB = mxGetN(prhs[1]);
  size_t mC, nC;
  if (nrhs == 3) // A·(B⊗C) is to be computed.
    {
      mC = mxGetM(prhs[2]);
      nC = mxGetN(prhs[2]);
      if (mB*mC != nA)
        mexErrMsgTxt("Input dimension error!");
    }
  else // A·(B⊗B) is to be computed.
    {
      if (mB*mB != nA)
        mexErrMsgTxt("Input dimension error!");
    }
  // Get input matrices:
  const double *A = mxGetPr(prhs[0]);
  const double *B = mxGetPr(prhs[1]);
  const double *C{nullptr};
  if (nrhs == 3)
    C = mxGetPr(prhs[2]);

  // Initialization of the ouput:
  if (nrhs == 3)
    plhs[0] = mxCreateDoubleMatrix(mA, nB*nC, mxREAL);
  else
    plhs[0] = mxCreateDoubleMatrix(mA, nB*nB, mxREAL);
  double *D = mxGetPr(plhs[0]);

  // Computational part:
  if (nrhs == 2)
    full_A_times_kronecker_B_B(A, B, D, mA, nA, mB, nB);
  else
    full_A_times_kronecker_B_C(A, B, C, D, mA, nA, mB, nB, mC, nC);
}
