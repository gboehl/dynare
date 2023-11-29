/*
 * Copyright © 2007-2022 Dynare Team
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

/*
 * This mex file computes A·(B⊗C) or A·(B⊗B) without explicitly building B⊗C or B⊗B, so that
 * one can consider large matrices A, B and/or C, and assuming that A is a the hessian of a DSGE
 * model (dynare format). This mex file should not be used outside dyn_second_order_solver.m.
 */

#include <algorithm>

#include <dynmex.h>

#include <omp.h>

#define DEBUG_OMP 0

void
sparse_hessian_times_B_kronecker_B(const mwIndex* isparseA, const mwIndex* jsparseA,
                                   const double* vsparseA, const double* B, double* D, size_t mA,
                                   size_t nA, size_t mB, size_t nB, int number_of_threads)
{
  /*
  **   Loop over the columns of B⊗B (or of the result matrix D).
  **   This loop is splitted into two nested loops because we use the
  **   symmetric pattern of the hessian matrix.
  */
#pragma omp parallel for num_threads(number_of_threads)
  for (mwIndex j1B = 0; j1B < static_cast<mwIndex>(nB); j1B++)
    {
#if DEBUG_OMP
      mexPrintf("%d thread number is %d (%d).\n", j1B, omp_get_thread_num(), omp_get_num_threads());
#endif
      for (mwIndex j2B = j1B; j2B < static_cast<mwIndex>(nB); j2B++)
        {
          mwIndex jj = j1B * nB + j2B; // column of B⊗B index.
          mwIndex iv = 0;
          int nz_in_column_ii_of_A = 0;
          mwIndex k1 = 0;
          mwIndex k2 = 0;
          /*
          ** Loop over the rows of B⊗B (column jj).
          */
          for (mwIndex ii = 0; ii < static_cast<mwIndex>(nA); ii++)
            {
              k1 = jsparseA[ii];
              k2 = jsparseA[ii + 1];
              if (k1 < k2) // otherwise column ii of A does not have non zero elements (and there is
                           // nothing to compute).
                {
                  ++nz_in_column_ii_of_A;
                  mwIndex i1B = ii / mB;
                  mwIndex i2B = ii % mB;
                  double bb = B[j1B * mB + i1B] * B[j2B * mB + i2B];
                  /*
                  ** Loop over the non zero entries of A(:,ii).
                  */
                  for (mwIndex k = k1; k < k2; k++)
                    {
                      mwIndex kk = isparseA[k];
                      D[jj * mA + kk] = D[jj * mA + kk] + bb * vsparseA[iv];
                      iv++;
                    }
                }
            }
          if (nz_in_column_ii_of_A > 0)
            std::copy_n(&D[jj * mA], mA, &D[(j2B * nB + j1B) * mA]);
        }
    }
}

void
sparse_hessian_times_B_kronecker_C(const mwIndex* isparseA, const mwIndex* jsparseA,
                                   const double* vsparseA, const double* B, const double* C,
                                   double* D, size_t mA, size_t nA, size_t mB, size_t nB, size_t mC,
                                   size_t nC, int number_of_threads)
{
  /*
  **   Loop over the columns of B⊗B (or of the result matrix D).
  */
#pragma omp parallel for num_threads(number_of_threads)
  for (mwIndex jj = 0; jj < static_cast<mwIndex>(nB * nC); jj++) // column of B⊗C index.
    {
      // Uncomment the following line to check if all processors are used.
#if DEBUG_OMP
      mexPrintf("%d thread number is %d (%d).\n", jj, omp_get_thread_num(), omp_get_num_threads());
#endif
      mwIndex jB = jj / nC;
      mwIndex jC = jj % nC;
      mwIndex k1 = 0;
      mwIndex k2 = 0;
      mwIndex iv = 0;
      int nz_in_column_ii_of_A = 0;
      /*
      ** Loop over the rows of B⊗C (column jj).
      */
      for (mwIndex ii = 0; ii < static_cast<mwIndex>(nA); ii++)
        {
          k1 = jsparseA[ii];
          k2 = jsparseA[ii + 1];
          if (k1 < k2) // otherwise column ii of A does not have non zero elements (and there is
                       // nothing to compute).
            {
              ++nz_in_column_ii_of_A;
              mwIndex iC = ii % mC;
              mwIndex iB = ii / mC;
              double cb = C[jC * mC + iC] * B[jB * mB + iB];
              /*
              ** Loop over the non zero entries of A(:,ii).
              */
              for (mwIndex k = k1; k < k2; k++)
                {
                  mwIndex kk = isparseA[k];
                  D[jj * mA + kk] = D[jj * mA + kk] + cb * vsparseA[iv];
                  iv++;
                }
            }
        }
    }
}

void
mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
  // Check input and output:
  if (nrhs > 4 || nrhs < 3 || nlhs != 1)
    {
      mexErrMsgTxt("sparse_hessian_times_B_kronecker_C takes 3 or 4 input arguments and provides 1 "
                   "output argument.");
      return; // Needed to shut up some GCC warnings
    }

  if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || !mxIsSparse(prhs[0]))
    mexErrMsgTxt("sparse_hessian_times_B_kronecker_C: First input must be a real sparse matrix.");
  if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || mxIsSparse(prhs[1]))
    mexErrMsgTxt("sparse_hessian_times_B_kronecker_C: Second input must be a real dense matrix.");
  if (nrhs == 4 && (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) || mxIsSparse(prhs[2])))
    mexErrMsgTxt("sparse_hessian_times_B_kronecker_C: Third input must be a real dense matrix.");

  // Get & Check dimensions (columns and rows):
  size_t mA = mxGetM(prhs[0]);
  size_t nA = mxGetN(prhs[0]);
  size_t mB = mxGetM(prhs[1]);
  size_t nB = mxGetN(prhs[1]);
  size_t mC, nC;
  if (nrhs == 4) // A·(B⊗C) is to be computed.
    {
      mC = mxGetM(prhs[2]);
      nC = mxGetN(prhs[2]);
      if (mB * mC != nA)
        mexErrMsgTxt("Input dimension error!");
    }
  else // A·(B⊗B) is to be computed.
    {
      if (mB * mB != nA)
        mexErrMsgTxt("Input dimension error!");
    }
  // Get input matrices:
  int numthreads;
  const double* B = mxGetPr(prhs[1]);
  const double* C;
  const mxArray* numthreads_mx;
  if (nrhs == 4)
    {
      C = mxGetPr(prhs[2]);
      numthreads_mx = prhs[3];
    }
  else
    numthreads_mx = prhs[2];

  if (!(mxIsScalar(numthreads_mx) && mxIsNumeric(numthreads_mx)))
    mexErrMsgTxt("sparse_hessian_times_B_kronecker_C: Last input must be a numeric scalar.");
  numthreads = static_cast<int>(mxGetScalar(numthreads_mx));
  if (numthreads <= 0)
    mexErrMsgTxt("sparse_hessian_times_B_kronecker_C: Last input must be a positive integer.");

  // Sparse (dynare) hessian matrix.
  const mwIndex* isparseA = mxGetIr(prhs[0]);
  const mwIndex* jsparseA = mxGetJc(prhs[0]);
  const double* vsparseA = mxGetPr(prhs[0]);
  // Initialization of the ouput:
  if (nrhs == 4)
    plhs[0] = mxCreateDoubleMatrix(mA, nB * nC, mxREAL);
  else
    plhs[0] = mxCreateDoubleMatrix(mA, nB * nB, mxREAL);
  double* D = mxGetPr(plhs[0]);

  // Computational part:
  if (nrhs == 3)
    sparse_hessian_times_B_kronecker_B(isparseA, jsparseA, vsparseA, B, D, mA, nA, mB, nB,
                                       numthreads);
  else
    sparse_hessian_times_B_kronecker_C(isparseA, jsparseA, vsparseA, B, C, D, mA, nA, mB, nB, mC,
                                       nC, numthreads);
}
