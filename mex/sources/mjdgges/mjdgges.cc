/*
 * Copyright © 2006-2019 Dynare Team
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

#include <cmath>
#include <algorithm>
#include <memory>

#include <dynmex.h>
#include <dynlapack.h>

double criterium;

lapack_int
my_criteria(const double *alpha_r, const double *alpha_i, const double *beta)
{
  return *alpha_r * *alpha_r + *alpha_i * *alpha_i < criterium * criterium * *beta * *beta;
}

/* MATLAB interface */
void
mexFunction(int nlhs, mxArray *plhs[],
            int nrhs, const mxArray *prhs[])
{
  /* Check for proper number of arguments */
  if (nrhs < 2 || nrhs > 4 || nlhs == 0 || nlhs > 7)
    DYN_MEX_FUNC_ERR_MSG_TXT("MJDGGES: takes 2, 3 or 4 input arguments and between 1 and 7 output arguments.");

  /* Check that A and B are real matrices of the same dimension.*/
  size_t m1 = mxGetM(prhs[0]);
  size_t n1 = mxGetN(prhs[0]);
  size_t m2 = mxGetM(prhs[1]);
  size_t n2 = mxGetN(prhs[1]);
  if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])
      || !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])
      || m1 != n1 || m2 != n1 || m2 != n2)
    DYN_MEX_FUNC_ERR_MSG_TXT("MJDGGES requires two square real matrices of the same dimension.");

  /* Create a matrix for the return argument */
  plhs[1] = mxCreateDoubleMatrix(n1, n1, mxREAL);
  plhs[2] = mxCreateDoubleMatrix(n1, n1, mxREAL);
  plhs[3] = mxCreateDoubleMatrix(n1, n1, mxREAL);
  plhs[4] = mxCreateDoubleMatrix(1, 1, mxREAL);
  plhs[5] = mxCreateDoubleMatrix(n1, 1, mxCOMPLEX);
  plhs[6] = mxCreateDoubleMatrix(1, 1, mxREAL);

  /* Assign pointers to the various parameters */
  double *s = mxGetPr(plhs[1]);
  double *t = mxGetPr(plhs[2]);
  double *z = mxGetPr(plhs[3]);
  double *sdim = mxGetPr(plhs[4]);
#if defined(MATLAB_MEX_FILE) && MATLAB_VERSION >= 0x0904
  mxComplexDouble *gev = mxGetComplexDoubles(plhs[5]);
#else
  double *gev_r = mxGetPr(plhs[5]);
  double *gev_i = mxGetPi(plhs[5]);
#endif
  double *info = mxGetPr(plhs[6]);

  const double *a = mxGetPr(prhs[0]);
  const double *b = mxGetPr(prhs[1]);

  /* set criterium for stable eigenvalues */
  if (nrhs >= 3 && mxGetM(prhs[2]) > 0)
    criterium = *mxGetPr(prhs[2]);
  else
    criterium = 1+1e-6;

  /* set criterium for 0/0 generalized eigenvalues */
  double zhreshold;
  if (nrhs == 4 && mxGetM(prhs[3]) > 0)
    zhreshold = *mxGetPr(prhs[3]);
  else
    zhreshold = 1e-6;

  /* keep a and b intact */
  std::copy_n(a, n1*n1, s);
  std::copy_n(b, n1*n1, t);

  lapack_int i_n = static_cast<lapack_int>(n1);
  auto alpha_r = std::make_unique<double[]>(n1);
  auto alpha_i = std::make_unique<double[]>(n1);
  auto beta = std::make_unique<double[]>(n1);
  lapack_int lwork = 16*i_n+16;
  auto work = std::make_unique<double[]>(lwork);
  auto bwork = std::make_unique<lapack_int[]>(i_n);
  lapack_int i_info, i_sdim;

  dgges("N", "V", "S", my_criteria, &i_n, s, &i_n, t, &i_n, &i_sdim, alpha_r.get(), alpha_i.get(),
        beta.get(), nullptr, &i_n, z, &i_n, work.get(), &lwork, bwork.get(), &i_info);

  *sdim = static_cast<double>(i_sdim);
  *info = static_cast<double>(i_info);

  for (size_t i = 0; i < n1; i++)
    {
      if (std::abs(alpha_r[i]) > zhreshold || std::abs(beta[i]) > zhreshold)
#if defined(MATLAB_MEX_FILE) && MATLAB_VERSION >= 0x0904
        gev[i].real = alpha_r[i] / beta[i];
#else
        gev_r[i] = alpha_r[i] / beta[i];
#endif
      else
        {
          /* the ratio is too close to 0/0;
             returns specific error number only if no other error */
          if (*info == 0)
            *info = -30;
        }
      if (alpha_i[i] == 0.0 && beta[i] == 0.0)
#if defined(MATLAB_MEX_FILE) && MATLAB_VERSION >= 0x0904
        gev[i].imag = 0.0;
#else
        gev_i[i] = 0.0;
#endif
      else
#if defined(MATLAB_MEX_FILE) && MATLAB_VERSION >= 0x0904
        gev[i].imag = alpha_i[i] / beta[i];
#else
        gev_i[i] = alpha_i[i] / beta[i];
#endif
    }

  plhs[0] = mxCreateDoubleScalar(0);
}
