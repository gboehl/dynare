/*
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

#include <string>
#include <memory>
#include <algorithm>

#include <dynmex.h>

#include "DynamicModelCaller.hh"

void
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  if (nlhs < 1 || nlhs > 2 || nrhs != 19)
    mexErrMsgTxt("Must have 19 input arguments and 1 or 2 output arguments");
  bool compute_jacobian = nlhs == 2;

  // Give explicit names to input arguments
  const mxArray *y_mx = prhs[0];
  const mxArray *basename_mx = prhs[1];
  const mxArray *ntt_mx = prhs[2];
  const mxArray *y0_mx = prhs[3];
  const mxArray *yT_mx = prhs[4];
  const mxArray *exo_path_mx = prhs[5];
  const mxArray *params_mx = prhs[6];
  const mxArray *steady_state_mx = prhs[7];
  const mxArray *periods_mx = prhs[8];
  const mxArray *ny_mx = prhs[9];
  const mxArray *maximum_lag_mx = prhs[10];
  const mxArray *maximum_endo_lag_mx = prhs[11];
  const mxArray *lead_lag_incidence_mx = prhs[12];
  const mxArray *nzij_pred_mx = prhs[13];
  const mxArray *nzij_current_mx = prhs[14];
  const mxArray *nzij_fwrd_mx = prhs[15];
  const mxArray *has_external_function_mx = prhs[16];
  const mxArray *use_dll_mx = prhs[17];
  const mxArray *num_threads_mx = prhs[18];

  // Check input and map it to local variables
  if (!(mxIsChar(basename_mx) && mxGetM(basename_mx) == 1))
    mexErrMsgTxt("basename should be a character string");
  char *basename = mxArrayToString(basename_mx);

  if (!(mxIsScalar(ntt_mx) && mxIsNumeric(ntt_mx)))
    mexErrMsgTxt("ntt should be a numeric scalar");
  size_t ntt = mxGetScalar(ntt_mx);

  if (!(mxIsScalar(periods_mx) && mxIsNumeric(periods_mx)))
    mexErrMsgTxt("periods should be a numeric scalar");
  mwIndex periods = static_cast<mwIndex>(mxGetScalar(periods_mx));

  if (!(mxIsScalar(ny_mx) && mxIsNumeric(ny_mx)))
    mexErrMsgTxt("ny should be a numeric scalar");
  mwIndex ny = static_cast<mwIndex>(mxGetScalar(ny_mx));

  if (!(mxIsScalar(maximum_lag_mx) && mxIsNumeric(maximum_lag_mx)))
    mexErrMsgTxt("maximum_lag should be a numeric scalar");
  mwIndex maximum_lag = static_cast<mwIndex>(mxGetScalar(maximum_lag_mx));

  if (!(mxIsScalar(maximum_endo_lag_mx) && mxIsNumeric(maximum_endo_lag_mx)))
    mexErrMsgTxt("maximum_endo_lag should be a numeric scalar");
  mwIndex maximum_endo_lag = static_cast<mwIndex>(mxGetScalar(maximum_endo_lag_mx));

  if (!(mxIsDouble(y_mx) && mxGetM(y_mx) == static_cast<size_t>(ny*periods) && mxGetN(y_mx) == 1))
    mexErrMsgTxt("y should be a double precision column-vector of ny*periods elements");
  const double *y = mxGetPr(y_mx);

  if (!(mxIsDouble(y0_mx) && mxGetM(y0_mx) == static_cast<size_t>(ny) && mxGetN(y0_mx) == 1))
    mexErrMsgTxt("y0 should be a double precision column-vector of ny elements");
  const double *y0 = mxGetPr(y0_mx);

  if (!(mxIsDouble(yT_mx) && mxGetM(yT_mx) == static_cast<size_t>(ny) && mxGetN(yT_mx) == 1))
    mexErrMsgTxt("yT should be a double precision column-vector of ny elements");
  const double *yT = mxGetPr(yT_mx);

  if (!(mxIsDouble(exo_path_mx) && mxGetM(exo_path_mx) >= static_cast<size_t>(periods+maximum_lag)))
    mexErrMsgTxt("exo_path should be a double precision matrix with at least periods+maximum_lag rows");
  mwIndex nx = static_cast<mwIndex>(mxGetN(exo_path_mx));
  size_t nb_row_x = mxGetM(exo_path_mx);
  const double *exo_path = mxGetPr(exo_path_mx);

  if (!(mxIsDouble(params_mx) && mxGetN(params_mx) == 1))
    mexErrMsgTxt("params should be a double precision column-vector");
  const double *params = mxGetPr(params_mx);

  if (!(mxIsDouble(steady_state_mx) && mxGetN(steady_state_mx) == 1))
    mexErrMsgTxt("steady_state should be a double precision column-vector");
  const double *steady_state = mxGetPr(steady_state_mx);

  if (!(mxIsDouble(lead_lag_incidence_mx) && mxGetM(lead_lag_incidence_mx) == static_cast<size_t>(2+maximum_endo_lag)
        && mxGetN(lead_lag_incidence_mx) == static_cast<size_t>(ny)))
    mexErrMsgTxt("lead_lag_incidence should be a double precision matrix with 2+maximum_endo_lag rows and endo_nbr columns");
  const double *lead_lag_incidence = mxGetPr(lead_lag_incidence_mx);

  if (!(mxIsInt32(nzij_pred_mx) && mxGetN(nzij_pred_mx) == 2))
    mexErrMsgTxt("nzij_pred should be an int32 matrix with 2 columns");
  size_t nnz_pred = mxGetM(nzij_pred_mx);
#if MX_HAS_INTERLEAVED_COMPLEX
  const int32_T *nzij_pred = mxGetInt32s(nzij_pred_mx);
#else
  const int32_T *nzij_pred = static_cast<const int32_T *>(mxGetData(nzij_pred_mx));
#endif

  if (!(mxIsInt32(nzij_current_mx) && mxGetN(nzij_current_mx) == 2))
    mexErrMsgTxt("nzij_current should be an int32 matrix with 2 columns");
  size_t nnz_current = mxGetM(nzij_current_mx);
#if MX_HAS_INTERLEAVED_COMPLEX
  const int32_T *nzij_current = mxGetInt32s(nzij_current_mx);
#else
  const int32_T *nzij_current = static_cast<const int32_T *>(mxGetData(nzij_current_mx));
#endif

  if (!(mxIsInt32(nzij_fwrd_mx) && mxGetN(nzij_fwrd_mx) == 2))
    mexErrMsgTxt("nzij_fwrd should be an int32 matrix with 2 columns");
  size_t nnz_fwrd = mxGetM(nzij_fwrd_mx);
#if MX_HAS_INTERLEAVED_COMPLEX
  const int32_T *nzij_fwrd = mxGetInt32s(nzij_fwrd_mx);
#else
  const int32_T *nzij_fwrd = static_cast<const int32_T *>(mxGetData(nzij_fwrd_mx));
#endif

  if (!(mxIsLogicalScalar(has_external_function_mx)))
    mexErrMsgTxt("has_external_function should be a logical scalar");
  bool has_external_function = static_cast<bool>(mxGetScalar(has_external_function_mx));

  if (!(mxIsLogicalScalar(use_dll_mx)))
    mexErrMsgTxt("use_dll should be a logical scalar");
  bool use_dll = static_cast<bool>(mxGetScalar(use_dll_mx));

  if (!(mxIsScalar(num_threads_mx) && mxIsNumeric(num_threads_mx)))
    mexErrMsgTxt("num_threads should be a numeric scalar");
  int num_threads = static_cast<int>(mxGetScalar(num_threads_mx));

  // Allocate output matrices
  plhs[0] = mxCreateDoubleMatrix(periods*ny, 1, mxREAL);
  double *stacked_residual = mxGetPr(plhs[0]);

  mwIndex nzmax = periods*nnz_current+(periods-1)*(nnz_pred+nnz_fwrd);

  double *stacked_jacobian = nullptr;
  mwIndex *ir = nullptr, *jc = nullptr;
  if (compute_jacobian)
    {
      plhs[1] = mxCreateSparse(periods*ny, periods*ny, nzmax, mxREAL);
      stacked_jacobian = mxGetPr(plhs[1]);
      ir = mxGetIr(plhs[1]);
      jc = mxGetJc(plhs[1]);

      /* Create the index vectors (IR, JC) of the sparse stacked jacobian. This
         makes parallelization across periods possible when evaluating the model,
         since all indices are known ex ante. */
      mwIndex k = 0;
      jc[0] = 0;
      for (mwIndex T = 0; T < periods; T++)
        {
          size_t row_pred = 0, row_current = 0, row_fwrd = 0;
          for (int32_T j = 0; j < static_cast<int32_T>(ny); j++)
            {
              if (T != 0)
                while (row_fwrd < nnz_fwrd && nzij_fwrd[row_fwrd+nnz_fwrd]-1 == j)
                  ir[k++] = (T-1)*ny + nzij_fwrd[row_fwrd++]-1;
              while (row_current < nnz_current && nzij_current[row_current+nnz_current]-1 == j)
                ir[k++] = T*ny + nzij_current[row_current++]-1;
              if (T != periods-1)
                while (row_pred < nnz_pred && nzij_pred[row_pred+nnz_pred]-1 == j)
                  ir[k++] = (T+1)*ny + nzij_pred[row_pred++]-1;
              jc[T*ny+j+1] = k;
            }
        }
    }

  size_t ndynvars = static_cast<size_t>(*std::max_element(lead_lag_incidence, lead_lag_incidence+(maximum_endo_lag+2)*ny));

  if (use_dll)
    DynamicModelDllCaller::load_dll(basename);

  DynamicModelCaller::error_msg.clear();

  /* Parallelize the main loop, if use_dll and no external function (to avoid
     parallel calls to MATLAB) */
#pragma omp parallel num_threads(num_threads) if (use_dll && !has_external_function)
  {
    // Allocate (thread-private) model evaluator (which allocates space for temporaries)
    std::unique_ptr<DynamicModelCaller> m;
    if (use_dll)
      m = std::make_unique<DynamicModelDllCaller>(ntt, nx, ny, ndynvars, exo_path, nb_row_x, params, steady_state, compute_jacobian);
    else
      m = std::make_unique<DynamicModelMatlabCaller>(basename, ntt, ndynvars, exo_path_mx, params_mx, steady_state_mx, compute_jacobian);

    // Main computing loop
#pragma omp for
    for (mwIndex T = 0; T < periods; T++)
      {
        // Fill vector of dynamic variables
        for (mwIndex j = 0; j < maximum_endo_lag+2; j++)
          for (mwIndex i = 0; i < ny; i++)
            {
              int idx = static_cast<int>(lead_lag_incidence[j+i*(2+maximum_endo_lag)])-1;
              if (idx != -1)
                {
                  if (T+j == maximum_endo_lag-1)
                    m->y(idx) = y0[i];
                  else if (T+j == maximum_endo_lag+periods)
                    m->y(idx) = yT[i];
                  else
                    m->y(idx) = y[i+(T+j-maximum_endo_lag)*ny];
                }
            }

        // Compute the residual and Jacobian, and fill the stacked residual
        m->eval(T+maximum_lag, stacked_residual+T*ny);

        if (compute_jacobian)
          {
            // Fill the stacked jacobian
            for (mwIndex col = T > maximum_endo_lag ? (T-maximum_endo_lag)*ny : 0; // We can't use std::max() here, because mwIndex is unsigned under MATLAB
                 col < std::min(periods*ny, (T+2)*ny); col++)
              {
                mwIndex k = jc[col];
                while (k < jc[col+1])
                  {
                    if (ir[k] < T*ny)
                      {
                        k++;
                        continue;
                      }
                    if (ir[k] >= (T+1)*ny)
                      break;

                    mwIndex eq = ir[k]-T*ny;
                    mwIndex lli_row = col/ny-(T-maximum_endo_lag); // 0, 1 or 2
                    mwIndex lli_col = col%ny;
                    mwIndex dynvar = static_cast<mwIndex>(lead_lag_incidence[lli_row+lli_col*(2+maximum_endo_lag)])-1;
                    stacked_jacobian[k] = m->jacobian(eq+dynvar*ny);
                    k++;
                  }
              }
          }
      }
  }

  /* Mimic a try/catch using a global string, since exceptions are not allowed
     to cross OpenMP boundary */
  if (!DynamicModelCaller::error_msg.empty())
    mexErrMsgTxt(DynamicModelCaller::error_msg.c_str());

  if (compute_jacobian)
    {
      /* The stacked jacobian so far constructed can still reference some zero
         elements, because some expressions that are symbolically different from
         zero may still evaluate to zero for some endogenous/parameter values. The
         following code further compresses the Jacobian by removing the references
         to those extra zeros. This makes a significant speed difference when
         inversing the Jacobian for some large models. */
      mwIndex k_orig = 0, k_new = 0;
      for (mwIndex col = 0; col < periods*ny; col++)
        {
          while (k_orig < jc[col+1])
            {
              if (stacked_jacobian[k_orig] != 0.0)
                {
                  if (k_new != k_orig)
                    {
                      stacked_jacobian[k_new] = stacked_jacobian[k_orig];
                      ir[k_new] = ir[k_orig];
                    }
                  k_new++;
                }
              k_orig++;
            }
          jc[col+1] = k_new;
        }
    }

  if (use_dll)
    DynamicModelDllCaller::unload_dll();
}
