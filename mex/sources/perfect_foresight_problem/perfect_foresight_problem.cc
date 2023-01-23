/*
 * Copyright © 2019-2023 Dynare Team
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

#include <string>
#include <memory>
#include <algorithm>

#include <dynmex.h>

#include "DynamicModelCaller.hh"

void
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  if (nlhs < 1 || nlhs > 2 || nrhs != 9)
    mexErrMsgTxt("Must have 9 input arguments and 1 or 2 output arguments");
  bool compute_jacobian = nlhs == 2;

  // Give explicit names to input arguments
  const mxArray *y_mx = prhs[0];
  const mxArray *y0_mx = prhs[1];
  const mxArray *yT_mx = prhs[2];
  const mxArray *exo_path_mx = prhs[3];
  const mxArray *params_mx = prhs[4];
  const mxArray *steady_state_mx = prhs[5];
  const mxArray *periods_mx = prhs[6];
  const mxArray *M_mx = prhs[7];
  const mxArray *options_mx = prhs[8];

  // Extract various fields from M_
  const mxArray *basename_mx = mxGetField(M_mx, 0, "fname");
  if (!(basename_mx && mxIsChar(basename_mx) && mxGetM(basename_mx) == 1))
    mexErrMsgTxt("M_.fname should be a character string");
  std::string basename{mxArrayToString(basename_mx)};

  const mxArray *endo_nbr_mx = mxGetField(M_mx, 0, "endo_nbr");
  if (!(endo_nbr_mx && mxIsScalar(endo_nbr_mx) && mxIsNumeric(endo_nbr_mx)))
    mexErrMsgTxt("M_.endo_nbr should be a numeric scalar");
  mwIndex ny = static_cast<mwIndex>(mxGetScalar(endo_nbr_mx));

  const mxArray *maximum_lag_mx = mxGetField(M_mx, 0, "maximum_lag");
  if (!(maximum_lag_mx && mxIsScalar(maximum_lag_mx) && mxIsNumeric(maximum_lag_mx)))
    mexErrMsgTxt("M_.maximum_lag should be a numeric scalar");
  mwIndex maximum_lag = static_cast<mwIndex>(mxGetScalar(maximum_lag_mx));

  const mxArray *dynamic_tmp_nbr_mx = mxGetField(M_mx, 0, "dynamic_tmp_nbr");
  if (!(dynamic_tmp_nbr_mx && mxIsDouble(dynamic_tmp_nbr_mx) && mxGetNumberOfElements(dynamic_tmp_nbr_mx) >= 2)
      || mxIsComplex(dynamic_tmp_nbr_mx) || mxIsSparse(dynamic_tmp_nbr_mx))
    mexErrMsgTxt("M_.dynamic_tmp_nbr should be a real dense array of at least 2 elements");
  size_t ntt {static_cast<size_t>(mxGetPr(dynamic_tmp_nbr_mx)[0]) +
    (compute_jacobian ? static_cast<size_t>(mxGetPr(dynamic_tmp_nbr_mx)[1]) : 0)};

  const mxArray *has_external_function_mx = mxGetField(M_mx, 0, "has_external_function");
  if (!(has_external_function_mx && mxIsLogicalScalar(has_external_function_mx)))
    mexErrMsgTxt("M_.has_external_function should be a logical scalar");
  bool has_external_function = static_cast<bool>(mxGetScalar(has_external_function_mx));

  // Extract various fields from options_
  const mxArray *use_dll_mx = mxGetField(options_mx, 0, "use_dll");
  if (!(use_dll_mx && mxIsLogicalScalar(use_dll_mx)))
    mexErrMsgTxt("options_.use_dll should be a logical scalar");
  bool use_dll = static_cast<bool>(mxGetScalar(use_dll_mx));

  const mxArray *linear_mx = mxGetField(options_mx, 0, "linear");
  if (!(linear_mx && mxIsLogicalScalar(linear_mx)))
    mexErrMsgTxt("options_.linear should be a logical scalar");
  bool linear = static_cast<bool>(mxGetScalar(linear_mx));

  const mxArray *threads_mx = mxGetField(options_mx, 0, "threads");
  if (!threads_mx)
    mexErrMsgTxt("Can't find field options_.threads");
  const mxArray *num_threads_mx = mxGetField(threads_mx, 0, "perfect_foresight_problem");
  if (!(num_threads_mx && mxIsScalar(num_threads_mx) && mxIsNumeric(num_threads_mx)))
    mexErrMsgTxt("options_.threads.perfect_foresight_problem should be a numeric scalar");
  int num_threads = static_cast<int>(mxGetScalar(num_threads_mx));

  // Check other input and map it to local variables
  if (!(mxIsScalar(periods_mx) && mxIsNumeric(periods_mx)))
    mexErrMsgTxt("periods should be a numeric scalar");
  mwIndex periods = static_cast<mwIndex>(mxGetScalar(periods_mx));

  if (!(mxIsDouble(y_mx) && mxGetM(y_mx) == static_cast<size_t>(ny*periods) && mxGetN(y_mx) == 1))
    mexErrMsgTxt("y should be a double precision column-vector of M_.endo_nbr*periods elements");
  const double *y = mxGetPr(y_mx);

  if (!(mxIsDouble(y0_mx) && mxGetM(y0_mx) == static_cast<size_t>(ny) && mxGetN(y0_mx) == 1))
    mexErrMsgTxt("y0 should be a double precision column-vector of M_.endo_nbr elements");
  const double *y0 = mxGetPr(y0_mx);

  if (!(mxIsDouble(yT_mx) && mxGetM(yT_mx) == static_cast<size_t>(ny) && mxGetN(yT_mx) == 1))
    mexErrMsgTxt("yT should be a double precision column-vector of M_.endo_nbr elements");
  const double *yT = mxGetPr(yT_mx);

  if (!(mxIsDouble(exo_path_mx) && mxGetM(exo_path_mx) >= static_cast<size_t>(periods+maximum_lag)))
    mexErrMsgTxt("exo_path should be a double precision matrix with at least periods+M_.maximum_lag rows");
  mwIndex nx = static_cast<mwIndex>(mxGetN(exo_path_mx));
  size_t nb_row_x = mxGetM(exo_path_mx);
  const double *exo_path = mxGetPr(exo_path_mx);

  const mxArray *g1_sparse_rowval_mx {mxGetField(M_mx, 0, "dynamic_g1_sparse_rowval")};
  if (!(mxIsInt32(g1_sparse_rowval_mx)))
    mexErrMsgTxt("M_.dynamic_g1_sparse_rowval should be an int32 vector");
#if MX_HAS_INTERLEAVED_COMPLEX
  const int32_T *g1_sparse_rowval {mxGetInt32s(g1_sparse_rowval_mx)};
#else
  const int32_T *g1_sparse_rowval {static_cast<const int32_T *>(mxGetData(g1_sparse_rowval_mx))};
#endif

  const mxArray *g1_sparse_colval_mx {mxGetField(M_mx, 0, "dynamic_g1_sparse_colval")};
  if (!(mxIsInt32(g1_sparse_colval_mx)))
    mexErrMsgTxt("M_.dynamic_g1_sparse_colval should be an int32 vector");
  if (mxGetNumberOfElements(g1_sparse_colval_mx) != mxGetNumberOfElements(g1_sparse_rowval_mx))
    mexErrMsgTxt("M_.dynamic_g1_sparse_colval should have the same length as M_.dynamic_g1_sparse_rowval");

  const mxArray *g1_sparse_colptr_mx {mxGetField(M_mx, 0, "dynamic_g1_sparse_colptr")};
  if (!(mxIsInt32(g1_sparse_colptr_mx) && mxGetNumberOfElements(g1_sparse_colptr_mx) != 3*static_cast<size_t>(ny)+1))
    mexErrMsgTxt(("M_.dynamic_g1_sparse_colptr should be an int32 vector with " + std::to_string(3*ny+1) + " elements").c_str());
#if MX_HAS_INTERLEAVED_COMPLEX
  const int32_T *g1_sparse_colptr {mxGetInt32s(g1_sparse_colptr_mx)};
#else
  const int32_T *g1_sparse_colptr {static_cast<const int32_T *>(mxGetData(g1_sparse_colptr_mx))};
#endif
  if (static_cast<size_t>(g1_sparse_colptr[3*ny+nx])-1 != mxGetNumberOfElements(g1_sparse_rowval_mx))
    mexErrMsgTxt("The size of M_.dynamic_g1_sparse_rowval is not consistent with the last element of M_.dynamic_g1_sparse_colptr");

  if (!(mxIsDouble(params_mx) && mxGetN(params_mx) == 1))
    mexErrMsgTxt("params should be a double precision column-vector");
  const double *params = mxGetPr(params_mx);

  if (!(mxIsDouble(steady_state_mx) && mxGetN(steady_state_mx) == 1))
    mexErrMsgTxt("steady_state should be a double precision column-vector");
  const double *steady_state = mxGetPr(steady_state_mx);

  // Allocate output matrices
  plhs[0] = mxCreateDoubleMatrix(periods*ny, 1, mxREAL);
  double *stacked_residual = mxGetPr(plhs[0]);

  /* Number of non-zero values in the stacked Jacobian.
     Contemporaneous derivatives appear at all periods, while lag or lead
     derivatives appear at all periods except one. */
  mwSize nzmax {static_cast<mwSize>((g1_sparse_colptr[3*ny]-1)*(periods-1) + (g1_sparse_colptr[2*ny]-g1_sparse_colptr[ny]))};

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
        for (mwIndex j {0}; j < ny; j++) // Column within the period (i.e. variable)
          {
            if (T > 0)
              for (int32_T idx {g1_sparse_colptr[2*ny+j]-1}; idx < g1_sparse_colptr[2*ny+j+1]-1; idx++)
                ir[k++] = (T-1)*ny + g1_sparse_rowval[idx]-1; // Derivatives w.r.t. y_{t+1} in T-1
            for (int32_T idx {g1_sparse_colptr[ny+j]-1}; idx < g1_sparse_colptr[ny+j+1]-1; idx++)
              ir[k++] = T*ny + g1_sparse_rowval[idx]-1; // Derivatives w.r.t. y_t in T
            if (T < periods-1)
              for (int32_T idx {g1_sparse_colptr[j]-1}; idx < g1_sparse_colptr[j+1]-1; idx++)
                ir[k++] = (T+1)*ny + g1_sparse_rowval[idx]-1; // Derivatives w.r.t. y_{t-1} in T+1
            jc[T*ny+j+1] = k;
          }
    }

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
      m = std::make_unique<DynamicModelDllCaller>(ntt, ny, nx, params, steady_state, g1_sparse_colptr, linear, compute_jacobian);
    else
      m = std::make_unique<DynamicModelMatlabCaller>(basename, ny, nx, params_mx, steady_state_mx, g1_sparse_rowval_mx, g1_sparse_colval_mx, g1_sparse_colptr_mx, linear, compute_jacobian);

    // Main computing loop
#pragma omp for
    for (mwIndex T = 0; T < periods; T++)
      {
        // Fill vector of dynamic variables
        if (T > 0 && T < periods-1)
          std::copy_n(y+(T-1)*ny, 3*ny, m->y());
        else if (T > 0) // Last simulation period
          {
            std::copy_n(y+(T-1)*ny, 2*ny, m->y());
            std::copy_n(yT, ny, m->y() + 2*ny);
          }
        else if (T < periods-1) // First simulation period
          {
            std::copy_n(y0, ny, m->y());
            std::copy_n(y+T*ny, 2*ny, m->y() + ny);
          }
        else // Special case: periods=1 (and so T=0)
          {
            std::copy_n(y0, ny, m->y());
            std::copy_n(y, ny, m->y() + ny);
            std::copy_n(yT, ny, m->y() + 2*ny);
          }

        // Fill exogenous
        for (mwIndex j {0}; j < nx; j++)
          m->x()[j] = exo_path[T+maximum_lag + nb_row_x*j];

        // Compute the residual and Jacobian, and fill the stacked residual
        m->eval(stacked_residual+T*ny);

        if (compute_jacobian)
          {
            // Fill the stacked jacobian
            for (mwIndex col { T > 1 ? (T-1)*ny : 0 }; // We can't use std::max() here, because mwIndex is unsigned under MATLAB
                 col < std::min(periods*ny, (T+2)*ny); col++)
              {
                mwIndex k = jc[col];
                while (k < jc[col+1])
                  {
                    if (ir[k] >= (T+1)*ny)
                      break; // Nothing to copy for this column
                    if (ir[k] >= T*ny)
                      {
                        /* Within the current column, this is the first line of
                           the stacked Jacobian that contains elements from the
                           (small) Jacobian just computed, so copy the whole
                           column of the latter to the former. */
                        m->copy_jacobian_column(col - (T-1)*ny, stacked_jacobian + k);
                        break;
                      }
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
