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

#include "DynamicModelCaller.hh"

#include <algorithm>
#include <filesystem>

using namespace std::literals::string_literals;

std::string DynamicModelCaller::error_msg;

#if !defined(_WIN32) && !defined(__CYGWIN32__)
void* DynamicModelDllCaller::resid_mex {nullptr};
void* DynamicModelDllCaller::g1_mex {nullptr};
#else
HINSTANCE DynamicModelDllCaller::resid_mex {nullptr};
HINSTANCE DynamicModelDllCaller::g1_mex {nullptr};
#endif
DynamicModelDllCaller::dynamic_tt_fct DynamicModelDllCaller::residual_tt_fct {nullptr},
    DynamicModelDllCaller::g1_tt_fct {nullptr};
DynamicModelDllCaller::dynamic_fct DynamicModelDllCaller::residual_fct {nullptr},
    DynamicModelDllCaller::g1_fct {nullptr};

void
DynamicModelDllCaller::load_dll(const std::string& basename)
{
  // Load symbols from dynamic MEX
  const std::filesystem::path sparse_dir {"+" + basename + "/+sparse/"};
  const std::filesystem::path resid_mex_name {sparse_dir / ("dynamic_resid"s + MEXEXT)},
      g1_mex_name {sparse_dir / ("dynamic_g1"s + MEXEXT)};
#if !defined(__CYGWIN32__) && !defined(_WIN32)
  resid_mex = dlopen(resid_mex_name.c_str(), RTLD_NOW);
  g1_mex = dlopen(g1_mex_name.c_str(), RTLD_NOW);
#else
  resid_mex = LoadLibraryW(resid_mex_name.c_str());
  g1_mex = LoadLibraryW(g1_mex_name.c_str());
#endif
  if (!resid_mex)
    mexErrMsgTxt("Can't load dynamic_resid MEX file");
  if (!g1_mex)
    mexErrMsgTxt("Can't load dynamic_g1 MEX file");

#if !defined(__CYGWIN32__) && !defined(_WIN32)
  residual_tt_fct = reinterpret_cast<dynamic_tt_fct>(dlsym(resid_mex, "dynamic_resid_tt"));
  residual_fct = reinterpret_cast<dynamic_fct>(dlsym(resid_mex, "dynamic_resid"));
  g1_tt_fct = reinterpret_cast<dynamic_tt_fct>(dlsym(g1_mex, "dynamic_g1_tt"));
  g1_fct = reinterpret_cast<dynamic_fct>(dlsym(g1_mex, "dynamic_g1"));
#else
# pragma GCC diagnostic push
# pragma GCC diagnostic ignored "-Wcast-function-type"
  residual_tt_fct = reinterpret_cast<dynamic_tt_fct>(GetProcAddress(resid_mex, "dynamic_resid_tt"));
  residual_fct = reinterpret_cast<dynamic_fct>(GetProcAddress(resid_mex, "dynamic_resid"));
  g1_tt_fct = reinterpret_cast<dynamic_tt_fct>(GetProcAddress(g1_mex, "dynamic_g1_tt"));
  g1_fct = reinterpret_cast<dynamic_fct>(GetProcAddress(g1_mex, "dynamic_g1"));
# pragma GCC diagnostic pop
#endif
  if (!residual_tt_fct || !residual_fct || !g1_tt_fct || !g1_fct)
    mexErrMsgTxt("Can't load functions in dynamic MEX file");
}

void
DynamicModelDllCaller::unload_dll()
{
#if !defined(__CYGWIN32__) && !defined(_WIN32)
  dlclose(resid_mex);
  dlclose(g1_mex);
#else
  FreeLibrary(resid_mex);
  FreeLibrary(g1_mex);
#endif
}

DynamicModelDllCaller::DynamicModelDllCaller(size_t ntt, mwIndex ny, mwIndex nx,
                                             const double* params_arg,
                                             const double* steady_state_arg,
                                             const int32_T* g1_sparse_colptr_arg, bool linear_arg,
                                             bool compute_jacobian_arg) :
    DynamicModelCaller {linear_arg, compute_jacobian_arg},
    params {params_arg},
    steady_state {steady_state_arg},
    g1_sparse_colptr {g1_sparse_colptr_arg}
{
  tt = std::make_unique<double[]>(ntt);
  y_p = std::make_unique<double[]>(3 * ny);
  x_p = std::make_unique<double[]>(nx);
  if (compute_jacobian)
    jacobian_p = std::make_unique<double[]>(g1_sparse_colptr[3 * ny + nx] - 1);
}

void
DynamicModelDllCaller::copy_jacobian_column(mwIndex col, double* dest) const
{
  std::copy_n(jacobian_p.get() + g1_sparse_colptr[col] - 1,
              g1_sparse_colptr[col + 1] - g1_sparse_colptr[col], dest);
}

void
DynamicModelDllCaller::eval(double* resid)
{
  residual_tt_fct(y_p.get(), x_p.get(), params, steady_state, tt.get());
  residual_fct(y_p.get(), x_p.get(), params, steady_state, tt.get(), resid);
  if (compute_jacobian)
    {
      g1_tt_fct(y_p.get(), x_p.get(), params, steady_state, tt.get());
      g1_fct(y_p.get(), x_p.get(), params, steady_state, tt.get(), jacobian_p.get());

      if (linear)
        compute_jacobian = false; // If model is linear, no need to recompute Jacobian later
    }
}

DynamicModelMatlabCaller::DynamicModelMatlabCaller(std::string basename_arg, mwIndex ny, mwIndex nx,
                                                   const mxArray* params_mx_arg,
                                                   const mxArray* steady_state_mx_arg,
                                                   const mxArray* g1_sparse_rowval_mx_arg,
                                                   const mxArray* g1_sparse_colval_mx_arg,
                                                   const mxArray* g1_sparse_colptr_mx_arg,
                                                   bool linear_arg, bool compute_jacobian_arg) :
    DynamicModelCaller {linear_arg, compute_jacobian_arg},
    basename {std::move(basename_arg)},
    y_mx {mxCreateDoubleMatrix(3 * ny, 1, mxREAL)},
    x_mx {mxCreateDoubleMatrix(nx, 1, mxREAL)},
    jacobian_mx {nullptr},
    params_mx {mxDuplicateArray(params_mx_arg)},
    steady_state_mx {mxDuplicateArray(steady_state_mx_arg)},
    g1_sparse_rowval_mx {mxDuplicateArray(g1_sparse_rowval_mx_arg)},
    g1_sparse_colval_mx {mxDuplicateArray(g1_sparse_colval_mx_arg)},
    g1_sparse_colptr_mx {mxDuplicateArray(g1_sparse_colptr_mx_arg)}
{
}

DynamicModelMatlabCaller::~DynamicModelMatlabCaller()
{
  mxDestroyArray(y_mx);
  mxDestroyArray(x_mx);
  if (jacobian_mx)
    mxDestroyArray(jacobian_mx);
  mxDestroyArray(params_mx);
  mxDestroyArray(steady_state_mx);
  mxDestroyArray(g1_sparse_rowval_mx);
  mxDestroyArray(g1_sparse_colval_mx);
  mxDestroyArray(g1_sparse_colptr_mx);
}

void
DynamicModelMatlabCaller::copy_jacobian_column(mwIndex col, double* dest) const
{
  if (jacobian_mx)
    {
#if MX_HAS_INTERLEAVED_COMPLEX
      const int32_T* g1_sparse_rowval {mxGetInt32s(g1_sparse_rowval_mx)};
      const int32_T* g1_sparse_colptr {mxGetInt32s(g1_sparse_colptr_mx)};
#else
      const int32_T* g1_sparse_rowval {static_cast<const int32_T*>(mxGetData(g1_sparse_rowval_mx))};
      const int32_T* g1_sparse_colptr {static_cast<const int32_T*>(mxGetData(g1_sparse_colptr_mx))};
#endif

      /* We cannot assume that jacobian_mx internally uses
         g1_sparse_{rowval,colval,colptr}, because the call to sparse() in
         dynamic_g1.m may have further compressed the matrix by removing
         elements that are numerically zero, despite being symbolically
         non-zero. */
      mwIndex *ir {mxGetIr(jacobian_mx)}, *jc {mxGetJc(jacobian_mx)};
      mwIndex isrc {jc[col]}; // Index in value array of source Jacobian
      for (mwIndex idest {0}; // Index in value array of destination Jacobian
           idest < static_cast<mwIndex>(g1_sparse_colptr[col + 1] - g1_sparse_colptr[col]); idest++)
        {
          mwIndex row {
              static_cast<mwIndex>(g1_sparse_rowval[idest + g1_sparse_colptr[col] - 1] - 1)};
          while (isrc < jc[col + 1] && ir[isrc] < row)
            isrc++;
          if (isrc < jc[col + 1] && ir[isrc] == row)
            dest[idest] = mxGetPr(jacobian_mx)[isrc];
          else
            dest[idest] = 0.0;
        }
    }
}

void
DynamicModelMatlabCaller::eval(double* resid)
{
  mxArray *T_order_mx, *T_mx;

  {
    // Compute residuals
    std::string funcname {basename + ".sparse.dynamic_resid"};
    mxArray *plhs[3], *prhs[] = {y_mx, x_mx, params_mx, steady_state_mx};

    mxArray* exception {mexCallMATLABWithTrap(std::extent_v<decltype(plhs)>, plhs,
                                              std::extent_v<decltype(prhs)>, prhs,
                                              funcname.c_str())};
    if (exception)
      {
        error_msg = "An error occurred when calling " + funcname;
        return; // Avoid manipulating null pointers in plhs, see #1832
      }

    if (!mxIsDouble(plhs[0]) || mxIsSparse(plhs[0]))
      {
        error_msg = "Residuals should be a dense array of double floats";
        return;
      }

    if (mxIsComplex(plhs[0]))
      plhs[0] = cmplxToReal<false>(plhs[0]);

    std::copy_n(mxGetPr(plhs[0]), mxGetNumberOfElements(plhs[0]), resid);
    mxDestroyArray(plhs[0]);

    T_order_mx = plhs[1];
    T_mx = plhs[2];
  }

  if (compute_jacobian)
    {
      // Compute Jacobian
      std::string funcname {basename + ".sparse.dynamic_g1"};
      mxArray *plhs[1], *prhs[] = {y_mx,
                                   x_mx,
                                   params_mx,
                                   steady_state_mx,
                                   g1_sparse_rowval_mx,
                                   g1_sparse_colval_mx,
                                   g1_sparse_colptr_mx,
                                   T_order_mx,
                                   T_mx};

      mxArray* exception {mexCallMATLABWithTrap(std::extent_v<decltype(plhs)>, plhs,
                                                std::extent_v<decltype(prhs)>, prhs,
                                                funcname.c_str())};
      if (exception)
        {
          error_msg = "An error occurred when calling " + funcname;
          return; // Avoid manipulating null pointers in plhs, see #1832
        }

      if (jacobian_mx)
        {
          mxDestroyArray(jacobian_mx);
          jacobian_mx = nullptr;
        }

      if (!mxIsDouble(plhs[0]) || !mxIsSparse(plhs[0]))
        {
          error_msg = "Jacobian should be a sparse array of double floats";
          return;
        }

      if (mxIsComplex(plhs[0]))
        jacobian_mx = cmplxToReal<true>(plhs[0]);
      else
        jacobian_mx = plhs[0];

      if (linear)
        compute_jacobian = false; // If model is linear, no need to recompute Jacobian later
    }

  mxDestroyArray(T_order_mx);
  mxDestroyArray(T_mx);
}
