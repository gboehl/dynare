/*
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

#include <dynmex.h>

#include "DynamicModelCaller.hh"

#include <algorithm>

std::string DynamicModelCaller::error_msg;

#if defined(_WIN32) || defined(__CYGWIN32__)
HINSTANCE DynamicModelDllCaller::dynamic_mex{nullptr};
#else
void *DynamicModelDllCaller::dynamic_mex{nullptr};
#endif
DynamicModelDllCaller::dynamic_tt_fct DynamicModelDllCaller::residual_tt_fct{nullptr}, DynamicModelDllCaller::g1_tt_fct{nullptr};
DynamicModelDllCaller::dynamic_deriv_fct DynamicModelDllCaller::residual_fct{nullptr}, DynamicModelDllCaller::g1_fct{nullptr};

void
DynamicModelDllCaller::load_dll(const std::string &basename)
{
  // Load symbols from dynamic MEX
  std::string mex_name;
#if !defined(__CYGWIN32__) && !defined(_WIN32)
  mex_name = "./";
#endif
  mex_name += "+" + basename + "/dynamic" + MEXEXT;
#if !defined(__CYGWIN32__) && !defined(_WIN32)
  dynamic_mex = dlopen(mex_name.c_str(), RTLD_NOW);
#else
  dynamic_mex = LoadLibrary(mex_name.c_str());
#endif
  if (!dynamic_mex)
    mexErrMsgTxt("Can't load dynamic MEX file");

#if !defined(__CYGWIN32__) && !defined(_WIN32)
  residual_tt_fct = reinterpret_cast<dynamic_tt_fct>(dlsym(dynamic_mex, "dynamic_resid_tt"));
  residual_fct = reinterpret_cast<dynamic_deriv_fct>(dlsym(dynamic_mex, "dynamic_resid"));
  g1_tt_fct = reinterpret_cast<dynamic_tt_fct>(dlsym(dynamic_mex, "dynamic_g1_tt"));
  g1_fct = reinterpret_cast<dynamic_deriv_fct>(dlsym(dynamic_mex, "dynamic_g1"));
#else
  residual_tt_fct = reinterpret_cast<dynamic_tt_fct>(GetProcAddress(dynamic_mex, "dynamic_resid_tt"));
  residual_fct = reinterpret_cast<dynamic_deriv_fct>(GetProcAddress(dynamic_mex, "dynamic_resid"));
  g1_tt_fct = reinterpret_cast<dynamic_tt_fct>(GetProcAddress(dynamic_mex, "dynamic_g1_tt"));
  g1_fct = reinterpret_cast<dynamic_deriv_fct>(GetProcAddress(dynamic_mex, "dynamic_g1"));
#endif
  if (!residual_tt_fct || !residual_fct || !g1_tt_fct || !g1_fct)
    mexErrMsgTxt("Can't load functions in dynamic MEX file");
}

void
DynamicModelDllCaller::unload_dll()
{
#if !defined(__CYGWIN32__) && !defined(_WIN32)
  dlclose(dynamic_mex);
#else
  FreeLibrary(dynamic_mex);
#endif
}

DynamicModelDllCaller::DynamicModelDllCaller(size_t ntt, mwIndex nx, mwIndex ny, size_t ndynvars, const double *x_arg, size_t nb_row_x_arg, const double *params_arg, const double *steady_state_arg, bool linear_arg, bool compute_jacobian_arg) :
  DynamicModelCaller{linear_arg, compute_jacobian_arg},
  nb_row_x{nb_row_x_arg}, x{x_arg}, params{params_arg}, steady_state{steady_state_arg}
{
  tt = std::make_unique<double[]>(ntt);
  y_p = std::make_unique<double[]>(ndynvars);
  if (compute_jacobian)
    jacobian_p = std::make_unique<double[]>((ndynvars+nx)*ny);
}

void
DynamicModelDllCaller::eval(int it, double *resid)
{
  residual_tt_fct(y_p.get(), x, nb_row_x, params, steady_state, it, tt.get());
  residual_fct(y_p.get(), x, nb_row_x, params, steady_state, it, tt.get(), resid);
  if (compute_jacobian)
    {
      g1_tt_fct(y_p.get(), x, nb_row_x, params, steady_state, it, tt.get());
      g1_fct(y_p.get(), x, nb_row_x, params, steady_state, it, tt.get(), jacobian_p.get());

      if (linear)
        compute_jacobian = false; // If model is linear, no need to recompute Jacobian later
    }
}

DynamicModelMatlabCaller::DynamicModelMatlabCaller(std::string basename_arg, size_t ntt, size_t ndynvars, const mxArray *x_mx_arg, const mxArray *params_mx_arg, const mxArray *steady_state_mx_arg, bool linear_arg, bool compute_jacobian_arg) :
  DynamicModelCaller{linear_arg, compute_jacobian_arg},
  basename{std::move(basename_arg)},
  T_mx{mxCreateDoubleMatrix(ntt, 1, mxREAL)},
  y_mx{mxCreateDoubleMatrix(ndynvars, 1, mxREAL)},
  it_mx{mxCreateDoubleScalar(0)},
  T_flag_mx{mxCreateLogicalScalar(false)},
  jacobian_mx{nullptr},
  x_mx{mxDuplicateArray(x_mx_arg)},
  params_mx{mxDuplicateArray(params_mx_arg)},
  steady_state_mx{mxDuplicateArray(steady_state_mx_arg)}
{
}

DynamicModelMatlabCaller::~DynamicModelMatlabCaller()
{
  mxDestroyArray(T_mx);
  mxDestroyArray(y_mx);
  mxDestroyArray(it_mx);
  mxDestroyArray(T_flag_mx);
  if (jacobian_mx)
    mxDestroyArray(jacobian_mx);
  mxDestroyArray(x_mx);
  mxDestroyArray(params_mx);
  mxDestroyArray(steady_state_mx);
}

void
DynamicModelMatlabCaller::eval(int it, double *resid)
{
  *mxGetPr(it_mx) = it + 1;

  {
    // Compute temporary terms
    std::string funcname = basename + (compute_jacobian ? ".dynamic_g1_tt" : ".dynamic_resid_tt");
    mxArray *plhs[1], *prhs[] = { T_mx, y_mx, x_mx, params_mx, steady_state_mx, it_mx };

    mxArray *exception = mexCallMATLABWithTrap(1, plhs, 6, prhs, funcname.c_str());
    if (exception)
      {
        error_msg = std::string{"An error occurred when calling "} + funcname;
        return; // Avoid manipulating null pointers in plhs, see #1832
      }

    mxDestroyArray(T_mx);
    T_mx = plhs[0];
  }

  {
    // Compute residuals
    std::string funcname = basename + ".dynamic_resid";
    mxArray *plhs[1], *prhs[] = { T_mx, y_mx, x_mx, params_mx, steady_state_mx, it_mx, T_flag_mx };

    mxArray *exception = mexCallMATLABWithTrap(1, plhs, 7, prhs, funcname.c_str());
    if (exception)
      {
        error_msg = std::string{"An error occurred when calling "} + funcname;
        return; // Avoid manipulating null pointers in plhs, see #1832
      }

    if (!mxIsDouble(plhs[0]) || mxIsSparse(plhs[0]))
      {
        error_msg = "Residuals should be a dense array of double floats";
        return;
      }

    if (mxIsComplex(plhs[0]))
      plhs[0] = cmplxToReal(plhs[0]);

    std::copy_n(mxGetPr(plhs[0]), mxGetNumberOfElements(plhs[0]), resid);
    mxDestroyArray(plhs[0]);
  }

  if (compute_jacobian)
    {
      // Compute Jacobian
      std::string funcname = basename + ".dynamic_g1";
      mxArray *plhs[1], *prhs[] = { T_mx, y_mx, x_mx, params_mx, steady_state_mx, it_mx, T_flag_mx };

      mxArray *exception = mexCallMATLABWithTrap(1, plhs, 7, prhs, funcname.c_str());
      if (exception)
        {
          error_msg = std::string{"An error occurred when calling "} + funcname;
          return; // Avoid manipulating null pointers in plhs, see #1832
        }

      if (jacobian_mx)
        {
          mxDestroyArray(jacobian_mx);
          jacobian_mx = nullptr;
        }

      if (!mxIsDouble(plhs[0]) || mxIsSparse(plhs[0]))
        {
          error_msg = "Jacobian should be a dense array of double floats";
          return;
        }

      if (mxIsComplex(plhs[0]))
        jacobian_mx = cmplxToReal(plhs[0]);
      else
        jacobian_mx = plhs[0];

      if (linear)
        compute_jacobian = false; // If model is linear, no need to recompute Jacobian later
    }
}

/* NB: This is a duplicate of DynamicModelMFile::cmplxToReal() in
   k_order_perturbation MEX */
mxArray *
DynamicModelMatlabCaller::cmplxToReal(mxArray *cmplx_mx)
{
  mxArray *real_mx = mxCreateDoubleMatrix(mxGetM(cmplx_mx), mxGetN(cmplx_mx), mxREAL);

#if MX_HAS_INTERLEAVED_COMPLEX
  mxComplexDouble *cmplx = mxGetComplexDoubles(cmplx_mx);
#else
  double *cmplx_real = mxGetPr(cmplx_mx);
  double *cmplx_imag = mxGetPi(cmplx_mx);
#endif
  double *real = mxGetPr(real_mx);

  for (size_t i = 0; i < mxGetNumberOfElements(cmplx_mx); i++)
#if MX_HAS_INTERLEAVED_COMPLEX
    if (cmplx[i].imag == 0.0)
      real[i] = cmplx[i].real;
#else
    if (cmplx_imag[i] == 0.0)
      real[i] = cmplx_real[i];
#endif
    else
      real[i] = std::numeric_limits<double>::quiet_NaN();

  mxDestroyArray(cmplx_mx);
  return real_mx;
}
