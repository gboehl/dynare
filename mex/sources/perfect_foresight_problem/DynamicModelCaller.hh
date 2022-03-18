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

#include <string>
#include <memory>
#include <limits>

#if defined(_WIN32) || defined(__CYGWIN32__)
# ifndef NOMINMAX
#  define NOMINMAX // Do not define "min" and "max" macros
# endif
# include <windows.h>
#else
# include <dlfcn.h> // unix/linux DLL (.so) handling routines
#endif

class DynamicModelCaller
{
public:
  const bool linear;
  bool compute_jacobian; // Not constant, because will be changed from true to false for linear models after first evaluation

  // Used to store error messages (as exceptions cannot cross the OpenMP boundary)
  static std::string error_msg;

  DynamicModelCaller(bool linear_arg, bool compute_jacobian_arg) : linear{linear_arg}, compute_jacobian{compute_jacobian_arg}
  {
  };
  virtual ~DynamicModelCaller() = default;
  virtual double &y(size_t i) const = 0;
  virtual double jacobian(size_t i) const = 0;
  virtual void eval(int it, double *resid) = 0;
};

class DynamicModelDllCaller : public DynamicModelCaller
{
private:
#if defined(_WIN32) || defined(__CYGWIN32__)
  static HINSTANCE dynamic_mex;
#else
  static void *dynamic_mex;
#endif
  using dynamic_tt_fct = void (*)(const double *y, const double *x, int nb_row_x, const double *params, const double *steady_state, int it_, double *T);
  using dynamic_deriv_fct = void (*)(const double *y, const double *x, int nb_row_x, const double *params, const double *steady_state, int it_, const double *T, double *deriv);
  static dynamic_tt_fct residual_tt_fct, g1_tt_fct;
  static dynamic_deriv_fct residual_fct, g1_fct;
  size_t nb_row_x;
  const double *x, *params, *steady_state;
  std::unique_ptr<double[]> tt, y_p, jacobian_p;

public:
  DynamicModelDllCaller(size_t ntt, mwIndex nx, mwIndex ny, size_t ndynvars, const double *x_arg, size_t nb_row_x_arg, const double *params_arg, const double *steady_state_arg, bool linear_arg, bool compute_jacobian_arg);
  virtual ~DynamicModelDllCaller() = default;
  double &
  y(size_t i) const override
  {
    return y_p[i];
  };
  double
  jacobian(size_t i) const override
  {
    return jacobian_p[i];
  };
  void eval(int it, double *resid) override;
  static void load_dll(const std::string &basename);
  static void unload_dll();
};

class DynamicModelMatlabCaller : public DynamicModelCaller
{
private:
  std::string basename;
  mxArray *T_mx, *y_mx, *it_mx, *T_flag_mx, *jacobian_mx, *x_mx, *params_mx, *steady_state_mx;
  /* Given a complex dense matrix (of double floats), returns a real dense matrix of same size.
     Real elements of the original matrix are copied as-is to the new one.
     Complex elements are replaced by NaNs.
     Destroys the original matrix. */
  static mxArray *cmplxToReal(mxArray *m);
public:
  DynamicModelMatlabCaller(std::string basename_arg, size_t ntt, size_t ndynvars, const mxArray *x_mx_arg, const mxArray *params_mx_arg, const mxArray *steady_state_mx_arg, bool linear_arg, bool compute_jacobian_arg);
  ~DynamicModelMatlabCaller() override;
  double &
  y(size_t i) const override
  {
    return mxGetPr(y_mx)[i];
  };
  double
  jacobian(size_t i) const override
  {
    return jacobian_mx ? mxGetPr(jacobian_mx)[i] : std::numeric_limits<double>::quiet_NaN();
  };
  void eval(int it, double *resid) override;
  class Exception
  {
  public:
    const std::string msg;
    Exception(std::string msg_arg) : msg{std::move(msg_arg)}
    {
    };
  };
};
