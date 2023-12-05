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

#include <algorithm>
#include <limits>
#include <memory>
#include <string>
#include <type_traits>

#include <dynmex.h>

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
  bool compute_jacobian; // Not constant, because will be changed from true to false for linear
                         // models after first evaluation

  // Used to store error messages (as exceptions cannot cross the OpenMP boundary)
  static std::string error_msg;

  DynamicModelCaller(bool linear_arg, bool compute_jacobian_arg) :
      linear {linear_arg}, compute_jacobian {compute_jacobian_arg}
  {
  }
  virtual ~DynamicModelCaller() = default;
  [[nodiscard]] virtual double* y() const = 0;
  [[nodiscard]] virtual double* x() const = 0;
  /* Copy a column of the Jacobian to dest.
     Only copies non-zero elements, according to g1_sparse_{rowval,colval,colptr}. */
  virtual void copy_jacobian_column(mwIndex col, double* dest) const = 0;
  virtual void eval(double* resid) = 0;
};

class DynamicModelDllCaller : public DynamicModelCaller
{
private:
#if defined(_WIN32) || defined(__CYGWIN32__)
  static HINSTANCE resid_mex, g1_mex;
#else
  static void *resid_mex, *g1_mex;
#endif
  using dynamic_tt_fct = void (*)(const double* y, const double* x, const double* params,
                                  const double* steady_state, double* T);
  using dynamic_fct = void (*)(const double* y, const double* x, const double* params,
                               const double* steady_state, const double* T, double* value);
  static dynamic_tt_fct residual_tt_fct, g1_tt_fct;
  static dynamic_fct residual_fct, g1_fct;
  const double *params, *steady_state;
  std::unique_ptr<double[]> tt, y_p, x_p, jacobian_p;
  const int32_T* g1_sparse_colptr;

public:
  DynamicModelDllCaller(size_t ntt, mwIndex ny, mwIndex nx, const double* params_arg,
                        const double* steady_state_arg, const int32_T* g1_sparse_colptr_arg,
                        bool linear_arg, bool compute_jacobian_arg);
  virtual ~DynamicModelDllCaller() = default;
  [[nodiscard]] double*
  y() const override
  {
    return y_p.get();
  }
  [[nodiscard]] double*
  x() const override
  {
    return x_p.get();
  }
  void copy_jacobian_column(mwIndex col, double* dest) const override;
  void eval(double* resid) override;
  static void load_dll(const std::string& basename);
  static void unload_dll();
};

class DynamicModelMatlabCaller : public DynamicModelCaller
{
private:
  std::string basename;
  mxArray *y_mx, *x_mx, *jacobian_mx, *params_mx, *steady_state_mx, *g1_sparse_rowval_mx,
      *g1_sparse_colval_mx, *g1_sparse_colptr_mx;
  /* Given a complex matrix (of double floats), returns a sparse matrix of same size.
     Real elements of the original matrix are copied as-is to the new one.
     Complex elements are replaced by NaNs.
     Destroys the original matrix.
     There are two versions, one for dense matrices, another for sparse
     matrices. */
  template<bool sparse>
  static mxArray* cmplxToReal(mxArray* cmplx_mx);

public:
  DynamicModelMatlabCaller(std::string basename_arg, mwIndex ny, mwIndex nx,
                           const mxArray* params_mx_arg, const mxArray* steady_state_mx_arg,
                           const mxArray* g1_sparse_rowval_mx_arg,
                           const mxArray* g1_sparse_colval_mx_arg,
                           const mxArray* g1_sparse_colptr_mx_arg, bool linear_arg,
                           bool compute_jacobian_arg);
  ~DynamicModelMatlabCaller() override;
  [[nodiscard]] double*
  y() const override
  {
    return mxGetPr(y_mx);
  }
  [[nodiscard]] double*
  x() const override
  {
    return mxGetPr(x_mx);
  }
  void copy_jacobian_column(mwIndex col, double* dest) const override;
  void eval(double* resid) override;
  class Exception
  {
  public:
    const std::string msg;
    Exception(std::string msg_arg) : msg {std::move(msg_arg)}
    {
    }
  };
};

template<bool sparse>
mxArray*
DynamicModelMatlabCaller::cmplxToReal(mxArray* cmplx_mx)
{
  mxArray* real_mx {
      sparse ? mxCreateSparse(mxGetM(cmplx_mx), mxGetN(cmplx_mx), mxGetNzmax(cmplx_mx), mxREAL)
             : mxCreateDoubleMatrix(mxGetM(cmplx_mx), mxGetN(cmplx_mx), mxREAL)};

  if constexpr (sparse)
    {
      std::copy_n(mxGetIr(cmplx_mx), mxGetNzmax(cmplx_mx), mxGetIr(real_mx));
      std::copy_n(mxGetJc(cmplx_mx), mxGetN(cmplx_mx) + 1, mxGetJc(real_mx));
    }

#if MX_HAS_INTERLEAVED_COMPLEX
  mxComplexDouble* cmplx {mxGetComplexDoubles(cmplx_mx)};
#else
  double* cmplx_real {mxGetPr(cmplx_mx)};
  double* cmplx_imag {mxGetPi(cmplx_mx)};
#endif
  double* real {mxGetPr(real_mx)};
  for (std::conditional_t<sparse, mwSize, size_t> i {0};
       i <
       [&] {
         if constexpr (sparse)
           return mxGetNzmax(cmplx_mx);
         else
           return mxGetNumberOfElements(cmplx_mx);
       }(); // Use a lambda instead of the ternary operator to have the right type (there is no
            // constexpr ternary operator)
       i++)
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
