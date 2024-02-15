/*
 * Copyright © 2008-2024 Dynare Team
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

#ifndef DYNAMIC_DLL_HH
#define DYNAMIC_DLL_HH

#if defined(_WIN32) || defined(__CYGWIN32__)
# ifndef NOMINMAX
#  define NOMINMAX // Do not define "min" and "max" macros
# endif
# include <windows.h>
#else
# include <dlfcn.h> // unix/linux DLL (.so) handling routines
#endif

#include <string>
#include <utility>
#include <vector>

#include "dynare_exception.hh"

#include "dynamic_abstract_class.hh"

using dynamic_tt_fct = void (*)(const double* y, const double* x, const double* params,
                                const double* steady_state, double* T);
using dynamic_fct = void (*)(const double* y, const double* x, const double* params,
                             const double* steady_state, const double* T, double* values);

/**
 * creates pointer to Dynamic function inside <model>_dynamic.dll
 * and handles calls to it.
 **/
class DynamicModelDLL : public DynamicModelAC
{
private:
  using mex_handle_t =
#if defined(_WIN32) || defined(__CYGWIN32__)
      HINSTANCE
#else
      void*
#endif
      ;
  std::vector<mex_handle_t> mex_handles;
  std::vector<dynamic_fct> dynamic_fcts; // Index 0 is resid
  std::vector<dynamic_tt_fct> dynamic_fcts_tt;
  std::vector<double> tt_tmp; // Vector of temporary terms
  std::vector<std::vector<double>> derivatives_tmp;

public:
  // construct and load Dynamic model DLL
  DynamicModelDLL(const std::string& fname, int order_arg,
                  const mxArray* dynamic_g1_sparse_rowval_mx_arg,
                  const mxArray* dynamic_g1_sparse_colval_mx_arg,
                  const mxArray* dynamic_g1_sparse_colptr_mx_arg,
                  const std::vector<const mxArray*> dynamic_gN_sparse_indices_arg, int ntt);
  ~DynamicModelDLL() override;

  void eval(const Vector& y, const Vector& x, const Vector& params, const Vector& ySteady,
            Vector& residual, const std::map<int, int>& dynToDynpp,
            TensorContainer<FSSparseTensor>& derivatives) noexcept(false) override;

private:
  // Returns a non-null pointer on success
  static mex_handle_t load_mex(const std::string& mex_filename);
  static void unload_mex(mex_handle_t handle);
  // Returns non-null pointers on success
  static std::pair<dynamic_fct, dynamic_tt_fct> getSymbolsFromDLL(const std::string& func_name,
                                                                  mex_handle_t handle);
};
#endif
