/*
 * Copyright © 2008-2019 Dynare Team
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

#ifndef _DYNAMIC_DLL_HH
#define _DYNAMIC_DLL_HH

#if defined(_WIN32) || defined(__CYGWIN32__)
# ifndef NOMINMAX
#  define NOMINMAX // Do not define "min" and "max" macros
# endif
# include <windows.h>
#else
# include <dlfcn.h> // unix/linux DLL (.so) handling routines
#endif

#include <string>
#include <memory>

#include "dynamic_abstract_class.hh"

using dynamic_tt_fct = void (*)(const double *y, const double *x, int nb_row_x, const double *params, const double *steady_state, int it_, double *T);
using dynamic_deriv_fct = void (*) (const double *y, const double *x, int nb_row_x, const double *params, const double *steady_state, int it_, const double *T, double *deriv);

/**
 * creates pointer to Dynamic function inside <model>_dynamic.dll
 * and handles calls to it.
 **/
class DynamicModelDLL : public DynamicModelAC
{
private:
  std::vector<dynamic_tt_fct> dynamic_tt;
  std::vector<dynamic_deriv_fct> dynamic_deriv;
#if defined(_WIN32) || defined(__CYGWIN32__)
  HINSTANCE dynamicHinstance;  // DLL instance pointer in Windows
#else
  void *dynamicHinstance; // and in Linux or Mac
#endif
  std::unique_ptr<double[]> tt; // Vector of temporary terms

public:
  // construct and load Dynamic model DLL
  explicit DynamicModelDLL(const std::string &fname, int ntt_arg, int order);
  virtual ~DynamicModelDLL();

  void eval(const Vector &y, const Vector &x, const Vector &params, const Vector &ySteady,
            Vector &residual, std::vector<TwoDMatrix> &md) override;
};
#endif
