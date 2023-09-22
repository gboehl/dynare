/*
 * Copyright © 2008-2023 Dynare Team
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
#include <utility>

#include "dynare_exception.hh"

#include "dynamic_abstract_class.hh"

using dynamic_tt_fct = void (*)(const double *y, const double *x, int nb_row_x, const double *params, const double *steady_state, int it_, double *T);
using dynamic_resid_or_g1_fct = void (*)(const double *y, const double *x, int nb_row_x, const double *params, const double *steady_state, int it_, const double *T, double *resid_or_g1);
using dynamic_higher_deriv_fct = void (*)(const double *y, const double *x, int nb_row_x, const double *params, const double *steady_state, int it_, const double *T, double *g_i, double *g_j, double *g_v);

/**
 * creates pointer to Dynamic function inside <model>_dynamic.dll
 * and handles calls to it.
 **/
class DynamicModelDLL : public DynamicModelAC
{
private:
  std::vector<dynamic_tt_fct> dynamic_tt;
  dynamic_resid_or_g1_fct dynamic_resid, dynamic_g1;
  std::vector<dynamic_higher_deriv_fct> dynamic_higher_deriv; // Index 0 is g2
#if defined(_WIN32) || defined(__CYGWIN32__)
  HINSTANCE dynamicHinstance; // DLL instance pointer in Windows
#else
  void *dynamicHinstance; // and in Linux or Mac
#endif
  std::unique_ptr<double[]> tt; // Vector of temporary terms

  template<typename T>
  std::pair<T, dynamic_tt_fct>
  getSymbolsFromDLL(const std::string &funcname, const std::string &fName)
  {
    dynamic_tt_fct tt;
    T deriv;
#if defined(__CYGWIN32__) || defined(_WIN32)
# pragma GCC diagnostic push
# pragma GCC diagnostic ignored "-Wcast-function-type"
    deriv = reinterpret_cast<T>(GetProcAddress(dynamicHinstance, funcname.c_str()));
    tt = reinterpret_cast<dynamic_tt_fct>(GetProcAddress(dynamicHinstance, (funcname + "_tt").c_str()));
# pragma GCC diagnostic pop
#else
      deriv = reinterpret_cast<T>(dlsym(dynamicHinstance, funcname.c_str()));
      tt = reinterpret_cast<dynamic_tt_fct>(dlsym(dynamicHinstance, (funcname + "_tt").c_str()));
#endif
      if (!deriv || !tt)
        {
#if defined(__CYGWIN32__) || defined(_WIN32)
          FreeLibrary(dynamicHinstance);
#else
          dlclose(dynamicHinstance);
#endif
          throw DynareException(__FILE__, __LINE__, "Error when loading symbols from " + fName
#if !defined(__CYGWIN32__) && !defined(_WIN32)
                                + ": " + dlerror()
#endif
                                );
        }
      return { deriv, tt };
  }
public:
  // construct and load Dynamic model DLL
  explicit DynamicModelDLL(const std::string &fname, int ntt_arg, int order);
  virtual ~DynamicModelDLL();

  void eval(const Vector &y, const Vector &x, const Vector &params, const Vector &ySteady,
            Vector &residual, std::vector<TwoDMatrix> &md) override;
};
#endif
