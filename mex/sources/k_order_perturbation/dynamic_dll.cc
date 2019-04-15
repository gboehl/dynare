/*
 * Copyright (C) 2008-2019 Dynare Team
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

#include "dynamic_dll.hh"

#include "dynare_exception.hh"

#include <iostream>
#include <cassert>

DynamicModelDLL::DynamicModelDLL(const std::string &modName, int ntt_arg, int order)
  : DynamicModelAC(ntt_arg)
{
  std::string fName;
#if !defined(__CYGWIN32__) && !defined(_WIN32)
  fName = "./";
#endif
  fName += "+" + modName + "/dynamic" + MEXEXT;

#if defined(__CYGWIN32__) || defined(_WIN32)
  dynamicHinstance = LoadLibrary(fName.c_str());
#else // GNU/Linux or Mac
  dynamicHinstance = dlopen(fName.c_str(), RTLD_NOW);
#endif
  if (!dynamicHinstance)
    throw DynareException(__FILE__, __LINE__, "Error when loading " + fName
#if !defined(__CYGWIN32__) && !defined(_WIN32)
                          + ": " + dlerror()
#endif
                          );

  for (int i = 0; i <= order; i++)
    {
      std::string funcname = "dynamic_" + (i == 0 ? "resid" : "g" + std::to_string(i));
      void *deriv, *tt;
#if defined(__CYGWIN32__) || defined(_WIN32)
      deriv = GetProcAddress(dynamicHinstance, funcname.c_str());
      tt = GetProcAddress(dynamicHinstance, (funcname + "_tt").c_str());
#else
      deriv = dlsym(dynamicHinstance, funcname.c_str());
      tt = dlsym(dynamicHinstance, (funcname + "_tt").c_str());
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
      dynamic_deriv.push_back(reinterpret_cast<dynamic_deriv_fct>(deriv));
      dynamic_tt.push_back(reinterpret_cast<dynamic_tt_fct>(tt));
    }

  tt = std::make_unique<double[]>(ntt);
}

DynamicModelDLL::~DynamicModelDLL()
{
#if defined(__CYGWIN32__) || defined(_WIN32)
  auto result = FreeLibrary(dynamicHinstance);
  if (result == 0)
    {
      std::cerr << "Can't free the *_dynamic DLL" << std::endl;
      exit(EXIT_FAILURE);
    }
#else
  dlclose(dynamicHinstance);
#endif
}

void
DynamicModelDLL::eval(const Vector &y, const Vector &x, const Vector &modParams, const Vector &ySteady,
                      Vector &residual, std::vector<TwoDMatrix> &md) noexcept(false)
{
  assert(md.size() == dynamic_deriv.size()-1);

  for (size_t i = 0; i < dynamic_deriv.size(); i++)
    {
      dynamic_tt[i](y.base(), x.base(), 1, modParams.base(), ySteady.base(), 0, tt.get());
      dynamic_deriv[i](y.base(), x.base(), 1, modParams.base(), ySteady.base(), 0, tt.get(), i == 0 ? residual.base() : md[i-1].base());
    }
}
