/*
 * Copyright © 2008-2020 Dynare Team
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

#include "dynamic_dll.hh"

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

  dynamic_tt.resize(order+1);

  std::tie(dynamic_resid, dynamic_tt[0]) = getSymbolsFromDLL<dynamic_resid_or_g1_fct>("dynamic_resid", fName);
  std::tie(dynamic_g1, dynamic_tt[1]) = getSymbolsFromDLL<dynamic_resid_or_g1_fct>("dynamic_g1", fName);

  dynamic_higher_deriv.resize(std::max(0, order-1));
  for (int i = 2; i <= order; i++)
    std::tie(dynamic_higher_deriv[i-2], dynamic_tt[i]) = getSymbolsFromDLL<dynamic_higher_deriv_fct>("dynamic_g" + std::to_string(i), fName);

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
  assert(md.size() == dynamic_tt.size()-1);

  for (size_t i = 0; i < dynamic_tt.size(); i++)
    {
      dynamic_tt[i](y.base(), x.base(), 1, modParams.base(), ySteady.base(), 0, tt.get());
      if (i == 0)
        dynamic_resid(y.base(), x.base(), 1, modParams.base(), ySteady.base(), 0, tt.get(), residual.base());
      else if (i == 1)
        dynamic_g1(y.base(), x.base(), 1, modParams.base(), ySteady.base(), 0, tt.get(), md[0].base());
      else
        dynamic_higher_deriv[i-2](y.base(), x.base(), 1, modParams.base(), ySteady.base(), 0, tt.get(), &md[i-1].get(0, 0), &md[i-1].get(0, 1), &md[i-1].get(0, 2));
    }
}
