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

DynamicModelDLL::DynamicModelDLL(const std::string &modName, int ntt_arg)
  : DynamicModelAC(ntt_arg)
{
  std::string fName;
#if !defined(__CYGWIN32__) && !defined(_WIN32)
  fName = "./";
#endif
  fName += "+" + modName + "/dynamic" + MEXEXT;

#if defined(__CYGWIN32__) || defined(_WIN32)
  dynamicHinstance = LoadLibrary(fName.c_str());
  if (!dynamicHinstance)
    throw DynareException(__FILE__, __LINE__, "Error when loading " + fName + ": can't dynamically load the file");
  dynamic_resid_tt = reinterpret_cast<dynamic_tt_fct>(GetProcAddress(dynamicHinstance, "dynamic_resid_tt"));
  dynamic_resid = reinterpret_cast<dynamic_resid_fct>(GetProcAddress(dynamicHinstance, "dynamic_resid"));
  dynamic_g1_tt = reinterpret_cast<dynamic_tt_fct>(GetProcAddress(dynamicHinstance, "dynamic_g1_tt"));
  dynamic_g1 = reinterpret_cast<dynamic_g1_fct>(GetProcAddress(dynamicHinstance, "dynamic_g1"));
  dynamic_g2_tt = reinterpret_cast<dynamic_tt_fct>(GetProcAddress(dynamicHinstance, "dynamic_g2_tt"));
  dynamic_g2 = reinterpret_cast<dynamic_g2_fct>(GetProcAddress(dynamicHinstance, "dynamic_g2"));
  dynamic_g3_tt = reinterpret_cast<dynamic_tt_fct>(GetProcAddress(dynamicHinstance, "dynamic_g3_tt"));
  dynamic_g3 = reinterpret_cast<dynamic_g3_fct>(GetProcAddress(dynamicHinstance, "dynamic_g3"));
  if (!dynamic_resid_tt || !dynamic_resid
      || !dynamic_g1_tt || !dynamic_g1
      || !dynamic_g2_tt || !dynamic_g2
      || !dynamic_g3_tt || !dynamic_g3)
    {
      FreeLibrary(dynamicHinstance); // Free the library
      throw DynareException(__FILE__, __LINE__, "Error when loading " + fName + ": can't locate the relevant dynamic symbols within the MEX file");
    }
#else // Linux or Mac
  dynamicHinstance = dlopen(fName.c_str(), RTLD_NOW);
  if (!dynamicHinstance)
    throw DynareException(__FILE__, __LINE__, "Error when loading " + fName + ": " + dlerror());
  dynamic_resid_tt = reinterpret_cast<dynamic_tt_fct>(dlsym(dynamicHinstance, "dynamic_resid_tt"));
  dynamic_resid = reinterpret_cast<dynamic_resid_fct>(dlsym(dynamicHinstance, "dynamic_resid"));
  dynamic_g1_tt = reinterpret_cast<dynamic_tt_fct>(dlsym(dynamicHinstance, "dynamic_g1_tt"));
  dynamic_g1 = reinterpret_cast<dynamic_g1_fct>(dlsym(dynamicHinstance, "dynamic_g1"));
  dynamic_g2_tt = reinterpret_cast<dynamic_tt_fct>(dlsym(dynamicHinstance, "dynamic_g2_tt"));
  dynamic_g2 = reinterpret_cast<dynamic_g2_fct>(dlsym(dynamicHinstance, "dynamic_g2"));
  dynamic_g3_tt = reinterpret_cast<dynamic_tt_fct>(dlsym(dynamicHinstance, "dynamic_g3_tt"));
  dynamic_g3 = reinterpret_cast<dynamic_g3_fct>(dlsym(dynamicHinstance, "dynamic_g3"));
  if (!dynamic_resid_tt || !dynamic_resid
      || !dynamic_g1_tt || !dynamic_g1
      || !dynamic_g2_tt || !dynamic_g2
      || !dynamic_g3_tt || !dynamic_g3)
    {
      dlclose(dynamicHinstance); // Free the library
      throw DynareException(__FILE__, __LINE__, "Error when loading " + fName + ": " + dlerror());
    }
#endif

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
  dynamic_resid_tt(y.base(), x.base(), 1, modParams.base(), ySteady.base(), 0, tt.get());
  dynamic_resid(y.base(), x.base(), 1, modParams.base(), ySteady.base(), 0, tt.get(), residual.base());
  dynamic_g1_tt(y.base(), x.base(), 1, modParams.base(), ySteady.base(), 0, tt.get());
  dynamic_g1(y.base(), x.base(), 1, modParams.base(), ySteady.base(), 0, tt.get(), md[0].base());
  if (md.size() >= 2)
    {
      dynamic_g2_tt(y.base(), x.base(), 1, modParams.base(), ySteady.base(), 0, tt.get());
      dynamic_g2(y.base(), x.base(), 1, modParams.base(), ySteady.base(), 0, tt.get(), md[1].base());
    }
  if (md.size() >= 3)
    {
      dynamic_g3_tt(y.base(), x.base(), 1, modParams.base(), ySteady.base(), 0, tt.get());
      dynamic_g3(y.base(), x.base(), 1, modParams.base(), ySteady.base(), 0, tt.get(), md[2].base());
    }
}
