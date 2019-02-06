/*
 * Copyright (C) 2008-2012 Dynare Team
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

#include <sstream>

DynamicModelDLL::DynamicModelDLL(const std::string &modName) noexcept(false)
{
  std::string fName;
#if !defined(__CYGWIN32__) && !defined(_WIN32)
  fName = "./";
#endif
  fName += "+" + modName + "/dynamic" + MEXEXT;

  try
    {
#if defined(__CYGWIN32__) || defined(_WIN32)
      dynamicHinstance = LoadLibrary(fName.c_str());
      if (dynamicHinstance == nullptr)
        throw 1;
      ntt = (int *) GetProcAddress(dynamicHinstance, "ntt");
      dynamic_resid_tt = (dynamic_tt_fct) GetProcAddress(dynamicHinstance, "dynamic_resid_tt");
      dynamic_resid = (dynamic_resid_fct) GetProcAddress(dynamicHinstance, "dynamic_resid");
      dynamic_g1_tt = (dynamic_tt_fct) GetProcAddress(dynamicHinstance, "dynamic_g1_tt");
      dynamic_g1 = (dynamic_g1_fct) GetProcAddress(dynamicHinstance, "dynamic_g1");
      dynamic_g2_tt = (dynamic_tt_fct) GetProcAddress(dynamicHinstance, "dynamic_g2_tt");
      dynamic_g2 = (dynamic_g2_fct) GetProcAddress(dynamicHinstance, "dynamic_g2");
      dynamic_g3_tt = (dynamic_tt_fct) GetProcAddress(dynamicHinstance, "dynamic_g3_tt");
      dynamic_g3 = (dynamic_g3_fct) GetProcAddress(dynamicHinstance, "dynamic_g3");
      if (ntt == nullptr
          || dynamic_resid_tt == nullptr || dynamic_resid == nullptr
          || dynamic_g1_tt == nullptr || dynamic_g1 == nullptr
          || dynamic_g2_tt == nullptr || dynamic_g2 == nullptr
          || dynamic_g3_tt == nullptr || dynamic_g3 == nullptr)
        {
          FreeLibrary(dynamicHinstance); // Free the library
          throw 2;
        }
#else // Linux or Mac
      dynamicHinstance = dlopen(fName.c_str(), RTLD_NOW);
      if (dynamicHinstance == nullptr)
        {
          std::cerr << dlerror() << std::endl;
          throw 1;
        }
      ntt = (int *) dlsym(dynamicHinstance, "ntt");
      dynamic_resid_tt = (dynamic_tt_fct) dlsym(dynamicHinstance, "dynamic_resid_tt");
      dynamic_resid = (dynamic_resid_fct) dlsym(dynamicHinstance, "dynamic_resid");
      dynamic_g1_tt = (dynamic_tt_fct) dlsym(dynamicHinstance, "dynamic_g1_tt");
      dynamic_g1 = (dynamic_g1_fct) dlsym(dynamicHinstance, "dynamic_g1");
      dynamic_g2_tt = (dynamic_tt_fct) dlsym(dynamicHinstance, "dynamic_g2_tt");
      dynamic_g2 = (dynamic_g2_fct) dlsym(dynamicHinstance, "dynamic_g2");
      dynamic_g3_tt = (dynamic_tt_fct) dlsym(dynamicHinstance, "dynamic_g3_tt");
      dynamic_g3 = (dynamic_g3_fct) dlsym(dynamicHinstance, "dynamic_g3");
      if (ntt == nullptr
          || dynamic_resid_tt == nullptr || dynamic_resid == nullptr
          || dynamic_g1_tt == nullptr || dynamic_g1 == nullptr
          || dynamic_g2_tt == nullptr || dynamic_g2 == nullptr
          || dynamic_g3_tt == nullptr || dynamic_g3 == nullptr)
        {
          dlclose(dynamicHinstance); // Free the library
          std::cerr << dlerror() << std::endl;
          throw 2;
        }
#endif

    }
  catch (int i)
    {
      std::ostringstream msg;
      msg << "Error when loading " << fName << " (";
      if (i == 1)
        msg << "can't dynamically load the file";
      if (i == 2)
        msg << "can't locate the relevant dynamic symbols within the MEX file";
      msg << ")";
      throw DynareException(__FILE__, __LINE__, msg.str());
    }
  catch (...)
    {
      throw DynareException(__FILE__, __LINE__, "Can't find the relevant dynamic symbols in " + fName);
    }
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
                      Vector &residual, TwoDMatrix *g1, TwoDMatrix *g2, TwoDMatrix *g3) noexcept(false)
{
  double *T = (double *) malloc(sizeof(double) * (*ntt));
  dynamic_resid_tt(y.base(), x.base(), 1, modParams.base(), ySteady.base(), 0, T);
  dynamic_resid(y.base(), x.base(), 1, modParams.base(), ySteady.base(), 0, T, residual.base());
  if (g1 || g2 || g3)
    dynamic_g1_tt(y.base(), x.base(), 1, modParams.base(), ySteady.base(), 0, T);
  if (g1)
    dynamic_g1(y.base(), x.base(), 1, modParams.base(), ySteady.base(), 0, T, g1->base());
  if (g2 || g3)
    dynamic_g2_tt(y.base(), x.base(), 1, modParams.base(), ySteady.base(), 0, T);
  if (g2)
    dynamic_g2(y.base(), x.base(), 1, modParams.base(), ySteady.base(), 0, T, g2->base());
  if (g3)
    {
      dynamic_g3_tt(y.base(), x.base(), 1, modParams.base(), ySteady.base(), 0, T);
      dynamic_g3(y.base(), x.base(), 1, modParams.base(), ySteady.base(), 0, T, g3->base());
    }
  free(T);
}
