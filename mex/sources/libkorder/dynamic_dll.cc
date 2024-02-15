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

#include "dynamic_dll.hh"

#include <cassert>
#include <iostream>

DynamicModelDLL::DynamicModelDLL(const std::string& modName, int order_arg,
                                 const mxArray* dynamic_g1_sparse_rowval_mx_arg,
                                 const mxArray* dynamic_g1_sparse_colval_mx_arg,
                                 const mxArray* dynamic_g1_sparse_colptr_mx_arg,
                                 const std::vector<const mxArray*> dynamic_gN_sparse_indices_arg,
                                 int ntt) :
    DynamicModelAC {order_arg, dynamic_g1_sparse_rowval_mx_arg, dynamic_g1_sparse_colval_mx_arg,
                    dynamic_g1_sparse_colptr_mx_arg, dynamic_gN_sparse_indices_arg},
    dynamic_fcts(order + 1),
    dynamic_fcts_tt(order + 1),
    tt_tmp(ntt),
    derivatives_tmp(order)
{
  using namespace std::string_literals;

  for (int o {0}; o <= order; o++)
    {
      std::string func_name {"dynamic_"s + (o == 0 ? "resid" : "g"s + std::to_string(o))};
      std::string mex_filename
      {
#if !defined(__CYGWIN32__) && !defined(_WIN32)
        "./"s +
#endif
            "+"s + modName + "/+sparse/" + func_name + MEXEXT
      };

      mex_handle_t handle {load_mex(mex_filename)};
      if (handle)
        mex_handles.push_back(handle);
      else
        {
          std::ranges::for_each(mex_handles, unload_mex);

          throw DynareException(__FILE__, __LINE__,
                                "Error when loading " + mex_filename
#if !defined(__CYGWIN32__) && !defined(_WIN32)
                                    + ": " + dlerror()
#endif
          );
        }

      std::tie(dynamic_fcts[o], dynamic_fcts_tt[o])
          = getSymbolsFromDLL(func_name, mex_handles.back());
      if (!dynamic_fcts[o] || !dynamic_fcts_tt[o])
        {
          std::ranges::for_each(mex_handles, unload_mex);

          throw DynareException(__FILE__, __LINE__,
                                "Error when loading symbols from " + mex_filename
#if !defined(__CYGWIN32__) && !defined(_WIN32)
                                    + ": " + dlerror()
#endif
          );
        }
    }

  derivatives_tmp[0].resize(mxGetNumberOfElements(dynamic_g1_sparse_rowval_mx));
  for (int o {2}; o <= order; o++)
    derivatives_tmp[o - 1].resize(mxGetM(dynamic_gN_sparse_indices[o - 2]));
}

DynamicModelDLL::~DynamicModelDLL()
{
  std::ranges::for_each(mex_handles, unload_mex);
}

DynamicModelDLL::mex_handle_t
DynamicModelDLL::load_mex(const std::string& mex_filename)
{
#if defined(__CYGWIN32__) || defined(_WIN32)
  return LoadLibrary(mex_filename.c_str());
#else // GNU/Linux or Mac
  return dlopen(mex_filename.c_str(), RTLD_NOW);
#endif
}

void
DynamicModelDLL::unload_mex(mex_handle_t handle)
{
#if defined(__CYGWIN32__) || defined(_WIN32)
  FreeLibrary(handle);
#else
  dlclose(handle);
#endif
}

std::pair<dynamic_fct, dynamic_tt_fct>
DynamicModelDLL::getSymbolsFromDLL(const std::string& func_name, mex_handle_t handle)
{
#if defined(__CYGWIN32__) || defined(_WIN32)
# pragma GCC diagnostic push
# pragma GCC diagnostic ignored "-Wcast-function-type"
  return {reinterpret_cast<dynamic_fct>(GetProcAddress(handle, func_name.c_str())),
          reinterpret_cast<dynamic_tt_fct>(GetProcAddress(handle, (func_name + "_tt").c_str()))};
# pragma GCC diagnostic pop
#else
  return {reinterpret_cast<dynamic_fct>(dlsym(handle, func_name.c_str())),
          reinterpret_cast<dynamic_tt_fct>(dlsym(handle, (func_name + "_tt").c_str()))};
#endif
}

void
DynamicModelDLL::eval(const Vector& y, const Vector& x, const Vector& modParams,
                      const Vector& ySteady, Vector& residual, const std::map<int, int>& dynToDynpp,
                      TensorContainer<FSSparseTensor>& derivatives) noexcept(false)
{
  for (int o {0}; o <= order; o++)
    {
      dynamic_fcts_tt[o](y.base(), x.base(), modParams.base(), ySteady.base(), tt_tmp.data());
      dynamic_fcts[o](y.base(), x.base(), modParams.base(), ySteady.base(), tt_tmp.data(),
                      o == 0 ? residual.base() : derivatives_tmp[o - 1].data());

      if (o >= 1)
        {
          IntSequence s(o, 0);
          auto tensor = std::make_unique<FSSparseTensor>(o, dynToDynpp.size(), residual.length());
          size_t nnz {derivatives_tmp[o - 1].size()};

          for (size_t k {0}; k < nnz; k++)
            {
              int rowval;
              if (o == 1)
                {
#if MX_HAS_INTERLEAVED_COMPLEX
                  const int32_T* sparse_rowval {mxGetInt32s(dynamic_g1_sparse_rowval_mx)};
                  const int32_T* sparse_colval {mxGetInt32s(dynamic_g1_sparse_colval_mx)};
#else
                  const int32_T* sparse_rowval {
                      static_cast<const int32_T*>(mxGetData(dynamic_g1_sparse_rowval_mx))};
                  const int32_T* sparse_colval {
                      static_cast<const int32_T*>(mxGetData(dynamic_g1_sparse_colval_mx))};
#endif
                  s[0] = dynToDynpp.at(sparse_colval[k] - 1);
                  rowval = sparse_rowval[k] - 1;
                }
              else
                {
                  const mxArray* sparse_indices_mx {dynamic_gN_sparse_indices[o - 2]};
#if MX_HAS_INTERLEAVED_COMPLEX
                  const int32_T* sparse_indices {mxGetInt32s(sparse_indices_mx)};
#else
                  const int32_T* sparse_indices {
                      static_cast<const int32_T*>(mxGetData(sparse_indices_mx))};
#endif
                  for (int i {0}; i < o; i++)
                    s[i] = dynToDynpp.at(sparse_indices[k + (i + 1) * nnz] - 1);
                  s.sort();
                  rowval = sparse_indices[k] - 1;
                }
              tensor->insert(s, rowval, derivatives_tmp[o - 1][k]);
            }

          derivatives.insert(std::move(tensor));
        }
    }
}
