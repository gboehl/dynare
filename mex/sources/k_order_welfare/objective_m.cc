/*
 * Copyright © 2021-2024 Dynare Team
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
#include <array>
#include <cassert>
#include <type_traits>

#include "dynare_exception.hh"

#include "objective_m.hh"

ObjectiveMFile::ObjectiveMFile(const std::string& modName, int kOrder_arg,
                               const mxArray* objective_g1_sparse_rowval_mx_arg,
                               const mxArray* objective_g1_sparse_colval_mx_arg,
                               const mxArray* objective_g1_sparse_colptr_mx_arg,
                               const std::vector<const mxArray*> objective_gN_sparse_indices_arg) :
    ObjectiveMFilename {modName + ".objective.sparse.static"},
    kOrder {kOrder_arg},
    objective_g1_sparse_rowval_mx {objective_g1_sparse_rowval_mx_arg},
    objective_g1_sparse_colval_mx {objective_g1_sparse_colval_mx_arg},
    objective_g1_sparse_colptr_mx {objective_g1_sparse_colptr_mx_arg},
    objective_gN_sparse_indices {objective_gN_sparse_indices_arg}
{
}

void
ObjectiveMFile::eval(const Vector& y, const Vector& x, const Vector& modParams, Vector& residual,
                     const std::vector<int>& dynToDynpp,
                     TensorContainer<FSSparseTensor>& derivatives) const
{
  mxArray* y_mx = mxCreateDoubleMatrix(y.length(), 1, mxREAL);
  std::copy_n(y.base(), y.length(), mxGetPr(y_mx));

  mxArray* x_mx = mxCreateDoubleMatrix(1, x.length(), mxREAL);
  std::copy_n(x.base(), x.length(), mxGetPr(x_mx));

  mxArray* params_mx = mxCreateDoubleMatrix(modParams.length(), 1, mxREAL);
  std::copy_n(modParams.base(), modParams.length(), mxGetPr(params_mx));

  mxArray *T_order_mx, *T_mx;

  {
    // Compute residuals
    std::string funcname = ObjectiveMFilename + "_resid";
    std::array<mxArray*, 3> plhs;
    std::array prhs {y_mx, x_mx, params_mx};

    int retVal
        = mexCallMATLAB(plhs.size(), plhs.data(), prhs.size(), prhs.data(), funcname.c_str());
    if (retVal != 0)
      throw DynareException(__FILE__, __LINE__, "Trouble calling " + funcname);

    residual = Vector {plhs[0]};
    mxDestroyArray(plhs[0]);

    T_order_mx = plhs[1];
    T_mx = plhs[2];
  }

  {
    // Compute Jacobian
    std::string funcname = ObjectiveMFilename + "_g1";
    std::array<mxArray*, 3> plhs;
    std::array prhs {y_mx,
                     x_mx,
                     params_mx,
                     const_cast<mxArray*>(objective_g1_sparse_rowval_mx),
                     const_cast<mxArray*>(objective_g1_sparse_colval_mx),
                     const_cast<mxArray*>(objective_g1_sparse_colptr_mx),
                     T_order_mx,
                     T_mx};

    int retVal
        = mexCallMATLAB(plhs.size(), plhs.data(), prhs.size(), prhs.data(), funcname.c_str());
    if (retVal != 0)
      throw DynareException(__FILE__, __LINE__, "Trouble calling " + funcname);

    assert(static_cast<int>(mxGetN(plhs[0])) == y.length());
    double* g1_v {mxGetPr(plhs[0])};
    mwIndex* g1_ir {mxGetIr(plhs[0])};
    mwIndex* g1_jc {mxGetJc(plhs[0])};

    IntSequence s(1, 0);
    auto tensor = std::make_unique<FSSparseTensor>(1, y.length(), 1);

    for (int j {0}; j < y.length(); j++)
      for (mwIndex k {g1_jc[j]}; k < g1_jc[j + 1]; k++)
        {
          s[0] = dynToDynpp[j];
          tensor->insert(s, g1_ir[k], g1_v[k]);
        }

    mxDestroyArray(plhs[0]);

    mxDestroyArray(T_order_mx);
    T_order_mx = plhs[1];

    mxDestroyArray(T_mx);
    T_mx = plhs[2];

    derivatives.insert(std::move(tensor));
  }

  for (int o {2}; o <= kOrder; o++)
    {
      // Compute higher derivatives
      std::string funcname = ObjectiveMFilename + "_g" + std::to_string(o);
      std::array<mxArray*, 3> plhs;
      std::array prhs {y_mx, x_mx, params_mx, T_order_mx, T_mx};

      int retVal
          = mexCallMATLAB(plhs.size(), plhs.data(), prhs.size(), prhs.data(), funcname.c_str());
      if (retVal != 0)
        throw DynareException(__FILE__, __LINE__, "Trouble calling " + funcname);

      const mxArray* sparse_indices_mx {objective_gN_sparse_indices[o - 2]};
      size_t nnz {mxGetM(sparse_indices_mx)};
#if MX_HAS_INTERLEAVED_COMPLEX
      const int32_T* sparse_indices {mxGetInt32s(sparse_indices_mx)};
#else
      const int32_T* sparse_indices {static_cast<const int32_T*>(mxGetData(sparse_indices_mx))};
#endif

      assert(mxGetNumberOfElements(plhs[0]) == nnz);
      double* gN_v {mxGetPr(plhs[0])};

      IntSequence s(o, 0);
      auto tensor = std::make_unique<FSSparseTensor>(o, y.length(), 1);

      for (size_t k {0}; k < nnz; k++)
        {
          for (int i {0}; i < o; i++)
            s[i] = dynToDynpp[sparse_indices[k + (i + 1) * nnz] - 1];
          assert(s.isSorted());
          tensor->insert(s, sparse_indices[k] - 1, gN_v[k]);
        }

      mxDestroyArray(plhs[0]);

      mxDestroyArray(T_order_mx);
      T_order_mx = plhs[1];

      mxDestroyArray(T_mx);
      T_mx = plhs[2];

      derivatives.insert(std::move(tensor));
    }

  mxDestroyArray(y_mx);
  mxDestroyArray(x_mx);
  mxDestroyArray(params_mx);
  mxDestroyArray(T_order_mx);
  mxDestroyArray(T_mx);
}
