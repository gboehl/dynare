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

ObjectiveMFile::ObjectiveMFile(const std::string& modName, int ntt_arg) :
    ntt {ntt_arg}, ObjectiveMFilename {modName + ".objective.static"}
{
}

void
ObjectiveMFile::unpackSparseMatrixAndCopyIntoTwoDMatData(mxArray* sparseMat, TwoDMatrix& tdm)
{
  int totalCols = mxGetN(sparseMat);
  mwIndex* rowIdxVector = mxGetIr(sparseMat);
  mwIndex* colIdxVector = mxGetJc(sparseMat);

  assert(tdm.ncols() == 3);
  /* Under MATLAB, the following check always holds at equality; under Octave,
     there may be an inequality, because Octave diminishes nzmax if one gives
     zeros in the values vector when calling sparse(). */
  assert(tdm.nrows() >= static_cast<int>(mxGetNzmax(sparseMat)));

  double* ptr = mxGetPr(sparseMat);

  int rind = 0;
  int output_row = 0;

  for (int i = 0; i < totalCols; i++)
    for (int j = 0; j < static_cast<int>((colIdxVector[i + 1] - colIdxVector[i])); j++, rind++)
      {
        tdm.get(output_row, 0) = rowIdxVector[rind] + 1;
        tdm.get(output_row, 1) = i + 1;
        tdm.get(output_row, 2) = ptr[rind];
        output_row++;
      }

  /* If there are less elements than expected (that might happen if some
     derivative is symbolically not zero but numerically zero at the evaluation
     point), then fill in the matrix with empty entries, that will be
     recognized as such by KordpDynare::populateDerivativesContainer() */
  while (output_row < tdm.nrows())
    {
      tdm.get(output_row, 0) = 0;
      tdm.get(output_row, 1) = 0;
      tdm.get(output_row, 2) = 0;
      output_row++;
    }
}

void
ObjectiveMFile::eval(const Vector& y, const Vector& x, const Vector& modParams, Vector& residual,
                     std::vector<TwoDMatrix>& md) noexcept(false)
{
  mxArray* T_m = mxCreateDoubleMatrix(ntt, 1, mxREAL);

  mxArray* y_m = mxCreateDoubleMatrix(y.length(), 1, mxREAL);
  std::copy_n(y.base(), y.length(), mxGetPr(y_m));

  mxArray* x_m = mxCreateDoubleMatrix(1, x.length(), mxREAL);
  std::copy_n(x.base(), x.length(), mxGetPr(x_m));

  mxArray* params_m = mxCreateDoubleMatrix(modParams.length(), 1, mxREAL);
  std::copy_n(modParams.base(), modParams.length(), mxGetPr(params_m));

  mxArray* T_flag_m = mxCreateLogicalScalar(false);

  {
    // Compute temporary terms (for all orders)
    std::string funcname = ObjectiveMFilename + "_g" + std::to_string(md.size()) + "_tt";
    std::array<mxArray*, 1> plhs;
    std::array prhs {T_m, y_m, x_m, params_m};

    int retVal
        = mexCallMATLAB(plhs.size(), plhs.data(), prhs.size(), prhs.data(), funcname.c_str());
    if (retVal != 0)
      throw DynareException(__FILE__, __LINE__, "Trouble calling " + funcname);

    mxDestroyArray(T_m);
    T_m = plhs[0];
  }

  {
    // Compute residuals
    std::string funcname = ObjectiveMFilename + "_resid";
    std::array<mxArray*, 1> plhs;
    std::array prhs {T_m, y_m, x_m, params_m, T_flag_m};

    int retVal
        = mexCallMATLAB(plhs.size(), plhs.data(), prhs.size(), prhs.data(), funcname.c_str());
    if (retVal != 0)
      throw DynareException(__FILE__, __LINE__, "Trouble calling " + funcname);

    residual = Vector {plhs[0]};
    mxDestroyArray(plhs[0]);
  }

  for (size_t i = 1; i <= md.size(); i++)
    {
      // Compute model derivatives
      std::string funcname = ObjectiveMFilename + "_g" + std::to_string(i);
      std::array<mxArray*, 1> plhs;
      std::array prhs {T_m, y_m, x_m, params_m, T_flag_m};

      int retVal
          = mexCallMATLAB(plhs.size(), plhs.data(), prhs.size(), prhs.data(), funcname.c_str());
      if (retVal != 0)
        throw DynareException(__FILE__, __LINE__, "Trouble calling " + funcname);

      if (i == 1)
        {
          assert(static_cast<int>(mxGetM(plhs[0])) == md[i - 1].nrows());
          assert(static_cast<int>(mxGetN(plhs[0])) == md[i - 1].ncols());
          std::copy_n(mxGetPr(plhs[0]), mxGetM(plhs[0]) * mxGetN(plhs[0]), md[i - 1].base());
        }
      else
        unpackSparseMatrixAndCopyIntoTwoDMatData(plhs[0], md[i - 1]);

      mxDestroyArray(plhs[0]);
    }

  mxDestroyArray(T_m);
  mxDestroyArray(y_m);
  mxDestroyArray(x_m);
  mxDestroyArray(params_m);
  mxDestroyArray(T_flag_m);
}
