/*
 * Copyright (C) 2010-2019 Dynare Team
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

#include <algorithm>
#include <cassert>

#include "dynare_exception.hh"

#include "dynamic_m.hh"

DynamicModelMFile::DynamicModelMFile(const std::string &modName, int ntt_arg) :
  DynamicModelAC(ntt_arg),
  DynamicMFilename{modName + ".dynamic"}
{
}

void
DynamicModelMFile::unpackSparseMatrixAndCopyIntoTwoDMatData(mxArray *sparseMat, TwoDMatrix &tdm)
{
  int totalCols = mxGetN(sparseMat);
  mwIndex *rowIdxVector = mxGetIr(sparseMat);
  mwSize sizeRowIdxVector = mxGetNzmax(sparseMat);
  mwIndex *colIdxVector = mxGetJc(sparseMat);

  assert(tdm.ncols() == 3);
  assert(tdm.nrows() == sizeRowIdxVector);

  double *ptr = mxGetPr(sparseMat);

  int rind = 0;
  int output_row = 0;

  for (int i = 0; i < totalCols; i++)
    for (int j = 0; j < static_cast<int>((colIdxVector[i+1]-colIdxVector[i])); j++, rind++)
      {
        tdm.get(output_row, 0) = rowIdxVector[rind] + 1;
        tdm.get(output_row, 1) = i + 1;
        tdm.get(output_row, 2) = ptr[rind];
        output_row++;
      }

  /* If there are less elements than Nzmax (that might happen if some
     derivative is symbolically not zero but numerically zero at the evaluation
     point), then fill in the matrix with empty entries, that will be
     recognized as such by KordpDynare::populateDerivativesContainer() */
  while (output_row < static_cast<int>(sizeRowIdxVector))
    {
      tdm.get(output_row, 0) = 0;
      tdm.get(output_row, 1) = 0;
      tdm.get(output_row, 2) = 0;
      output_row++;
    }
}

void
DynamicModelMFile::eval(const Vector &y, const Vector &x, const Vector &modParams, const Vector &ySteady,
                        Vector &residual, std::vector<TwoDMatrix> &md) noexcept(false)
{
  constexpr int nlhs_dynamic = 4, nrhs_dynamic = 5;
  mxArray *prhs[nrhs_dynamic], *plhs[nlhs_dynamic];

  prhs[0] = mxCreateDoubleMatrix(y.length(), 1, mxREAL);
  prhs[1] = mxCreateDoubleMatrix(1, x.length(), mxREAL);
  prhs[2] = mxCreateDoubleMatrix(modParams.length(), 1, mxREAL);
  prhs[3] = mxCreateDoubleMatrix(ySteady.length(), 1, mxREAL);
  prhs[4] = mxCreateDoubleScalar(1.0);

  std::copy_n(y.base(), y.length(), mxGetPr(prhs[0]));
  std::copy_n(x.base(), x.length(), mxGetPr(prhs[1]));
  std::copy_n(modParams.base(), modParams.length(), mxGetPr(prhs[2]));
  std::copy_n(ySteady.base(), ySteady.length(), mxGetPr(prhs[3]));

  int retVal = mexCallMATLAB(nlhs_dynamic, plhs, nrhs_dynamic, prhs, DynamicMFilename.c_str());
  if (retVal != 0)
    throw DynareException(__FILE__, __LINE__, "Trouble calling " + DynamicMFilename);

  residual = Vector{plhs[0]};

  assert(static_cast<int>(mxGetM(plhs[1])) == md[0].nrows());
  assert(static_cast<int>(mxGetN(plhs[1])) == md[0].ncols());
  std::copy_n(mxGetPr(plhs[1]), mxGetM(plhs[1])*mxGetN(plhs[1]), md[0].base());

  if (md.size() >= 2)
    unpackSparseMatrixAndCopyIntoTwoDMatData(plhs[2], md[1]);
  if (md.size() >= 3)
    unpackSparseMatrixAndCopyIntoTwoDMatData(plhs[3], md[2]);

  for (int i = 0; i < nrhs_dynamic; i++)
    mxDestroyArray(prhs[i]);
  for (int i = 0; i < nlhs_dynamic; i++)
    mxDestroyArray(plhs[i]);
}
