/*
 * Copyright (C) 2010-2014 Dynare Team
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

#include <assert.h>

#include "dynamic_abstract_class.hh"

DynamicModelAC::~DynamicModelAC()
{
}

void
DynamicModelAC::copyDoubleIntoTwoDMatData(double *dm, TwoDMatrix *tdm, int rows, int cols)
{
  assert(rows == tdm->nrows());
  assert(cols == tdm->ncols());

  int dmIdx = 0;
  for (int j = 0; j < cols; j++)
    for (int i = 0; i < rows; i++)
      tdm->get(i, j) = dm[dmIdx++];
}

void
DynamicModelAC::unpackSparseMatrixAndCopyIntoTwoDMatData(mxArray *sparseMat, TwoDMatrix *tdm)
{
  int totalCols = mxGetN(sparseMat);
  mwIndex *rowIdxVector = mxGetIr(sparseMat);
  mwSize sizeRowIdxVector = mxGetNzmax(sparseMat);
  mwIndex *colIdxVector = mxGetJc(sparseMat);

  assert(tdm->ncols() == 3);
  assert(tdm->nrows() == sizeRowIdxVector);

  double *ptr = mxGetPr(sparseMat);

  int rind = 0;
  int output_row = 0;

  for (int i = 0; i < totalCols; i++)
    for (int j = 0; j < (int) (colIdxVector[i+1]-colIdxVector[i]); j++, rind++)
      {
        tdm->get(output_row, 0) = rowIdxVector[rind] + 1;
        tdm->get(output_row, 1) = i + 1;
        tdm->get(output_row, 2) = ptr[rind];
        output_row++;
      }

  /* If there are less elements than Nzmax (that might happen if some
     derivative is symbolically not zero but numerically zero at the evaluation
     point), then fill in the matrix with empty entries, that will be
     recognized as such by KordpDynare::populateDerivativesContainer() */
  while (output_row < (int) sizeRowIdxVector)
    {
      tdm->get(output_row, 0) = 0;
      tdm->get(output_row, 1) = 0;
      tdm->get(output_row, 2) = 0;
      output_row++;
    }
}
