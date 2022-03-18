/*
 * Copyright © 2010-2022 Dynare Team
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
#include <cassert>

#include "dynare_exception.hh"

#include "dynamic_m.hh"

DynamicModelMFile::DynamicModelMFile(const std::string &modName, int ntt_arg) :
  DynamicModelAC(ntt_arg),
  DynamicMFilename{modName + ".dynamic"}
{
}

/* NB: This is a duplicate of DynamicModelMatlabCaller::cmplxToReal() in
   perfect_foresight_problem MEX */
mxArray *
DynamicModelMFile::cmplxToReal(mxArray *cmplx_mx)
{
  mxArray *real_mx = mxCreateDoubleMatrix(mxGetM(cmplx_mx), mxGetN(cmplx_mx), mxREAL);

#if MX_HAS_INTERLEAVED_COMPLEX
  mxComplexDouble *cmplx = mxGetComplexDoubles(cmplx_mx);
#else
  double *cmplx_real = mxGetPr(cmplx_mx);
  double *cmplx_imag = mxGetPi(cmplx_mx);
#endif
  double *real = mxGetPr(real_mx);

  for (size_t i = 0; i < mxGetNumberOfElements(cmplx_mx); i++)
#if MX_HAS_INTERLEAVED_COMPLEX
    if (cmplx[i].imag == 0.0)
      real[i] = cmplx[i].real;
#else
    if (cmplx_imag[i] == 0.0)
      real[i] = cmplx_real[i];
#endif
    else
      real[i] = std::numeric_limits<double>::quiet_NaN();

  mxDestroyArray(cmplx_mx);
  return real_mx;
}

void
DynamicModelMFile::unpackSparseMatrixAndCopyIntoTwoDMatData(mxArray *sparseMat, TwoDMatrix &tdm)
{
  int totalCols = mxGetN(sparseMat);
  mwIndex *rowIdxVector = mxGetIr(sparseMat);
  mwIndex *colIdxVector = mxGetJc(sparseMat);

  assert(tdm.ncols() == 3);
  /* Under MATLAB, the following check always holds at equality; under Octave,
     there may be an inequality, because Octave diminishes nzmax if one gives
     zeros in the values vector when calling sparse(). */
  assert(tdm.nrows() >= mxGetNzmax(sparseMat));

  int rind = 0;
  int output_row = 0;

  for (int i = 0; i < totalCols; i++)
    for (int j = 0; j < static_cast<int>((colIdxVector[i+1]-colIdxVector[i])); j++, rind++)
      {
        tdm.get(output_row, 0) = rowIdxVector[rind] + 1;
        tdm.get(output_row, 1) = i + 1;
        if (!mxIsComplex(sparseMat))
          tdm.get(output_row, 2) = mxGetPr(sparseMat)[rind];
        else
          {
            double real, imag;
#if MX_HAS_INTERLEAVED_COMPLEX
            mxComplexDouble cmplx = mxGetComplexDoubles(sparseMat)[rind];
            real = cmplx.real;
            imag = cmplx.imag;
#else
            real = mxGetPr(sparseMat)[rind];
            imag = mxGetPi(sparseMat)[rind];
#endif
            tdm.get(output_row, 2) = imag == 0.0 ? real : std::numeric_limits<double>::quiet_NaN();
          }
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
DynamicModelMFile::eval(const Vector &y, const Vector &x, const Vector &modParams, const Vector &ySteady,
                        Vector &residual, std::vector<TwoDMatrix> &md) noexcept(false)
{
  mxArray *T_m = mxCreateDoubleMatrix(ntt, 1, mxREAL);

  mxArray *y_m = mxCreateDoubleMatrix(y.length(), 1, mxREAL);
  std::copy_n(y.base(), y.length(), mxGetPr(y_m));

  mxArray *x_m = mxCreateDoubleMatrix(1, x.length(), mxREAL);
  std::copy_n(x.base(), x.length(), mxGetPr(x_m));

  mxArray *params_m = mxCreateDoubleMatrix(modParams.length(), 1, mxREAL);
  std::copy_n(modParams.base(), modParams.length(), mxGetPr(params_m));

  mxArray *steady_state_m = mxCreateDoubleMatrix(ySteady.length(), 1, mxREAL);
  std::copy_n(ySteady.base(), ySteady.length(), mxGetPr(steady_state_m));

  mxArray *it_m = mxCreateDoubleScalar(1.0);
  mxArray *T_flag_m = mxCreateLogicalScalar(false);

  {
    // Compute temporary terms (for all orders)
    std::string funcname = DynamicMFilename + "_g" + std::to_string(md.size()) + "_tt";
    mxArray *plhs[1], *prhs[] = { T_m, y_m, x_m, params_m, steady_state_m, it_m };

    int retVal = mexCallMATLAB(1, plhs, 6, prhs, funcname.c_str());
    if (retVal != 0)
      throw DynareException(__FILE__, __LINE__, "Trouble calling " + funcname);

    mxDestroyArray(T_m);
    T_m = plhs[0];
  }

  {
    // Compute residuals
    std::string funcname = DynamicMFilename + "_resid";
    mxArray *plhs[1], *prhs[] = { T_m, y_m, x_m, params_m, steady_state_m, it_m, T_flag_m };

    int retVal = mexCallMATLAB(1, plhs, 7, prhs, funcname.c_str());
    if (retVal != 0)
      throw DynareException(__FILE__, __LINE__, "Trouble calling " + funcname);

    if (!mxIsDouble(plhs[0]) || mxIsSparse(plhs[0]))
      throw DynareException(__FILE__, __LINE__, "Residual should be a dense array of double floats");

    if (mxIsComplex(plhs[0]))
      plhs[0] = cmplxToReal(plhs[0]);

    residual = Vector{plhs[0]};
    mxDestroyArray(plhs[0]);
  }

  for (size_t i = 1; i <= md.size(); i++)
    {
      // Compute model derivatives
      std::string funcname = DynamicMFilename + "_g" + std::to_string(i);
      mxArray *plhs[1], *prhs[] = { T_m, y_m, x_m, params_m, steady_state_m, it_m, T_flag_m };

      int retVal = mexCallMATLAB(1, plhs, 7, prhs, funcname.c_str());
      if (retVal != 0)
        throw DynareException(__FILE__, __LINE__, "Trouble calling " + funcname);

      if (!mxIsDouble(plhs[0]))
        throw DynareException(__FILE__, __LINE__, "Derivatives matrix at order " + std::to_string(i) + "should be an array of double floats");

      if (i == 1)
        {
          if (mxIsSparse(plhs[0]))
            throw DynareException(__FILE__, __LINE__, "Derivatives matrix at order " + std::to_string(i) + " should be dense");
          assert(static_cast<int>(mxGetM(plhs[0])) == md[i-1].nrows());
          assert(static_cast<int>(mxGetN(plhs[0])) == md[i-1].ncols());
          std::copy_n(mxGetPr(plhs[0]), mxGetM(plhs[0])*mxGetN(plhs[0]), md[i-1].base());
        }
      else
        {
          if (!mxIsSparse(plhs[0]))
            throw DynareException(__FILE__, __LINE__, "Derivatives matrix at order " + std::to_string(i) + " should be sparse");
          unpackSparseMatrixAndCopyIntoTwoDMatData(plhs[0], md[i-1]);
        }

      mxDestroyArray(plhs[0]);
    }

  mxDestroyArray(T_m);
  mxDestroyArray(y_m);
  mxDestroyArray(x_m);
  mxDestroyArray(params_m);
  mxDestroyArray(steady_state_m);
  mxDestroyArray(it_m);
  mxDestroyArray(T_flag_m);
}
