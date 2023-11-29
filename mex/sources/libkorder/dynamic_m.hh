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

#ifndef _DYNAMIC_M_HH
#define _DYNAMIC_M_HH

#include "dynamic_abstract_class.hh"

#include "mex.h"
#include <dynmex.h>

/**
 * handles calls to <model>_dynamic.m
 *
 **/
class DynamicModelMFile : public DynamicModelAC
{
private:
  const std::string DynamicMFilename;
  /* Unpack a sparse matrix (of double floats) into a TwoDMatrix object.
     Real elements of the original matrix are copied as-is to the new one.
     Complex elements are replaced by NaNs. */
  static void unpackSparseMatrixAndCopyIntoTwoDMatData(mxArray* sparseMat, TwoDMatrix& tdm);
  /* Given a complex dense matrix (of double floats), returns a real dense matrix of same size.
     Real elements of the original matrix are copied as-is to the new one.
     Complex elements are replaced by NaNs.
     Destroys the original matrix. */
  static mxArray* cmplxToReal(mxArray* m);

public:
  explicit DynamicModelMFile(const std::string& modName, int ntt_arg);
  virtual ~DynamicModelMFile() = default;
  void eval(const Vector& y, const Vector& x, const Vector& params, const Vector& ySteady,
            Vector& residual, std::vector<TwoDMatrix>& md) override;
};
#endif
