/*
 * Copyright © 2010-2024 Dynare Team
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

#ifndef DYNAMIC_M_HH
#define DYNAMIC_M_HH

#include "dynamic_abstract_class.hh"

/**
 * handles calls to <model>_dynamic.m
 *
 **/
class DynamicModelMFile : public DynamicModelAC
{
private:
  const std::string DynamicMFilename;
  /* Given a complex dense matrix (of double floats), returns a real dense matrix of same size.
     Real elements of the original matrix are copied as-is to the new one.
     Complex elements are replaced by NaNs.
     Destroys the original matrix. */
  static mxArray* cmplxToReal(mxArray* m);

public:
  DynamicModelMFile(const std::string& modName, int order_arg,
                    const mxArray* dynamic_g1_sparse_rowval_mx_arg,
                    const mxArray* dynamic_g1_sparse_colval_mx_arg,
                    const mxArray* dynamic_g1_sparse_colptr_mx_arg,
                    const std::vector<const mxArray*> dynamic_gN_sparse_indices_arg);
  void eval(const Vector& y, const Vector& x, const Vector& params, const Vector& ySteady,
            Vector& residual, const std::map<int, int>& dynToDynpp,
            TensorContainer<FSSparseTensor>& derivatives) noexcept(false) override;
};

#endif
