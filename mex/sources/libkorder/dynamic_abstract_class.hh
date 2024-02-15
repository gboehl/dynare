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

#ifndef DYNAMIC_ABSTRACT_CLASS_HH
#define DYNAMIC_ABSTRACT_CLASS_HH

#include <vector>

#include "dynmex.h"

#include "Vector.hh"
#include "sparse_tensor.hh"
#include "t_container.hh"

class DynamicModelAC
{
protected:
  const int order;
  const mxArray *const dynamic_g1_sparse_rowval_mx, *const dynamic_g1_sparse_colval_mx,
                                                        *const dynamic_g1_sparse_colptr_mx;
  // Stores M_.dynamic_gN_sparse_indices, starting from N=2
  const std::vector<const mxArray*> dynamic_gN_sparse_indices;

public:
  DynamicModelAC(int order_arg, const mxArray* dynamic_g1_sparse_rowval_mx_arg,
                 const mxArray* dynamic_g1_sparse_colval_mx_arg,
                 const mxArray* dynamic_g1_sparse_colptr_mx_arg,
                 const std::vector<const mxArray*> dynamic_gN_sparse_indices_arg) :
      order {order_arg},
      dynamic_g1_sparse_rowval_mx {dynamic_g1_sparse_rowval_mx_arg},
      dynamic_g1_sparse_colval_mx {dynamic_g1_sparse_colval_mx_arg},
      dynamic_g1_sparse_colptr_mx {dynamic_g1_sparse_colptr_mx_arg},
      dynamic_gN_sparse_indices {dynamic_gN_sparse_indices_arg} {};
  virtual ~DynamicModelAC() = default;
  virtual void eval(const Vector& y, const Vector& x, const Vector& params, const Vector& ySteady,
                    Vector& residual, const std::map<int, int>& dynToDynpp,
                    TensorContainer<FSSparseTensor>& derivatives) noexcept(false)
      = 0;
};

#endif
