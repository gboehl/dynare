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

#ifndef OBJECTIVE_M_HH
#define OBJECTIVE_M_HH

#include <vector>

#include "dynmex.h"

#include "Vector.hh"
#include "sparse_tensor.hh"
#include "t_container.hh"

// Handles calls to <model>/+objective/+sparse/static*.m
class ObjectiveMFile
{
private:
  const std::string ObjectiveMFilename;
  const int kOrder;
  const mxArray *const objective_g1_sparse_rowval_mx, *const objective_g1_sparse_colval_mx,
                                                          *const objective_g1_sparse_colptr_mx;
  // Stores M_.objective_gN_sparse_indices, starting from N=2
  const std::vector<const mxArray*> objective_gN_sparse_indices;

public:
  ObjectiveMFile(const std::string& modName, int kOrder_arg,
                 const mxArray* objective_g1_sparse_rowval_mx_arg,
                 const mxArray* objective_g1_sparse_colval_mx_arg,
                 const mxArray* objective_g1_sparse_colptr_mx_arg,
                 const std::vector<const mxArray*> objective_gN_sparse_indices_arg);
  void eval(const Vector& y, const Vector& x, const Vector& params, Vector& residual,
            const std::vector<int>& dynToDynpp, TensorContainer<FSSparseTensor>& derivatives) const
      noexcept(false);
};

#endif
