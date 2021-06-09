/*
 * Copyright © 2004 Ondra Kamenik
 * Copyright © 2019 Dynare Team
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

#ifndef FACTORY_H
#define FACTORY_H

#include <random>
#include <memory>

#include "symmetry.hh"
#include "int_sequence.hh"
#include "twod_matrix.hh"
#include "equivalence.hh"
#include "rfs_tensor.hh"
#include "t_container.hh"

class Factory
{
  std::mt19937 mtgen;
  std::uniform_real_distribution<> dis;

  void init(const Symmetry &s, const IntSequence &nvs);
  void init(int dim, int nv);
  void fillMatrix(TwoDMatrix &m);
public:
  double get();
  // This can be used with UGSTensor, FGSTensor
  template<class _Ttype>
  std::unique_ptr<_Ttype>
  make(int r, const Symmetry &s, const IntSequence &nvs)
  {
    auto res = std::make_unique<_Ttype>(r, TensorDimens(s, nvs));
    init(s, nvs);
    fillMatrix(*res);
    return res;
  }

  // This can be used with FFSTensor, UFSTensor, FRTensor, URTensor
  template<class _Ttype>
  std::unique_ptr<_Ttype>
  make(int r, int nv, int dim)
  {
    auto res = std::make_unique<_Ttype>(r, nv, dim);
    init(dim, nv);
    fillMatrix(*res);
    return res;
  }

  template<class _Ttype, class _Ctype>
  _Ctype
  makeCont(int r, const IntSequence &nvs, int maxdim)
  {
    int symnum = nvs.size();
    _Ctype res(symnum);
    for (int dim = 1; dim <= maxdim; dim++)
      if (symnum == 1)
        // Full symmetry
        res.insert(make<_Ttype>(r, Symmetry{dim}, nvs));
      else
        // General symmetry
        for (int i = 0; i <= dim; i++)
          res.insert(make<_Ttype>(r, Symmetry{i, dim-i}, nvs));
    return res;
  }

  template<class _Ttype, class _Ptype>
  _Ptype
  makePoly(int r, int nv, int maxdim)
  {
    _Ptype p(r, nv);
    for (int d = 1; d <= maxdim; d++)
      p.insert(make<_Ttype>(r, nv, d));
    return p;
  }

  Vector makeVector(int n);
};

#endif
