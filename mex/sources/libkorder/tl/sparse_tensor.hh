/*
 * Copyright © 2004 Ondra Kamenik
 * Copyright © 2019-2023 Dynare Team
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

// Sparse tensor.

/* Here we declare a sparse full and general symmetry tensors with the
   multidimensional index along columns. We implement them as a std::multimap
   associating to each sequence of coordinates IntSequence a set of pairs (row,
   number). This is very convenient but not optimal in terms of memory
   consumption. So the implementation can be changed.

   The current std::multimap implementation allows insertions. Another
   advantage of this approach is that we do not need to calculate column
   numbers from the IntSequence, since the column is accessed directly via the
   key which is IntSequence.

   The only operation we need to do with the full symmetry sparse tensor
   is a left multiplication of a row oriented single column tensor. The
   result of such an operation is a column of the same size as the sparse
   tensor. Other important operations are slicing operations. We need to
   do sparse and dense slices of full symmetry sparse tensors. In fact,
   the only constructor of general symmetry sparse tensor is slicing from
   the full symmetry sparse. */

#ifndef SPARSE_TENSOR_H
#define SPARSE_TENSOR_H

#include "Vector.hh"
#include "gs_tensor.hh"
#include "symmetry.hh"
#include "tensor.hh"

#include <map>

struct ltseq
{
  bool
  operator()(const IntSequence& s1, const IntSequence& s2) const
  {
    return s1 < s2;
  }
};

/* This is a super class of both full symmetry and general symmetry sparse
   tensors. It contains a std::multimap and implements insertions. It tracks
   maximum and minimum row, for which there is an item. */

class SparseTensor
{
public:
  using Map = std::multimap<IntSequence, std::pair<int, double>, ltseq>;

protected:
  Map m;
  int dim;
  int nr;
  int nc;
  int first_nz_row;
  int last_nz_row {-1};

public:
  SparseTensor(int d, int nnr, int nnc) : dim(d), nr(nnr), nc(nnc), first_nz_row(nr)
  {
  }
  virtual ~SparseTensor() = default;
  void insert(IntSequence s, int r, double c);
  [[nodiscard]] const Map&
  getMap() const
  {
    return m;
  }
  [[nodiscard]] int
  dimen() const
  {
    return dim;
  }
  [[nodiscard]] int
  nrows() const
  {
    return nr;
  }
  [[nodiscard]] int
  ncols() const
  {
    return nc;
  }
  [[nodiscard]] double
  getFillFactor() const
  {
    return static_cast<double>(m.size()) / nrows() / ncols();
  }
  [[nodiscard]] double getFoldIndexFillFactor() const;
  [[nodiscard]] double getUnfoldIndexFillFactor() const;
  [[nodiscard]] int
  getNumNonZero() const
  {
    return m.size();
  }
  [[nodiscard]] int
  getFirstNonZeroRow() const
  {
    return first_nz_row;
  }
  [[nodiscard]] int
  getLastNonZeroRow() const
  {
    return last_nz_row;
  }
  [[nodiscard]] virtual const Symmetry& getSym() const = 0;
  void print() const;
  [[nodiscard]] bool isFinite() const;
};

/* This is a full symmetry sparse tensor. It implements multColumnAndAdd() and,
   in addition to SparseTensor, it has ‘nv’ (number of variables) and symmetry
   (basically it is a dimension). */

class FSSparseTensor : public SparseTensor
{
private:
  int nv;
  Symmetry sym;

public:
  FSSparseTensor(int d, int nvar, int r);
  void insert(IntSequence s, int r, double c);
  void multColumnAndAdd(const Tensor& t, Vector& v) const;
  [[nodiscard]] const Symmetry&
  getSym() const override
  {
    return sym;
  }
  [[nodiscard]] int
  nvar() const
  {
    return nv;
  }
  void print() const;
};

/* This is a general symmetry sparse tensor. It has TensorDimens and can be
   constructed as a slice of the full symmetry sparse tensor. The slicing
   constructor takes the same form as the slicing FGSTensor constructor from
   full symmetry sparse tensor. */

class GSSparseTensor : public SparseTensor
{
private:
  TensorDimens tdims;

public:
  GSSparseTensor(const FSSparseTensor& t, const IntSequence& ss, const IntSequence& coor,
                 TensorDimens td);
  void insert(IntSequence s, int r, double c);
  [[nodiscard]] const Symmetry&
  getSym() const override
  {
    return tdims.getSym();
  }
  [[nodiscard]] const TensorDimens&
  getDims() const
  {
    return tdims;
  }
  void print() const;
};

#endif
