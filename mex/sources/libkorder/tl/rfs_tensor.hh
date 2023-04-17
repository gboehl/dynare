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

// Row-wise full symmetry tensor.

/* Here we define classes for full symmetry tensors with the multidimensional
   index identified with rows. The primary usage is for storage of data coming
   from (or from a sum of)
      ₗ
      ∏  [g_{s^|cₘ|}]_cₘ(α)^γₘ
     ᵐ⁼¹
   where α comes from a multidimensional index that goes through some set S,
   and cₘ is some equivalence class. So we model a tensor of the form:

     ⎡ ₗ                       ⎤
     ⎢ ∏  [g_{s^|cₘ|}]_cₘ(α)^γₘ⎥
     ⎣ᵐ⁼¹                      ⎦S^γ₁…γₗ

   Since all γ₁…γₗ correspond to the same variable, the tensor is fully
   symmetric. The set of indices S cannot be very large and sometimes it is
   only one element. This case is handled in a special subclass.

   We provide both folded and unfolded versions. Their logic is perfectly the
   same as in UFSTensor and FFSTensor with two exceptions. One has been already
   mentioned, the multidimensional index is along the rows. The second are
   conversions between the two types. Since this kind of tensor is used to
   multiply (from the right) a tensor whose multidimensional index is
   identified with columns, we will need a different way of a conversion. If
   the multiplication of two folded tensors is to be equivalent with
   multiplication of two unfolded, the folding of the right tensor must sum all
   equivalent elements since they are multiplied with the same number from the
   folded tensor. (Equivalent here means all elements of unfolded tensor
   corresponding to one element in folded tensor.) For this reason, it is
   necessary to calculate a column number from the given sequence, so we
   implement getOffset(). Process of unfolding is not used, so we implemented
   it so that unfolding and then folding a tensor would yield the same data. */

#ifndef RFS_TENSOR_H
#define RFS_TENSOR_H

#include "tensor.hh"
#include "fs_tensor.hh"
#include "symmetry.hh"

/* This is straightforward and very similar to UFSTensor. */

class FRTensor;
class URTensor : public UTensor
{
  int nv;
public:
  URTensor(int c, int nvar, int d)
    : UTensor(indor::along_row, IntSequence(d, nvar),
              UFSTensor::calcMaxOffset(nvar, d), c, d), nv(nvar)
  {
  }
  URTensor(const URTensor &) = default;
  URTensor(URTensor &&) = default;
  explicit URTensor(const FRTensor &ft);

  ~URTensor() override = default;

  void increment(IntSequence &v) const override;
  void decrement(IntSequence &v) const override;
  std::unique_ptr<FTensor> fold() const override;

  int getOffset(const IntSequence &v) const override;
  int
  nvar() const
  {
    return nv;
  }
  Symmetry
  getSym() const
  {
    return Symmetry{dimen()};
  }
};

/* This is straightforward and very similar to FFSTensor. */

class FRTensor : public FTensor
{
  int nv;
public:
  FRTensor(int c, int nvar, int d)
    : FTensor(indor::along_row, IntSequence(d, nvar),
              FFSTensor::calcMaxOffset(nvar, d), c, d), nv(nvar)
  {
  }
  FRTensor(const FRTensor &) = default;
  FRTensor(FRTensor &&) = default;
  explicit FRTensor(const URTensor &ut);

  ~FRTensor() override = default;

  void increment(IntSequence &v) const override;
  void decrement(IntSequence &v) const override;
  std::unique_ptr<UTensor> unfold() const override;

  int
  nvar() const
  {
    return nv;
  }
  int
  getOffset(const IntSequence &v) const override
  {
    return FTensor::getOffset(v, nv);
  }
  Symmetry
  getSym() const
  {
    return Symmetry{dimen()};
  }
};

/* The following class represents specialization of URTensor coming from
   Kronecker multiplication of a few vectors. So the resulting row-oriented
   tensor has one column. We provide two constructors, one constructs the
   tensor from a few vectors stored as std::vector<ConstVector>. The second
   makes the Kronecker power of one given vector. */

class URSingleTensor : public URTensor
{
public:
  URSingleTensor(int nvar, int d)
    : URTensor(1, nvar, d)
  {
  }
  URSingleTensor(const std::vector<ConstVector> &cols);
  URSingleTensor(const ConstVector &v, int d);
  URSingleTensor(const URSingleTensor &) = default;
  URSingleTensor(URSingleTensor &&) = default;
  ~URSingleTensor() override = default;
  std::unique_ptr<FTensor> fold() const override;
};

/* This class represents one column row-oriented tensor. The only way to
   construct it is from URSingleTensor or from scratch. The folding algorithm
   is the same as folding of general URTensor. Only its implementation is
   different, since we do not copy rows, but only elements. */

class FRSingleTensor : public FRTensor
{
public:
  FRSingleTensor(int nvar, int d)
    : FRTensor(1, nvar, d)
  {
  }
  explicit FRSingleTensor(const URSingleTensor &ut);
  FRSingleTensor(const FRSingleTensor &) = default;
  FRSingleTensor(FRSingleTensor &&) = default;
  ~FRSingleTensor() override = default;
};

#endif
