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

// Full symmetry tensor.

/* Here we define folded and unfolded tensors for full symmetry. All
   tensors from here are identifying the multidimensional index with
   columns. */

#ifndef FS_TENSOR_H
#define FS_TENSOR_H

#include "tensor.hh"
#include "symmetry.hh"
#include "int_power.hh"

class FGSTensor;
class UGSTensor;
class FRSingleTensor;
class FSSparseTensor;

/* Folded tensor with full symmetry maintains only information about
   number of symmetrical variables nv. Further, we implement what is
   left from the super class FTensor.

   We implement getOffset() which should be used with care since
   its complexity.

   We implement a method adding a given general symmetry tensor to the
   full symmetry tensor supposing the variables of the general symmetry
   tensor are stacked giving only one variable of the full symmetry
   tensor.             ⎛y⎞
   For instance, if x= ⎝u⎠, then we can add tensor [g_y²u] to tensor [g_x³].
   This is done in method addSubTensor(). Consult FGSTensor class declaration
   to know what is general symmetry tensor.

   Note that the past-the-end index is of the form (nv,…,nv), because
   of the specific implementation of FFSTensor::increment().
*/

class UFSTensor;
class FFSTensor : public FTensor
{
  int nv;
public:
  /* Constructs given the number of rows (explicit since the tensor is
     column-oriented), the number of variables in each dimension, and the
     number of dimensions */
  FFSTensor(int r, int nvar, int d)
    : FTensor(indor::along_col, IntSequence(d, nvar),
              r, calcMaxOffset(nvar, d), d),
      nv(nvar)
  {
  }

  /* Constructs a tensor by one-dimensional contraction from the higher
     dimensional tensor t. This is, it constructs a tensor

     [g_yⁿ]_α₁…αₙ = [t_yⁿ⁺¹]_α₁…αₙβ [x]^β

     See the implementation for details. */
  FFSTensor(const FFSTensor &t, const ConstVector &x);

  /* Converts from sparse tensor (which is fully symmetric and folded by
     nature). */
  explicit FFSTensor(const FSSparseTensor &t);

  FFSTensor(const FFSTensor &) = default;
  FFSTensor(FFSTensor &&) = default;

  // Constructs from unfolded fully symmetric
  explicit FFSTensor(const UFSTensor &ut);

  // Constructs a subtensor of selected rows
  FFSTensor(int first_row, int num, FFSTensor &t)
    : FTensor(first_row, num, t), nv(t.nv)
  {
  }

  void increment(IntSequence &v) const override;
  void decrement(IntSequence &v) const override;
  std::unique_ptr<UTensor> unfold() const override;
  Symmetry
  getSym() const
  {
    return Symmetry{dimen()};
  }

  int getOffset(const IntSequence &v) const override;
  void addSubTensor(const FGSTensor &t);
  int
  nvar() const
  {
    return nv;
  }
  static int calcMaxOffset(int nvar, int d);
};

/* Unfolded fully symmetric tensor is almost the same in structure as
   FFSTensor, but the method unfoldData(). It takes columns which also
   exist in folded version and copies them to all their symmetrical
   locations. This is useful when constructing unfolded tensor from
   folded one. */

class UFSTensor : public UTensor
{
  int nv;
public:
  UFSTensor(int r, int nvar, int d)
    : UTensor(indor::along_col, IntSequence(d, nvar),
              r, calcMaxOffset(nvar, d), d),
      nv(nvar)
  {
  }
  UFSTensor(const UFSTensor &t, const ConstVector &x);
  UFSTensor(const UFSTensor &) = default;
  UFSTensor(UFSTensor &&) = default;
  explicit UFSTensor(const FFSTensor &ft);
  UFSTensor(int first_row, int num, UFSTensor &t)
    : UTensor(first_row, num, t), nv(t.nv)
  {
  }

  void increment(IntSequence &v) const override;
  void decrement(IntSequence &v) const override;
  std::unique_ptr<FTensor> fold() const override;
  Symmetry
  getSym() const
  {
    return Symmetry{dimen()};
  }

  int getOffset(const IntSequence &v) const override;
  void addSubTensor(const UGSTensor &t);
  int
  nvar() const
  {
    return nv;
  }
  static int
  calcMaxOffset(int nvar, int d)
  {
    return power(nvar, d);
  }
private:
  void unfoldData();
};

#endif
