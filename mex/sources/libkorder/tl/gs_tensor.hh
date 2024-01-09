/*
 * Copyright © 2004 Ondra Kamenik
 * Copyright © 2019-2024 Dynare Team
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

// General symmetry tensor.

/* Here we define tensors for general symmetry. All tensors from here are
   identifying the multidimensional index with columns. Thus all
   symmetries regard to columns. The general symmetry here is not the most
   general. It captures all symmetries of indices which are given by
   continuous partitioning of indices. Two items are symmetric if they
   belong to the same group. The continuity implies that if two items
   belong to one group, then all items between them belong to that
   group. This continuous partitioning of indices is described by
   Symmetry class.

   The dimension of the tensors here are described (besides the symmetry)
   also by number of variables for each group. This is dealt in the class
   for tensor dimensions defined also here. */

#ifndef GS_TENSOR_HH
#define GS_TENSOR_HH

#include "fs_tensor.hh"
#include "rfs_tensor.hh"
#include "symmetry.hh"
#include "tensor.hh"

class FGSTensor;
class UGSTensor;
class FSSparseTensor;

/* This class encapsulates symmetry information for the general
   symmetry tensor. It maintains a vector of variable numbers ‘nvs’, and
   symmetry ‘sym’. For example, let the symmetry be y²u³, and
   variable numbers be 10 for y, and 5 for u. Then the ‘nvs’ is
   (10,5), and ‘sym’ is (2,3). Also it maintains ‘nvmax’ unfolded ‘nvs’ with
   respect to the symmetry, this is (10,10,5,5,5).

   The class is able to calculate number of offsets (columns or rows depending
   what matrix coordinate we describe) in unfolded and folded tensors
   with the given symmetry. */

class TensorDimens
{
protected:
  IntSequence nvs;
  Symmetry sym;
  IntSequence nvmax;

public:
  TensorDimens(Symmetry s, IntSequence nvars) :
      nvs(std::move(nvars)), sym(std::move(s)), nvmax(nvs.unfold(sym))
  {
  }
  // Full-symmetry special case
  TensorDimens(int nvar, int dimen) : nvs {nvar}, sym {dimen}, nvmax(dimen, nvar)
  {
  }
  // Constructs the tensor dimensions for slicing (see the implementation for details)
  TensorDimens(const IntSequence& ss, const IntSequence& coor);

  [[nodiscard]] bool
  operator==(const TensorDimens& td) const
  {
    return nvs == td.nvs && sym == td.sym;
  }

  [[nodiscard]] int
  dimen() const
  {
    return sym.dimen();
  }
  [[nodiscard]] int
  getNVX(int i) const
  {
    return nvmax[i];
  }
  [[nodiscard]] const IntSequence&
  getNVS() const
  {
    return nvs;
  }
  [[nodiscard]] const IntSequence&
  getNVX() const
  {
    return nvmax;
  }
  [[nodiscard]] const Symmetry&
  getSym() const
  {
    return sym;
  }

  [[nodiscard]] int calcUnfoldMaxOffset() const;
  [[nodiscard]] int calcFoldMaxOffset() const;
  [[nodiscard]] int calcFoldOffset(const IntSequence& v) const;
  void decrement(IntSequence& v) const;
};

/* Here is a class for folded general symmetry tensor. It only contains
   tensor dimensions, it defines types for indices, implement virtual
   methods of super class FTensor. */

class GSSparseTensor;
class FGSTensor : public FTensor
{
  friend class UGSTensor;

  const TensorDimens tdims;

public:
  FGSTensor(int r, TensorDimens td) :
      FTensor(indor::along_col, td.getNVX(), r, td.calcFoldMaxOffset(), td.dimen()),
      tdims(std::move(td))
  {
  }
  FGSTensor(const FGSTensor&) = default;
  FGSTensor(FGSTensor&&) = default;

  FGSTensor(int first_row, int num, FGSTensor& t) : FTensor(first_row, num, t), tdims(t.tdims)
  {
  }

  // Constructs a slice from a fully symmetric sparse tensor
  FGSTensor(const FSSparseTensor& t, const IntSequence& ss, const IntSequence& coor,
            TensorDimens td);

  // Constructs a slice from a fully symmetric dense tensor
  FGSTensor(const FFSTensor& t, const IntSequence& ss, const IntSequence& coor, TensorDimens td);

  // Converting constructors
  explicit FGSTensor(const UGSTensor& ut);
  explicit FGSTensor(const GSSparseTensor& sp);
  explicit FGSTensor(FFSTensor& t) : FTensor(0, t.nrows(), t), tdims(t.nvar(), t.dimen())
  {
  }

  ~FGSTensor() override = default;

  void increment(IntSequence& v) const override;
  void
  decrement(IntSequence& v) const override
  {
    tdims.decrement(v);
  }
  [[nodiscard]] std::unique_ptr<UTensor> unfold() const override;
  [[nodiscard]] const TensorDimens&
  getDims() const
  {
    return tdims;
  }
  [[nodiscard]] const Symmetry&
  getSym() const
  {
    return getDims().getSym();
  }

  /* Performs a contraction of one variable in the tensor. This is, for
     instance:

      [r_xⁱzᵏ]_α₁…αᵢγ₁…γₖ = [t_xⁱyʲzᵏ]_α₁…αᵢβ₁…βⱼγ₁…γₖ·[c]^β₁…βⱼ
  */
  void contractAndAdd(int i, FGSTensor& out, const FRSingleTensor& col) const;

  [[nodiscard]] int
  getOffset(const IntSequence& v) const override
  {
    return tdims.calcFoldOffset(v);
  }
};

/* Besides similar things that has FGSTensor, we have here also
   method unfoldData(), and helper method getFirstIndexOf()
   which corresponds to sorting coordinates in fully symmetric case (here
   the action is more complicated, so we put it to the method). */

class UGSTensor : public UTensor
{
  friend class FGSTensor;

  const TensorDimens tdims;

public:
  UGSTensor(int r, TensorDimens td) :
      UTensor(indor::along_col, td.getNVX(), r, td.calcUnfoldMaxOffset(), td.dimen()),
      tdims(std::move(td))
  {
  }

  UGSTensor(const UGSTensor&) = default;
  UGSTensor(UGSTensor&&) = default;

  UGSTensor(int first_row, int num, UGSTensor& t) : UTensor(first_row, num, t), tdims(t.tdims)
  {
  }

  // Constructs a slice from fully symmetric sparse tensor
  UGSTensor(const FSSparseTensor& t, const IntSequence& ss, const IntSequence& coor,
            TensorDimens td);

  // Constructs a slice from fully symmetric dense unfolded tensor
  UGSTensor(const UFSTensor& t, const IntSequence& ss, const IntSequence& coor, TensorDimens td);

  // Converting constructors
  explicit UGSTensor(const FGSTensor& ft);
  explicit UGSTensor(UFSTensor& t) : UTensor(0, t.nrows(), t), tdims(t.nvar(), t.dimen())
  {
  }

  ~UGSTensor() override = default;

  void increment(IntSequence& v) const override;
  void decrement(IntSequence& v) const override;
  [[nodiscard]] std::unique_ptr<FTensor> fold() const override;
  [[nodiscard]] const TensorDimens&
  getDims() const
  {
    return tdims;
  }
  [[nodiscard]] const Symmetry&
  getSym() const
  {
    return getDims().getSym();
  }

  void contractAndAdd(int i, UGSTensor& out, const URSingleTensor& col) const;
  [[nodiscard]] int getOffset(const IntSequence& v) const override;

private:
  void unfoldData();

public:
  [[nodiscard]] index getFirstIndexOf(const index& in) const;
};

#endif
