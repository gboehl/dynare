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

// Multiplying stacked tensor columns.

/* We need to calculate the following tensor product:
                   ⱼ                      ₗ
   [f_sʲ]_α₁…αⱼ =  ∑ [f_zˡ]_β₁…βₗ   ∑     ∏  [z_{s^|cₘ|}]_cₘ(α)^βₘ
                  ˡ⁼¹            c∈ℳₗ,ₖ ᵐ⁼¹
   where s=(y,u,u′,σ), and z is a composition of four variables, say (v,w,y,u).
   Note that z ends with y and u, and the only non-zero derivative of the
   trailing part of z involving y or u is the first derivative and is the unit
   matrix y_y=(1) or uᵤ=(1). Also, we suppose that the dependence of v, and w
   on s is such that whenever derivative of w is nonzero, then also of v. This
   means that there for any derivative and any index there is a continuous part
   of derivatives of v and optionally of w followed by column of zeros
   containing at most one 1.

   This structure can be modelled and exploited with some costs at
   programming. For example, let us consider the following product:

    [B_y²u³]_α₁α₂β₁β₂β₃ = … + [f_z³]_γ₁γ₂γ₃ [z_yu]^γ₁_α₁β₁ [z_y]^γ₂_α₂ [zᵤᵤ]^γ₃_β₂β₃ + …

   The term corresponds to equivalence { {0,2}, {1}, {3,4} }. For the fixed
   index α₁α₂β₁β₂β₃ we have to make a Kronecker product of the columns

    [z_yu]_α₁β₁ ⊗ [z_y]_α₂ ⊗ [zᵤᵤ]_β₂β₃

   which can be written as:

    ⎛[v_yu]_α₁β₁⎞ ⎛[v_y]_α₂⎞ ⎛[vᵤᵤ]_β₂β₃⎞
    ⎢[w_yu]_α₁β₁⎥ ⎢[w_y]_α₂⎥ ⎢[wᵤᵤ]_β₂β₃⎥
    ⎢     0     ⎥⊗⎢  1_α₂  ⎥⊗⎢     0    ⎥
    ⎝     0     ⎠ ⎝    0   ⎠ ⎝     0    ⎠
   where 1_α₂ is a column of zeros having the only 1 at α₂ index.

   This file develops the abstraction for this Kronecker product column
   without multiplication of the zeros at the top. Basically, it will be
   a column which is a Kronecker product of the columns without the
   zeros:

                  ⎛[v_y]_α₂⎞
   ⎛[v_yu]_α₁β₁⎞⊗ ⎢[w_y]_α₂⎥⊗⎛[vᵤᵤ]_β₂β₃⎞
   ⎝[w_yu]_α₁β₁⎠  ⎝    1   ⎠ ⎝[wᵤᵤ]_β₂β₃⎠

   The class will have a tensor infrastructure introducing ‘index’ which
   iterates over all items in the column with γ₁γ₂γ₃
   as coordinates in [f_z³]. The data of such a tensor is
   not suitable for any matrix operation and will have to be accessed
   only through the ‘index’. Note that this does not matter, since
   [f_zˡ] are sparse. */

#ifndef PYRAMID_PROD2_HH
#define PYRAMID_PROD2_HH

#include "permutation.hh"
#include "rfs_tensor.hh"
#include "stack_container.hh"
#include "tensor.hh"
#include "tl_exception.hh"

#include "Vector.hh"

#include <vector>

/* First we declare a helper class for the tensor. Its purpose is to
   gather the columns which are going to be Kronecker multiplied. The
   input of this helper class is StackProduct<FGSTensor> and coordinate
   ‘c’ of the column.

   It maintains ‘unit_flag’ array which says for what columns we must
   stack 1 below v and w. In this case, the value of ‘unit_flag’ is
   an index of the 1, otherwise the value of ‘unit_flag’ is -1.

   Also we have storage for the stacked columns ‘cols’. The object is
   responsible for memory management associated to this storage. That is
   why we do not allow any copy constructor, since we need to be sure
   that no accidental copies take place. We declare the copy constructor
   as private and not implement it. */

class IrregTensor;
class IrregTensorHeader
{
  friend class IrregTensor;
  int nv;
  IntSequence unit_flag;
  std::vector<std::unique_ptr<Vector>> cols;
  IntSequence end_seq;

public:
  IrregTensorHeader(const StackProduct<FGSTensor>& sp, const IntSequence& c);
  [[nodiscard]] int
  dimen() const
  {
    return unit_flag.size();
  }
  void increment(IntSequence& v) const;
  [[nodiscard]] int calcMaxOffset() const;
};

/* Here we declare the irregular tensor. There is no special logic here. We
   inherit from Tensor and we must implement three methods, increment(),
   decrement() and getOffset(). The last two are not implemented now, since
   they are not needed, and they raise an exception. The first just calls
   increment() of the header. Also we declare a method addTo() which adds this
   unfolded irregular single column tensor to folded (regular) single column
   tensor.

   The header IrregTensorHeader lives with an object by a reference. This is
   dangerous. However, we will use this class only in a simple loop and both
   IrregTensor and IrregTensorHeader will be destructed at the end of a block.
   Since the super class Tensor must be initialized before any member, we could
   do either a save copy of IrregTensorHeader, or relatively dangerous the
   reference member. For the reason above we chose the latter. */

class IrregTensor : public Tensor
{
  const IrregTensorHeader& header;

public:
  IrregTensor(const IrregTensorHeader& h);
  void addTo(FRSingleTensor& out) const;
  void
  increment(IntSequence& v) const override
  {
    header.increment(v);
  }
  void
  decrement([[maybe_unused]] IntSequence& v) const override
  {
    TL_RAISE("Not implemented error in IrregTensor::decrement");
  }
  [[nodiscard]] int
  getOffset([[maybe_unused]] const IntSequence& v) const override
  {
    TL_RAISE("Not implemented error in IrregTensor::getOffset");
  }
};

#endif
