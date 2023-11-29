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

// Multiplying tensor columns.

/* In here, we implement the Faà Di Bruno for folded
   tensors. Recall, that one step of the Faà Di Bruno is a formula:

                                 ₗ
    [B_sᵏ]_α₁…αₗ = [h_yˡ]_γ₁…γₗ  ∏  [g_{s^|cₘ|}]_cₘ(α)^γₘ
                                ᵐ⁼¹

   In contrast to unfolded implementation of UGSContainer::multAndAdd()
   with help of KronProdAll and UPSTensor, we take a completely
   different strategy. We cannot afford full instantiation of

             ₗ
       ∑     ∏  [g_{s^|cₘ|}]_cₘ(α)^γₘ
    c∈ℳₗ,ₖ ᵐ⁼¹

   and therefore we do it per partes. We select some number of columns,
   for instance 10, calculate 10 continuous iterators of tensor B. Then we
   form unfolded tensor

                  ⎡         ₗ                       ⎤
    [G]_S^γ₁…γₗ = ⎢   ∑     ∏  [g_{s^|cₘ|}]_cₘ(α)^γₘ⎥
                  ⎣c∈ℳₗ,ₖ ᵐ⁼¹                       ⎦_S

   where S is the selected set of 10 indices. This is done as Kronecker
   product of vectors corresponding to selected columns. Note that, in
   general, there is no symmetry in G, its type is special class for
   this purpose.

   If g is folded, then we have to form the folded version of G. There is
   no symmetry in G data, so we sum all unfolded indices corresponding
   to folded index together. This is perfectly OK, since we multiply
   these groups of (equivalent) items with the same number in fully
   symmetric g.

   After this, we perform ordinary matrix multiplication to obtain a
   selected set of columns of B.

   In here, we define a class for forming and representing [G]_S^γ₁…γₗ.
   Basically, this tensor is row-oriented (multidimensional index is along
   rows), and it is fully symmetric. So we inherit from URTensor. If we need
   its folded version, we simply use a suitable conversion. The new abstraction
   will have only a new constructor allowing a construction from the given set
   of indices S, and given set of tensors g. The rest of the process is
   implemented in FGSContainer::multAndAdd() unfolded code or
   FGSContainer::multAndAdd() folded code. */

#ifndef PYRAMID_PROD_H
#define PYRAMID_PROD_H

#include "gs_tensor.hh"
#include "int_sequence.hh"
#include "rfs_tensor.hh"
#include "t_container.hh"

#include <vector>

/* Here we define the new tensor for representing [G]_S^γ₁…γₗ. It allows a
   construction from container of folded general symmetry tensors ‘cont’, and
   set of indices ‘ts’. Also we have to supply dimensions of resulting tensor
   B, and dimensions of tensor h. */

class USubTensor : public URTensor
{
public:
  USubTensor(const TensorDimens& bdims, const TensorDimens& hdims, const FGSContainer& cont,
             const std::vector<IntSequence>& lst);
  void addKronColumn(int i, const std::vector<const FGSTensor*>& ts, const IntSequence& pindex);
};

#endif
