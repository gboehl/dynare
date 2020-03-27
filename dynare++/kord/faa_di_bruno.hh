/*
 * Copyright © 2005 Ondra Kamenik
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
 * along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
 */

// Faà Di Bruno evaluator

/* This defines a class which implements the Faà Di Bruno Formula

                                         ₗ
    [B_sᵏ]_α₁…αₗ = [f_zˡ]_β₁…βₗ    ∑     ∏  [z_{s^|cₘ|}]_cₘ(α)^βₘ
                                c∈ℳₗ,ₖ ᵐ⁼¹

   where sᵏ is a general symmetry of dimension k and z is a stack of
   functions. */

#ifndef FAA_DI_BRUNO_H
#define FAA_DI_BRUNO_H

#include "journal.hh"
#include "stack_container.hh"
#include "t_container.hh"
#include "sparse_tensor.hh"
#include "gs_tensor.hh"

#include <tuple>

class FaaDiBruno
{
  Journal &journal;
public:
  FaaDiBruno(Journal &jr)
    : journal(jr)
  {
  }
  void calculate(const StackContainer<FGSTensor> &cont, const TensorContainer<FSSparseTensor> &f,
                 FGSTensor &out);
  void calculate(const FoldedStackContainer &cont, const FGSContainer &g,
                 FGSTensor &out);
  void calculate(const StackContainer<UGSTensor> &cont, const TensorContainer<FSSparseTensor> &f,
                 UGSTensor &out);
  void calculate(const UnfoldedStackContainer &cont, const UGSContainer &g,
                 UGSTensor &out);
protected:
  std::tuple<int, int, int> estimRefinement(const TensorDimens &tdims, int nr, int l);

  // See FaaDiBruno::calculate() folded sparse code for why we have magic_mult
  constexpr static double magic_mult = 1.5;
};

#endif
