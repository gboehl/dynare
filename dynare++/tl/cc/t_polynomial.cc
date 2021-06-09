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

#include "t_polynomial.hh"
#include "kron_prod.hh"

// PowerProvider::getNext() unfolded code
/* This method constructs unfolded ‘ut’ of higher dimension, deleting
   the previous. */

const URSingleTensor &
PowerProvider::getNext(dummy<URSingleTensor>)
{
  if (ut)
    {
      auto ut_new = std::make_unique<URSingleTensor>(nv, ut->dimen()+1);
      KronProd::kronMult(ConstVector(origv), ConstVector(ut->getData()), ut_new->getData());
      ut = std::move(ut_new);
    }
  else
    {
      ut = std::make_unique<URSingleTensor>(nv, 1);
      ut->getData() = origv;
    }
  return *ut;
}

// PowerProvider::getNext() folded code
/* This method just constructs next unfolded ‘ut’ and creates folded ‘ft’. */

const FRSingleTensor &
PowerProvider::getNext(dummy<FRSingleTensor>)
{
  getNext<URSingleTensor>();
  ft = std::make_unique<FRSingleTensor>(*ut);
  return *ft;
}

UTensorPolynomial::UTensorPolynomial(const FTensorPolynomial &fp)
  : TensorPolynomial<UFSTensor, UGSTensor, URSingleTensor>(fp.nrows(), fp.nvars())
{
  for (const auto &it : fp)
    insert(std::make_unique<UFSTensor>(*(it.second)));
}

FTensorPolynomial::FTensorPolynomial(const UTensorPolynomial &up)
  : TensorPolynomial<FFSTensor, FGSTensor, FRSingleTensor>(up.nrows(), up.nvars())
{
  for (const auto &it : up)
    insert(std::make_unique<FFSTensor>(*(it.second)));
}
