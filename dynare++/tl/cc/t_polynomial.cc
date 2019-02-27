// Copyright 2004, Ondra Kamenik

#include "t_polynomial.hh"
#include "kron_prod.hh"

// |PowerProvider::getNext| unfolded code
/* This method constructs unfolded |ut| of higher dimension, deleting
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

// |PowerProvider::getNext| folded code
/* This method just constructs next unfolded |ut| and creates folded
   |ft|. */

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
