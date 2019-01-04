// Copyright 2005, Ondra Kamenik

// Faa Di Bruno evaluator

/* This defines a class which implements Faa Di Bruno Formula
   $$\left[B_{s^k}\right]_{\alpha_1\ldots\alpha_l}=\left[f_{z^l}\right]_{\beta_1\ldots\beta_l}
   \sum_{c\in M_{l,k}}\prod_{m=1}^l\left[z_{s^k(c_m)}\right]^{\beta_m}_{c_m(\alpha)}$$
   where $s^k$ is a general symmetry of dimension $k$ and $z$ is a stack of
   functions. */

#ifndef FAA_DI_BRUNO_H
#define FAA_DI_BRUNO_H

#include "journal.hh"
#include "stack_container.h"
#include "t_container.h"
#include "sparse_tensor.h"
#include "gs_tensor.h"

/* Nothing special here. See |@<|FaaDiBruno::calculate| folded sparse
   code@>| for reason of having |magic_mult|. */

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
  int estimRefinment(const TensorDimens &tdims, int nr, int l, int &avmem_mb, int &tmpmem_mb);
  static double magic_mult;
};

#endif
