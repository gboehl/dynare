// Copyright 2004, Ondra Kamenik

// Full symmetry tensor.

/* Here we define folded and unfolded tensors for full symmetry. All
   tensors from here are identifying the multidimensional index with
   columns. */

#ifndef FS_TENSOR_H
#define FS_TENSOR_H

#include "tensor.hh"
#include "symmetry.hh"

class FGSTensor;
class UGSTensor;
class FRSingleTensor;
class FSSparseTensor;

/* Folded tensor with full symmetry maintains only information about
   number of symmetrical variables |nv|. Further, we implement what is
   left from the super class |FTensor|.

   We implement |getOffset| which should be used with care since
   its complexity.

   We implement a method adding a given general symmetry tensor to the
   full symmetry tensor supposing the variables of the general symmetry
   tensor are stacked giving only one variable of the full symmetry
   tensor. For instance, if $x=[y^T, u^T]^T$, then we can add tensor
   $\left[g_{y^2u}\right]$ to tensor $g_{x^3}$. This is done in method
   |addSubTensor|. Consult |@<|FGSTensor| class declaration@>| to know
   what is general symmetry tensor. */

class UFSTensor;
class FFSTensor : public FTensor
{
  int nv;
public:
  /* Here are the constructors. The second constructor constructs a
     tensor by one-dimensional contraction from the higher dimensional
     tensor |t|. This is, it constructs a tensor
     $$\left[g_{y^n}\right]_{\alpha_1\ldots\alpha_n}=
     \left[t_{y^{n+1}}\right]_{\alpha_1\ldots\alpha_n\beta}[x]^\beta$$
     See implementation |@<|FFSTensor| contraction constructor@>| for details.

     The next constructor converts from sparse tensor (which is fully
     symmetric and folded by nature).

     The fourth constructs object from unfolded fully symmetric.

     The fifth constructs a subtensor of selected rows. */

  FFSTensor(int r, int nvar, int d)
    : FTensor(indor::along_col, IntSequence(d, nvar),
              r, calcMaxOffset(nvar, d), d), nv(nvar)
  {
  }
  FFSTensor(const FFSTensor &t, const ConstVector &x);
  FFSTensor(const FSSparseTensor &t);
  FFSTensor(const FFSTensor &ft)
     
  = default;
  FFSTensor(const UFSTensor &ut);
  FFSTensor(int first_row, int num, FFSTensor &t)
    : FTensor(first_row, num, t), nv(t.nv)
  {
  }

  void increment(IntSequence &v) const override;
  void decrement(IntSequence &v) const override;
  UTensor&unfold() const override;
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
   |FFSTensor|, but the method |unfoldData|. It takes columns which also
   exist in folded version and copies them to all their symmetrical
   locations. This is useful when constructing unfolded tensor from
   folded one. */

class UFSTensor : public UTensor
{
  int nv;
public:
  UFSTensor(int r, int nvar, int d)
    : UTensor(indor::along_col, IntSequence(d, nvar),
              r, calcMaxOffset(nvar, d), d), nv(nvar)
  {
  }
  UFSTensor(const UFSTensor &t, const ConstVector &x);
  UFSTensor(const UFSTensor &ut)
     
  = default;
  UFSTensor(const FFSTensor &ft);
  UFSTensor(int first_row, int num, UFSTensor &t)
    : UTensor(first_row, num, t), nv(t.nv)
  {
  }

  void increment(IntSequence &v) const override;
  void decrement(IntSequence &v) const override;
  FTensor&fold() const override;
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
