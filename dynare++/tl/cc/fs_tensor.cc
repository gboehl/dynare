// Copyright 2004, Ondra Kamenik

#include "fs_tensor.hh"
#include "gs_tensor.hh"
#include "sparse_tensor.hh"
#include "rfs_tensor.hh"
#include "tl_exception.hh"
#include "pascal_triangle.hh"

/* This constructs a fully symmetric tensor as given by the contraction:
   $$\left[g_{y^n}\right]_{\alpha_1\ldots\alpha_n}=
   \left[t_{y^{n+1}}\right]_{\alpha_1\ldots\alpha_n\beta}[x]^\beta$$

   We go through all columns of output tensor $[g]$ and for each column
   we cycle through all variables, insert a variable to the column
   coordinates obtaining a column of tensor $[t]$. the column is multiplied
   by an appropriate item of |x| and added to the column of $[g]$ tensor. */

FFSTensor::FFSTensor(const FFSTensor &t, const ConstVector &x)
  : FTensor(indor::along_col, IntSequence(t.dimen()-1, t.nvar()),
            t.nrows(), calcMaxOffset(t.nvar(), t.dimen()-1), t.dimen()-1),
    nv(t.nvar())
{
  TL_RAISE_IF(t.dimen() < 1,
              "Wrong dimension for tensor contraction of FFSTensor");
  TL_RAISE_IF(t.nvar() != x.length(),
              "Wrong number of variables for tensor contraction of FFSTensor");

  zeros();

  for (Tensor::index to = begin(); to != end(); ++to)
    for (int i = 0; i < nvar(); i++)
      {
        IntSequence from_ind(to.getCoor().insert(i));
        Tensor::index from(t, from_ind);
        addColumn(x[i], t, *from, *to);
      }
}

/* This returns number of indices for folded tensor with full
   symmetry. Let $n$ be a number of variables |nvar| and $d$ the
   dimension |dim|. Then the number of indices is $\pmatrix{n+d-1\cr d}$. */

int
FFSTensor::calcMaxOffset(int nvar, int d)
{
  if (nvar == 0 && d == 0)
    return 1;
  if (nvar == 0 && d > 0)
    return 0;
  return PascalTriangle::noverk(nvar + d - 1, d);
}

/* The conversion from sparse tensor is clear. We go through all the
   tensor and write to the dense what is found. */
FFSTensor::FFSTensor(const FSSparseTensor &t)
  : FTensor(indor::along_col, IntSequence(t.dimen(), t.nvar()),
            t.nrows(), calcMaxOffset(t.nvar(), t.dimen()), t.dimen()),
    nv(t.nvar())
{
  zeros();
  for (const auto & it : t.getMap())
    {
      index ind(*this, it.first);
      get(it.second.first, *ind) = it.second.second;
    }
}

/* The conversion from unfolded copies only columns of respective
   coordinates. So we go through all the columns in the folded tensor
   (this), make an index of the unfolded vector from coordinates, and
   copy the column. */

FFSTensor::FFSTensor(const UFSTensor &ut)
  : FTensor(indor::along_col, IntSequence(ut.dimen(), ut.nvar()),
            ut.nrows(), calcMaxOffset(ut.nvar(), ut.dimen()), ut.dimen()),
    nv(ut.nvar())
{
  for (index in = begin(); in != end(); ++in)
    {
      index src(ut, in.getCoor());
      copyColumn(ut, *src, *in);
    }
}

/* Here just make a new instance and return the reference. */
std::unique_ptr<UTensor>
FFSTensor::unfold() const
{
  return std::make_unique<UFSTensor>(*this);
}

/* Incrementing is easy. We have to increment by calling static method
   |UTensor::increment| first. In this way, we have coordinates of
   unfolded tensor. Then we have to skip to the closest folded index
   which corresponds to monotonizeing the integer sequence. */

void
FFSTensor::increment(IntSequence &v) const
{
  TL_RAISE_IF(v.size() != dimen(),
              "Wrong input/output vector size in FFSTensor::increment");

  UTensor::increment(v, nv);
  v.monotone();
}

/* Decrement calls static |FTensor::decrement|. */

void
FFSTensor::decrement(IntSequence &v) const
{
  TL_RAISE_IF(v.size() != dimen(),
              "Wrong input/output vector size in FFSTensor::decrement");

  FTensor::decrement(v, nv);
}

int
FFSTensor::getOffset(const IntSequence &v) const
{
  TL_RAISE_IF(v.size() != dimen(),
              "Wrong input vector size in FFSTensor::getOffset");

  return FTensor::getOffset(v, nv);
}

/* Here we add a general symmetry tensor to the (part of) full symmetry
   tensor provided that the unique variable of the full symmetry tensor
   is a stack of variables from the general symmetry tensor.

   We check for the dimensions and number of variables. Then we calculate
   a shift of coordinates when going from the general symmetry tensor to
   full symmetry (it corresponds to shift of coordinates induces by
   stacking the variables). Then we add the appropriate columns by going
   through the columns in general symmetry, adding the shift and sorting. */

void
FFSTensor::addSubTensor(const FGSTensor &t)
{
  TL_RAISE_IF(dimen() != t.getDims().dimen(),
              "Wrong dimensions for FFSTensor::addSubTensor");
  TL_RAISE_IF(nvar() != t.getDims().getNVS().sum(),
              "Wrong nvs for FFSTensor::addSubTensor");

  // set shift for |addSubTensor|
  /* Code shared with UFSTensor::addSubTensor() */
  IntSequence shift_pre(t.getSym().num(), 0);
  for (int i = 1; i < t.getSym().num(); i++)
    shift_pre[i] = shift_pre[i-1]+t.getDims().getNVS()[i-1];
  IntSequence shift(shift_pre.unfold(t.getSym()));

  for (Tensor::index ind = t.begin(); ind != t.end(); ++ind)
    {
      IntSequence c(ind.getCoor());
      c.add(1, shift);
      c.sort();
      Tensor::index tar(*this, c);
      addColumn(t, *ind, *tar);
    }
}

// |UFSTensor| contraction constructor
/* This is a bit more straightforward than |@<|FFSTensor| contraction constructor@>|.
   We do not add column by column but we do it by submatrices due to
   regularity of the unfolded tensor. */

UFSTensor::UFSTensor(const UFSTensor &t, const ConstVector &x)
  : UTensor(indor::along_col, IntSequence(t.dimen()-1, t.nvar()),
            t.nrows(), calcMaxOffset(t.nvar(), t.dimen()-1), t.dimen()-1),
    nv(t.nvar())
{
  TL_RAISE_IF(t.dimen() < 1,
              "Wrong dimension for tensor contraction of UFSTensor");
  TL_RAISE_IF(t.nvar() != x.length(),
              "Wrong number of variables for tensor contraction of UFSTensor");

  zeros();

  for (int i = 0; i < ncols(); i++)
    {
      ConstTwoDMatrix tpart(t, i *nvar(), nvar());
      Vector outcol{getCol(i)};
      tpart.multaVec(outcol, x);
    }
}

/* Here we convert folded full symmetry tensor to unfolded. We copy all
   columns of folded tensor, and then call |unfoldData()|. */

UFSTensor::UFSTensor(const FFSTensor &ft)
  : UTensor(indor::along_col, IntSequence(ft.dimen(), ft.nvar()),
            ft.nrows(), calcMaxOffset(ft.nvar(), ft.dimen()), ft.dimen()),
    nv(ft.nvar())
{
  for (index src = ft.begin(); src != ft.end(); ++src)
    {
      index in(*this, src.getCoor());
      copyColumn(ft, *src, *in);
    }
  unfoldData();
}

std::unique_ptr<FTensor>
UFSTensor::fold() const
{
  return std::make_unique<FFSTensor>(*this);
}

// |UFSTensor| increment and decrement
/* Here we just call |UTensor| respective static methods. */
void
UFSTensor::increment(IntSequence &v) const
{
  TL_RAISE_IF(v.size() != dimen(),
              "Wrong input/output vector size in UFSTensor::increment");

  UTensor::increment(v, nv);
}

void
UFSTensor::decrement(IntSequence &v) const
{
  TL_RAISE_IF(v.size() != dimen(),
              "Wrong input/output vector size in UFSTensor::decrement");

  UTensor::decrement(v, nv);
}

int
UFSTensor::getOffset(const IntSequence &v) const
{
  TL_RAISE_IF(v.size() != dimen(),
              "Wrong input vector size in UFSTensor::getOffset");

  return UTensor::getOffset(v, nv);
}

/* This is very similar to |@<|FFSTensor::addSubTensor| code@>|. The
   only difference is the addition. We go through all columns in the full
   symmetry tensor and cancel the shift. If the coordinates after the
   cancellation are positive, we find the column in the general symmetry
   tensor, and add it. */

void
UFSTensor::addSubTensor(const UGSTensor &t)
{
  TL_RAISE_IF(dimen() != t.getDims().dimen(),
              "Wrong dimensions for UFSTensor::addSubTensor");
  TL_RAISE_IF(nvar() != t.getDims().getNVS().sum(),
              "Wrong nvs for UFSTensor::addSubTensor");

  // set shift for |addSubTensor|
  /* Code shared with FFSTensor::addSubTensor() */
  IntSequence shift_pre(t.getSym().num(), 0);
  for (int i = 1; i < t.getSym().num(); i++)
    shift_pre[i] = shift_pre[i-1]+t.getDims().getNVS()[i-1];
  IntSequence shift(shift_pre.unfold(t.getSym()));

  for (Tensor::index tar = begin(); tar != end(); ++tar)
    {
      IntSequence c(tar.getCoor());
      c.sort();
      c.add(-1, shift);
      if (c.isPositive() && c.less(t.getDims().getNVX()))
        {
          Tensor::index from(t, c);
          addColumn(t, *from, *tar);
        }
    }
}

/* Here we go through all columns, find a column of folded index, and
   then copy the column data. Finding the index is done by sorting the
   integer sequence. */

void
UFSTensor::unfoldData()
{
  for (index in = begin(); in != end(); ++in)
    {
      IntSequence v(in.getCoor());
      v.sort();
      copyColumn(*index(*this, v), *in);
    }
}
