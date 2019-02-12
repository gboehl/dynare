// Copyright 2004, Ondra Kamenik

#include "rfs_tensor.hh"
#include "kron_prod.hh"
#include "tl_exception.hh"

// |FRTensor| conversion from unfolded
/* The conversion from unfolded to folded sums up all data from
   unfolded corresponding to one folded index. So we go through all the
   rows in the unfolded tensor |ut|, make an index of the folded tensor
   by sorting the coordinates, and add the row. */
FRTensor::FRTensor(const URTensor &ut)
  : FTensor(indor::along_row, IntSequence(ut.dimen(), ut.nvar()),
            FFSTensor::calcMaxOffset(ut.nvar(), ut.dimen()), ut.ncols(),
            ut.dimen()),
    nv(ut.nvar())
{
  zeros();
  for (index in = ut.begin(); in != ut.end(); ++in)
    {
      IntSequence vtmp(in.getCoor());
      vtmp.sort();
      index tar(*this, vtmp);
      addRow(ut, *in, *tar);
    }
}

/* Here just make a new instance and return the reference. */

UTensor &
FRTensor::unfold() const
{
  return *(new URTensor(*this));
}

/* Incrementing is easy. The same as for |FFSTensor|. */

void
FRTensor::increment(IntSequence &v) const
{
  TL_RAISE_IF(v.size() != dimen(),
              "Wrong input/output vector size in FRTensor::increment");

  UTensor::increment(v, nv);
  v.monotone();
}

/* Decrement calls static |FTensor::decrement|. */

void
FRTensor::decrement(IntSequence &v) const
{
  TL_RAISE_IF(v.size() != dimen(),
              "Wrong input/output vector size in FRTensor::decrement");

  FTensor::decrement(v, nv);
}

// |URTensor| conversion from folded
/* Here we convert folded full symmetry tensor to unfolded. We copy all
   columns of folded tensor to unfolded and leave other columns
   (duplicates) zero. In this way, if the unfolded tensor is folded back,
   we should get the same data. */
URTensor::URTensor(const FRTensor &ft)
  : UTensor(indor::along_row, IntSequence(ft.dimen(), ft.nvar()),
            UFSTensor::calcMaxOffset(ft.nvar(), ft.dimen()), ft.ncols(),
            ft.dimen()),
    nv(ft.nvar())
{
  zeros();
  for (index src = ft.begin(); src != ft.end(); ++src)
    {
      index in(*this, src.getCoor());
      copyRow(ft, *src, *in);
    }
}

/* Here we just return a reference to new instance of folded tensor. */

FTensor &
URTensor::fold() const
{
  return *(new FRTensor(*this));
}

/* Here we just call |UTensor| respective static methods. */

void
URTensor::increment(IntSequence &v) const
{
  TL_RAISE_IF(v.size() != dimen(),
              "Wrong input/output vector size in URTensor::increment");

  UTensor::increment(v, nv);
}

void
URTensor::decrement(IntSequence &v) const
{
  TL_RAISE_IF(v.size() != dimen(),
              "Wrong input/output vector size in URTensor::decrement");

  UTensor::decrement(v, nv);
}

int
URTensor::getOffset(const IntSequence &v) const
{
  TL_RAISE_IF(v.size() != dimen(),
              "Wrong input vector size in URTensor::getOffset");

  return UTensor::getOffset(v, nv);
}

/* Here we construct $v_1\otimes v_2\otimes\ldots\otimes v_n$, where
   $v_1,v_2,\ldots,v_n$ are stored in |vector<ConstVector>|. */

URSingleTensor::URSingleTensor(const std::vector<ConstVector> &cols)
  : URTensor(1, cols[0].length(), cols.size())
{
  if (dimen() == 1)
    {
      getData() = cols[0];
      return;
    }

  auto *last = new Vector(cols[cols.size()-1]);
  for (int i = cols.size()-2; i > 0; i--)
    {
      auto *newlast = new Vector(Tensor::power(nvar(), cols.size()-i));
      KronProd::kronMult(cols[i], ConstVector(*last), *newlast);
      delete last;
      last = newlast;
    }
  KronProd::kronMult(cols[0], ConstVector(*last), getData());
  delete last;
}

/* Here we construct $v\otimes\ldots\otimes v$, where the number of $v$
   copies is |d|. */

URSingleTensor::URSingleTensor(const ConstVector &v, int d)
  : URTensor(1, v.length(), d)
{
  if (d == 1)
    {
      getData() = v;
      return;
    }

  auto *last = new Vector(v);
  for (int i = d-2; i > 0; i--)
    {
      auto *newlast = new Vector(last->length()*v.length());
      KronProd::kronMult(v, ConstVector(*last), *newlast);
      delete last;
      last = newlast;
    }
  KronProd::kronMult(v, ConstVector(*last), getData());
  delete last;
}

/* Here we construct |FRSingleTensor| from |URSingleTensor| and return
   its reference. */

FTensor &
URSingleTensor::fold() const
{
  return *(new FRSingleTensor(*this));
}

// |FRSingleTensor| conversion from unfolded
/* The conversion from unfolded |URSingleTensor| to folded
   |FRSingleTensor| is completely the same as conversion from |URTensor|
   to |FRTensor|, only we do not copy rows but elements. */
FRSingleTensor::FRSingleTensor(const URSingleTensor &ut)
  : FRTensor(1, ut.nvar(), ut.dimen())
{
  zeros();
  for (index in = ut.begin(); in != ut.end(); ++in)
    {
      IntSequence vtmp(in.getCoor());
      vtmp.sort();
      index tar(*this, vtmp);
      get(*tar, 0) += ut.get(*in, 0);
    }
}
