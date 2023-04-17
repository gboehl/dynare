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

#include "ps_tensor.hh"
#include "fs_tensor.hh"
#include "tl_exception.hh"
#include "tl_static.hh"
#include "stack_container.hh"

#include <iostream>

/* Here we decide, what method for filling a slice in slicing constructor to
   use. A few experiments suggest, that if the tensor is more than 8% filled,
   the first method (fillFromSparseOne()) is better. For fill factors less than
   1%, the second can be 3 times quicker. */

UPSTensor::fill_method
UPSTensor::decideFillMethod(const FSSparseTensor &t)
{
  if (t.getFillFactor() > 0.08)
    return fill_method::first;
  else
    return fill_method::second;
}

/* Here we make a slice. We decide what fill method to use and set it. */

UPSTensor::UPSTensor(const FSSparseTensor &t, const IntSequence &ss,
                     const IntSequence &coor, PerTensorDimens ptd)
  : UTensor(indor::along_col, ptd.getNVX(),
            t.nrows(), ptd.calcUnfoldMaxOffset(), ptd.dimen()),
    tdims(std::move(ptd))
{
  TL_RAISE_IF(coor.size() != t.dimen(),
              "Wrong coordinates length of stacks for UPSTensor slicing constructor");
  TL_RAISE_IF(ss.sum() != t.nvar(),
              "Wrong length of stacks for UPSTensor slicing constructor");

  if (decideFillMethod(t) == fill_method::first)
    fillFromSparseOne(t, ss, coor);
  else
    fillFromSparseTwo(t, ss, coor);
}

void
UPSTensor::increment(IntSequence &v) const
{
  TL_RAISE_IF(v.size() != dimen(),
              "Wrong input/output vector size in UPSTensor::increment");

  UTensor::increment(v, tdims.getNVX());
}

void
UPSTensor::decrement(IntSequence &v) const
{
  TL_RAISE_IF(v.size() != dimen(),
              "Wrong input/output vector size in UPSTensor::decrement");

  UTensor::decrement(v, tdims.getNVX());
}

std::unique_ptr<FTensor>
UPSTensor::fold() const
{
  TL_RAISE("Never should come to this place in UPSTensor::fold");
}

int
UPSTensor::getOffset(const IntSequence &v) const
{
  TL_RAISE_IF(v.size() != dimen(),
              "Wrong input vector size in UPSTensor::getOffset");

  return UTensor::getOffset(v, tdims.getNVX());
}

void
UPSTensor::addTo(FGSTensor &out) const
{
  TL_RAISE_IF(out.getDims() != tdims,
              "Tensors have incompatible dimens in UPSTensor::addTo");
  for (index in = out.begin(); in != out.end(); ++in)
    {
      IntSequence vtmp(dimen());
      tdims.getPer().apply(in.getCoor(), vtmp);
      index tin(*this, vtmp);
      out.addColumn(*this, *tin, *in);
    }
}

/* In here, we have to add this permuted symmetry unfolded tensor to an
   unfolded not permuted tensor. One easy way would be to go through the
   target tensor, permute each index, and add the column.

   However, it may happen, that the permutation has some non-empty
   identity tail. In this case, we can add not only individual columns,
   but much bigger data chunks, which is usually more
   efficient. Therefore, the code is quite dirty, because we have not an
   iterator, which iterates over tensor at some higher levels. So we
   simulate it by the following code.

   First we set ‘cols’ to the length of the data chunk and ‘off’ to its
   dimension. Then we need a front part of ‘nvmax’ of ‘out’, which is
   ‘nvmax_part’. Our iterator here is an integer sequence ‘outrun’ with
   full length, and ‘outrun_part’ its front part. The ‘outrun’ is
   initialized to zeros. In each step we need to increment ‘outrun’
   ‘cols’-times, this is done by incrementing its prefix ‘outrun_part’.

   So we loop over all ‘cols’-wide partitions of ‘out’, permute ‘outrun’
   to obtain ‘perrun’ to obtain column of this matrix. (note that the
   trailing part of ‘perrun’ is the same as of ‘outrun’. Then we
   construct submatrices, add them, and increment ‘outrun’. */

void
UPSTensor::addTo(UGSTensor &out) const
{
  TL_RAISE_IF(out.getDims() != tdims,
              "Tensors have incompatible dimens in UPSTensor::addTo");
  int cols = tailIdentitySize();
  int off = tdims.tailIdentity();
  IntSequence outrun(out.dimen(), 0);
  IntSequence outrun_part(outrun, 0, out.dimen()-off);
  IntSequence nvmax_part(out.getDims().getNVX(), 0, out.dimen()-off);
  for (int out_col = 0; out_col < out.ncols(); out_col += cols)
    {
      // permute ‘outrun’
      IntSequence perrun(out.dimen());
      tdims.getPer().apply(outrun, perrun);
      index from(*this, perrun);
      // construct submatrices
      ConstTwoDMatrix subfrom(*this, *from, cols);
      TwoDMatrix subout(out, out_col, cols);
      // add
      subout.add(1, subfrom);
      // increment ‘outrun’ by cols
      UTensor::increment(outrun_part, nvmax_part);
    }
}

/* This returns a product of all items in ‘nvmax’ which make up the
   trailing identity part. */

int
UPSTensor::tailIdentitySize() const
{
  return tdims.getNVX().mult(dimen()-tdims.tailIdentity(), dimen());
}

/* This fill method is pretty dumb. We go through all columns in ‘this’
   tensor, translate coordinates to sparse tensor, sort them and find an
   item in the sparse tensor. There are many not successful lookups for
   really sparse tensor, that is why the second method works better for
   really sparse tensors. */

void
UPSTensor::fillFromSparseOne(const FSSparseTensor &t, const IntSequence &ss,
                             const IntSequence &coor)
{
  IntSequence cumtmp(ss.size());
  cumtmp[0] = 0;
  for (int i = 1; i < ss.size(); i++)
    cumtmp[i] = cumtmp[i-1] + ss[i-1];
  IntSequence cum(coor.size());
  for (int i = 0; i < coor.size(); i++)
    cum[i] = cumtmp[coor[i]];

  zeros();
  for (Tensor::index run = begin(); run != end(); ++run)
    {
      IntSequence c(run.getCoor());
      c.add(1, cum);
      c.sort();
      auto sl = t.getMap().lower_bound(c);
      if (sl != t.getMap().end())
        {
          auto su = t.getMap().upper_bound(c);
          for (auto srun = sl; srun != su; ++srun)
            get(srun->second.first, *run) = srun->second.second;
        }
    }
}

/* This is the second way of filling the slice. For instance, let the slice
   correspond to partitions “abac”. In here we first calculate lower and upper
   bounds for index of the sparse tensor for the slice. These are ‘lb_srt’ and
   ‘ub_srt’ respectively. They corresponds to ordering “aabc”. Then we go
   through that interval, and select items which are really between the bounds.
   Then we take the index, subtract the lower bound to get it to coordinates of
   the slice. We get something like (i_a,j_a,k_b,l_c). Then we apply the
   inverse permutation as of the sorting form abac ↦ aabc to get index
   (i_a,k_b,j_a,l_c). Recall that the slice is unfolded, so we have to apply
   all permutations preserving the stack coordinates “abac”. In our case we get
   list of indices (i_a,k_b,j_a,l_c) and (j_a,k_b,i_a,l_c). For all we copy the
   item of the sparse tensor to the appropriate column. */

void
UPSTensor::fillFromSparseTwo(const FSSparseTensor &t, const IntSequence &ss,
                             const IntSequence &coor)
{
  IntSequence coor_srt(coor);
  coor_srt.sort();
  IntSequence cum(ss.size());
  cum[0] = 0;
  for (int i = 1; i < ss.size(); i++)
    cum[i] = cum[i-1] + ss[i-1];
  IntSequence lb_srt(coor.size());
  IntSequence ub_srt(coor.size());
  for (int i = 0; i < coor.size(); i++)
    {
      lb_srt[i] = cum[coor_srt[i]];
      ub_srt[i] = cum[coor_srt[i]] + ss[coor_srt[i]] - 1;
    }

  const PermutationSet &pset = TLStatic::getPerm(coor.size());
  std::vector<Permutation> pp = pset.getPreserving(coor);

  Permutation unsort(coor);
  zeros();
  auto lbi = t.getMap().lower_bound(lb_srt);
  auto ubi = t.getMap().upper_bound(ub_srt);
  for (auto run = lbi; run != ubi; ++run)
    {
      if (lb_srt.lessEq(run->first) && run->first.lessEq(ub_srt))
        {
          IntSequence c(run->first);
          c.add(-1, lb_srt);
          unsort.apply(c);
          for (auto &i : pp)
            {
              IntSequence cp(coor.size());
              i.apply(c, cp);
              Tensor::index ind(*this, cp);
              TL_RAISE_IF(*ind < 0 || *ind >= ncols(),
                          "Internal error in slicing constructor of UPSTensor");
              get(run->second.first, *ind) = run->second.second;
            }
        }
    }
}

/* Here we calculate the maximum offsets in each folded dimension
   (dimension sizes, hence ‘ds’). */

void
PerTensorDimens2::setDimensionSizes()
{
  const IntSequence &nvs = getNVS();
  for (int i = 0; i < numSyms(); i++)
    {
      TensorDimens td(syms[i], nvs);
      ds[i] = td.calcFoldMaxOffset();
    }
}

/* If there are two folded dimensions, the offset in such a dimension is offset
   of the second plus offset of the first times the maximum offset of the
   second. If there are n+1 dimensions, the offset is a sum of offsets of the
   last dimension plus the offset in the first n dimensions multiplied by the
   maximum offset of the last dimension. This is exactly what the following
   code does. */

int
PerTensorDimens2::calcOffset(const IntSequence &coor) const
{
  TL_RAISE_IF(coor.size() != dimen(),
              "Wrong length of coordinates in PerTensorDimens2::calcOffset");
  IntSequence cc(coor);
  int ret = 0;
  int off = 0;
  for (int i = 0; i < numSyms(); i++)
    {
      TensorDimens td(syms[i], getNVS());
      IntSequence c(cc, off, off+syms[i].dimen());
      int a = td.calcFoldOffset(c);
      ret = ret*ds[i] + a;
      off += syms[i].dimen();
    }
  return ret;
}

void
PerTensorDimens2::print() const
{
  std::cout << "nvmax: ";
  nvmax.print();
  std::cout << "per:   ";
  per.print();
  std::cout << "syms:  ";
  syms.print();
  std::cout << "dims:  ";
  ds.print();
}

/* Here we increment the given integer sequence. It corresponds to
   UTensor::increment() of the whole sequence, and then partial monotonizing of
   the subsequences with respect to the symmetries of each dimension. */

void
FPSTensor::increment(IntSequence &v) const
{
  TL_RAISE_IF(v.size() != dimen(),
              "Wrong length of coordinates in FPSTensor::increment");
  UTensor::increment(v, tdims.getNVX());
  int off = 0;
  for (int i = 0; i < tdims.numSyms(); i++)
    {
      IntSequence c(v, off, off+tdims.getSym(i).dimen());
      c.pmonotone(tdims.getSym(i));
      off += tdims.getSym(i).dimen();
    }
}

void
FPSTensor::decrement(IntSequence &v) const
{
  TL_RAISE("FPSTensor::decrement not implemented");
}

std::unique_ptr<UTensor>
FPSTensor::unfold() const
{
  TL_RAISE("Unfolding of FPSTensor not implemented");
}

/* We only call calcOffset() on the PerTensorDimens2. */

int
FPSTensor::getOffset(const IntSequence &v) const
{
  return tdims.calcOffset(v);
}

/* Here we add the tensor to ‘out’. We go through all columns of the
   ‘out’, apply the permutation to get index in the tensor, and add the
   column. Note that if the permutation is identity, then the dimensions
   of the tensors might not be the same (since this tensor is partially
   folded). */

void
FPSTensor::addTo(FGSTensor &out) const
{
  for (index tar = out.begin(); tar != out.end(); ++tar)
    {
      IntSequence coor(dimen());
      tdims.getPer().apply(tar.getCoor(), coor);
      index src(*this, coor);
      out.addColumn(*this, *src, *tar);
    }
}

/* Here is the constructor which multiplies the Kronecker product with
   the general symmetry sparse tensor GSSparseTensor. The main idea is
   to go through items in the sparse tensor (each item selects rows in
   the matrices form the Kornecker product), then to Kronecker multiply
   the rows and multiply with the item, and to add the resulting row to
   the appropriate row of the resulting FPSTensor.

   The realization of this idea is a bit more complicated since we have
   to go through all items, and each item must be added as many times as
   it has its symmetric elements. Moreover, the permutations shuffle
   order of rows in their Kronecker product.

   So, we through all unfolded indices in a tensor with the same
   dimensions as the GSSparseTensor (sparse slice). For each such index
   we calculate its folded version (corresponds to ordering of
   subsequences within symmetries), we test if there is an item in the
   sparse slice with such coordinates, and if there is, we construct the
   Kronecker product of the rows, and go through all of items with the
   coordinates, and add to appropriate rows of ‘this’ tensor. */

FPSTensor::FPSTensor(const TensorDimens &td, const Equivalence &e, const Permutation &p,
                     const GSSparseTensor &a, const KronProdAll &kp)
  : FTensor(indor::along_col, PerTensorDimens(td, Permutation(e, p)).getNVX(),
            a.nrows(), kp.ncols(), td.dimen()),
    tdims(td, e, p)
{
  zeros();

  UGSTensor dummy(0, a.getDims());
  for (Tensor::index run = dummy.begin(); run != dummy.end(); ++run)
    {
      Tensor::index fold_ind = dummy.getFirstIndexOf(run);
      const IntSequence &c = fold_ind.getCoor();
      auto sl = a.getMap().lower_bound(c);
      if (sl != a.getMap().end())
        {
          auto row_prod = kp.multRows(run.getCoor());
          auto su = a.getMap().upper_bound(c);
          for (auto srun = sl; srun != su; ++srun)
            {
              Vector out_row{getRow(srun->second.first)};
              out_row.add(srun->second.second, *row_prod);
            }
        }
    }
}
