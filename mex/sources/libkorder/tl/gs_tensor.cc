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

#include "gs_tensor.hh"
#include "sparse_tensor.hh"
#include "tl_exception.hh"
#include "kron_prod.hh"
#include "pascal_triangle.hh"

/* Constructor used for slicing fully symmetric tensor. It constructs the
   dimensions from the partitioning of variables of fully symmetric tensor. Let
   the partitioning be, for instance, (a,b,c,d), where (n_a,n_b,n_c,n_d)
   are lengths of the partitions. Let one want to get a slice only of the part
   of the fully symmetric tensor corresponding to indices of the form b²d³.
   This corresponds to the symmetry a⁰b²c⁰d³. So, the dimension of the
   slice would be also (n_a,n_b,n_c,n_d) for number of variables and
   (0,2,0,3) for the symmetry. So we provide the constructor which takes
   sizes of partitions (n_a,n_b,n_c,n_d) as IntSequence, and indices of
   picked partitions, in our case (1,1,3,3,3), as IntSequence. */

TensorDimens::TensorDimens(const IntSequence &ss, const IntSequence &coor)
  : nvs(ss),
    sym(ss.size()),
    nvmax(coor.size(), 0)
{
  TL_RAISE_IF(!coor.isSorted(),
              "Coordinates not sorted in TensorDimens slicing constructor");
  TL_RAISE_IF(coor[0] < 0 || coor[coor.size()-1] >= ss.size(),
              "A coordinate out of stack range in TensorDimens slicing constructor");

  for (int i = 0; i < coor.size(); i++)
    {
      sym[coor[i]]++;
      nvmax[i] = ss[coor[i]];
    }
}

/* Number of unfold offsets is a product of all members of nvmax. */
int
TensorDimens::calcUnfoldMaxOffset() const
{
  return nvmax.mult();
}

/* Number of folded offsets is a product of all unfold offsets within
   each equivalence class of the symmetry. */

int
TensorDimens::calcFoldMaxOffset() const
{
  int res = 1;
  for (int i = 0; i < nvs.size(); i++)
    {
      if (nvs[i] == 0 && sym[i] > 0)
        return 0;
      if (sym[i] > 0)
        res *= PascalTriangle::noverk(nvs[i]+sym[i]-1, sym[i]);
    }
  return res;
}

/* Here we implement offset calculation for folded general symmetry
   tensor. The offset of a given sequence is calculated by breaking the
   sequence to subsequences according to the symmetry. The offset is
   orthogonal with respect to the blocks, this means that indexing within
   the blocks is independent. If there are two blocks, for instance, then
   the offset will be an offset within the outer block (the first)
   multiplied with all offsets of the inner block (last) plus an offset
   within the second block.

   Generally, the resulting offset will be
        ₛ      ₛ
   r =  ∑ rᵢ·  ∏  nⱼ
       ⁱ⁼¹   ʲ⁼ⁱ⁺¹
   where s is the number of blocks (getSym().num()), rᵢ is the offset
   within the i-th block, and nⱼ is the number of all offsets in the j-th
   block.

   In the code, we go from the innermost to the outermost, maintaining the
   product in ‘pow’. */

int
TensorDimens::calcFoldOffset(const IntSequence &v) const
{
  TL_RAISE_IF(v.size() != dimen(),
              "Wrong input vector size in TensorDimens::getFoldOffset");

  int res = 0;
  int pow = 1;
  int blstart = v.size();
  for (int ibl = getSym().num()-1; ibl >= 0; ibl--)
    {
      int bldim = getSym()[ibl];
      if (bldim > 0)
        {
          blstart -= bldim;
          int blnvar = getNVX(blstart);
          IntSequence subv(v, blstart, blstart+bldim);
          res += FTensor::getOffset(subv, blnvar)*pow;
          pow *= FFSTensor::calcMaxOffset(blnvar, bldim);
        }
    }
  TL_RAISE_IF(blstart != 0,
              "Error in tracing symmetry in TensorDimens::getFoldOffset");
  return res;
}

/* In order to find the predecessor of index within folded generally
   symmetric tensor, note, that a decrease action in i-th partition of
   symmetric indices can happen only if all indices in all subsequent
   partitions are zero. Then the decrease action of whole the index
   consists of decrease action of the first nonzero partition from the
   right, and setting these trailing zero partitions to their maximum
   indices.

   So we set ‘iblock’ to the number of last partitions. During the
   execution, ‘block_first’, and ‘block_last’ will point to the first
   element of ‘iblock’ and, first element of following block.

   Then we check for all trailing zero partitions, set them to their
   maximums and return ‘iblock’ to point to the first non-zero partition
   (or the first partition). Then for this partition, we decrease the
   index (fully symmetric within that partition). */

void
TensorDimens::decrement(IntSequence &v) const
{
  TL_RAISE_IF(getNVX().size() != v.size(),
              "Wrong size of input/output sequence in TensorDimens::decrement");

  int iblock = getSym().num()-1;
  int block_last = v.size();
  int block_first = block_last-getSym()[iblock];

  // check for zero trailing blocks
  while (iblock > 0 && v[block_last-1] == 0)
    {
      for (int i = block_first; i < block_last; i++)
        v[i] = getNVX(i); // equivalent to nvs[iblock]
      iblock--;
      block_last = block_first;
      block_first -= getSym()[iblock];
    }

  // decrease the non-zero block
  IntSequence vtmp(v, block_first, block_last);
  FTensor::decrement(vtmp, getNVX(block_first));
}

// FGSTensor conversion from UGSTensor
/* Here we go through columns of folded, calculate column of unfolded,
   and copy data. */
FGSTensor::FGSTensor(const UGSTensor &ut)
  : FTensor(indor::along_col, ut.tdims.getNVX(), ut.nrows(),
            ut.tdims.calcFoldMaxOffset(), ut.dimen()),
    tdims(ut.tdims)
{
  for (index ti = begin(); ti != end(); ++ti)
    {
      index ui(ut, ti.getCoor());
      copyColumn(ut, *ui, *ti);
    }
}

// FGSTensor slicing constructor from FSSparseTensor
/* Here is the code of slicing constructor from the sparse tensor. We
   first calculate coordinates of first and last index of the slice
   within the sparse tensor (these are ‘lb’ and ‘ub’), and then we
   iterate through all items between them (in lexicographical ordering of
   sparse tensor), and check whether an item is between the ‘lb’ and ‘ub’
   in Cartesian ordering (this corresponds to belonging to the
   slices). If it belongs, then we subtract the lower bound ‘lb’ to
   obtain coordinates in the ‘this’ tensor and we copy the item. */
FGSTensor::FGSTensor(const FSSparseTensor &t, const IntSequence &ss,
                     const IntSequence &coor, TensorDimens td)
  : FTensor(indor::along_col, td.getNVX(), t.nrows(),
            td.calcFoldMaxOffset(), td.dimen()),
    tdims(std::move(td))
{
  // set ‘lb’ and ‘ub’ to lower and upper bounds of indices
  /* Here we first set ‘s_offsets’ to offsets of partitions whose lengths
     are given by ‘ss’. So ‘s_offsets’ is a cumulative sum of ‘ss’.

     Then we create ‘lb’ to be coordinates of the possibly first index from
     the slice, and ‘ub’ to be coordinates of possibly last index of the
     slice. */
  IntSequence s_offsets(ss.size(), 0);
  for (int i = 1; i < ss.size(); i++)
    s_offsets[i] = s_offsets[i-1] + ss[i-1];

  IntSequence lb(coor.size());
  IntSequence ub(coor.size());
  for (int i = 0; i < coor.size(); i++)
    {
      lb[i] = s_offsets[coor[i]];
      ub[i] = s_offsets[coor[i]] + ss[coor[i]] - 1;
    }

  zeros();
  auto lbi = t.getMap().lower_bound(lb);
  auto ubi = t.getMap().upper_bound(ub);
  for (auto run = lbi; run != ubi; ++run)
    if (lb.lessEq(run->first) && run->first.lessEq(ub))
      {
        IntSequence c(run->first);
        c.add(-1, lb);
        Tensor::index ind(*this, c);
        TL_RAISE_IF(*ind < 0 || *ind >= ncols(),
                    "Internal error in slicing constructor of FGSTensor");
        get(run->second.first, *ind) = run->second.second;
      }
}

// FGSTensor slicing from FFSTensor
/* The code is similar to FGSTensor slicing from FSSparseTensor. */
FGSTensor::FGSTensor(const FFSTensor &t, const IntSequence &ss,
                     const IntSequence &coor, TensorDimens td)
  : FTensor(indor::along_col, td.getNVX(), t.nrows(),
            td.calcFoldMaxOffset(), td.dimen()),
    tdims(std::move(td))
{
  if (ncols() == 0)
    return;

  // set ‘lb’ and ‘ub’ to lower and upper bounds of indices
  /* Same code as in the previous converting constructor */
  IntSequence s_offsets(ss.size(), 0);
  for (int i = 1; i < ss.size(); i++)
    s_offsets[i] = s_offsets[i-1] + ss[i-1];

  IntSequence lb(coor.size());
  IntSequence ub(coor.size());
  for (int i = 0; i < coor.size(); i++)
    {
      lb[i] = s_offsets[coor[i]];
      ub[i] = s_offsets[coor[i]] + ss[coor[i]] - 1;
    }

  zeros();
  Tensor::index lbi(t, lb);
  Tensor::index ubi(t, ub);
  ++ubi;
  for (Tensor::index run = lbi; run != ubi; ++run)
    {
      if (lb.lessEq(run.getCoor()) && run.getCoor().lessEq(ub))
        {
          IntSequence c(run.getCoor());
          c.add(-1, lb);
          Tensor::index ind(*this, c);
          TL_RAISE_IF(*ind < 0 || *ind >= ncols(),
                      "Internal error in slicing constructor of FGSTensor");
          copyColumn(t, *run, *ind);
        }
    }
}

// FGSTensor conversion from GSSparseTensor
FGSTensor::FGSTensor(const GSSparseTensor &t)
  : FTensor(indor::along_col, t.getDims().getNVX(), t.nrows(),
            t.getDims().calcFoldMaxOffset(), t.dimen()), tdims(t.getDims())
{
  zeros();
  for (const auto &it : t.getMap())
    {
      index ind(*this, it.first);
      get(it.second.first, *ind) = it.second.second;
    }
}

/* First we increment as unfolded, then we must monotonize within
   partitions defined by the symmetry. This is done by
   IntSequence::pmonotone(). */

void
FGSTensor::increment(IntSequence &v) const
{
  TL_RAISE_IF(v.size() != dimen(),
              "Wrong input/output vector size in FGSTensor::increment");

  UTensor::increment(v, tdims.getNVX());
  v.pmonotone(tdims.getSym());
}

/* Return unfolded version of the tensor. */
std::unique_ptr<UTensor>
FGSTensor::unfold() const
{
  return std::make_unique<UGSTensor>(*this);
}

/* Here we implement the contraction

    [r_xⁱzᵏ]_α₁…αᵢγ₁…γₖ = [t_xⁱyʲzᵏ]_α₁…αᵢβ₁…βⱼγ₁…γₖ·[c]^β₁…βⱼ

   More generally, xⁱ and zᵏ can represent also general symmetries.

   The operation can be rewritten as a matrix product

    [t_xⁱyʲzᵏ]·(Iₗ⊗c⊗Iᵣ)

   where l is a number of columns in tensor with symmetry on the left
   (i.e. xⁱ), and r is a number of columns in tensor with a symmetry
   on the right (i.e. zᵏ). The code proceeds accordingly. We first
   form two symmetries ‘sym_left’ and ‘sym_right’, then calculate the
   number of columns dleft=l and dright=r, form the Kronecker
   product and multiply and add.

   The input parameter ‘i’ is the order of a variable being contracted
   starting from 0. */

void
FGSTensor::contractAndAdd(int i, FGSTensor &out,
                          const FRSingleTensor &col) const
{
  TL_RAISE_IF(i < 0 || i >= getSym().num(),
              "Wrong index for FGSTensor::contractAndAdd");

  TL_RAISE_IF(getSym()[i] != col.dimen() || tdims.getNVS()[i] != col.nvar(),
              "Wrong dimensions for FGSTensor::contractAndAdd");

  // set ‘sym_left’ and ‘sym_right’ to symmetries around ‘i’
  /* Here we have a symmetry of ‘this’ tensor and we have to set
     ‘sym_left’ to the subsymmetry left from the i-th variable and
     ‘sym_right’ to the subsymmetry right from the i-th variable. So we
     copy first all the symmetry and then put zeros to the left for
     ‘sym_right’ and to the right for ‘sym_left’. */
  Symmetry sym_left(getSym());
  Symmetry sym_right(getSym());
  for (int j = 0; j < getSym().num(); j++)
    {
      if (j <= i)
        sym_right[j] = 0;
      if (j >= i)
        sym_left[j] = 0;
    }

  int dleft = TensorDimens(sym_left, tdims.getNVS()).calcFoldMaxOffset();
  int dright = TensorDimens(sym_right, tdims.getNVS()).calcFoldMaxOffset();
  KronProdAll kp(3);
  kp.setUnit(0, dleft);
  kp.setMat(1, col);
  kp.setUnit(2, dright);
  FGSTensor tmp(out.nrows(), out.getDims());
  kp.mult(*this, tmp);
  out.add(1.0, tmp);
}

/* Here we go through folded tensor, and each index we convert to index
   of the unfolded tensor and copy the data to the unfolded. Then we
   unfold data within the unfolded tensor. */

UGSTensor::UGSTensor(const FGSTensor &ft)
  : UTensor(indor::along_col, ft.tdims.getNVX(), ft.nrows(),
            ft.tdims.calcUnfoldMaxOffset(), ft.dimen()),
    tdims(ft.tdims)
{
  for (index fi = ft.begin(); fi != ft.end(); ++fi)
    {
      index ui(*this, fi.getCoor());
      copyColumn(ft, *fi, *ui);
    }
  unfoldData();
}

// UGSTensor slicing from FSSparseTensor
/* This makes a folded slice from the sparse tensor and unfolds it. */
UGSTensor::UGSTensor(const FSSparseTensor &t, const IntSequence &ss,
                     const IntSequence &coor, TensorDimens td)
  : UTensor(indor::along_col, td.getNVX(), t.nrows(),
            td.calcUnfoldMaxOffset(), td.dimen()),
    tdims(std::move(td))
{
  if (ncols() == 0)
    return;

  FGSTensor ft(t, ss, coor, td);
  for (index fi = ft.begin(); fi != ft.end(); ++fi)
    {
      index ui(*this, fi.getCoor());
      copyColumn(ft, *fi, *ui);
    }
  unfoldData();
}

// UGSTensor slicing from UFSTensor
/* This makes a folded slice from dense and unfolds it. */
UGSTensor::UGSTensor(const UFSTensor &t, const IntSequence &ss,
                     const IntSequence &coor, TensorDimens td)
  : UTensor(indor::along_col, td.getNVX(), t.nrows(),
            td.calcUnfoldMaxOffset(), td.dimen()),
    tdims(std::move(td))
{
  FFSTensor folded(t);
  FGSTensor ft(folded, ss, coor, td);
  for (index fi = ft.begin(); fi != ft.end(); ++fi)
    {
      index ui(*this, fi.getCoor());
      copyColumn(ft, *fi, *ui);
    }
  unfoldData();
}

// UGSTensor increment and decrement codes
/* Clear, just call UTensor static methods. */
void
UGSTensor::increment(IntSequence &v) const
{
  TL_RAISE_IF(v.size() != dimen(),
              "Wrong input/output vector size in UGSTensor::increment");

  UTensor::increment(v, tdims.getNVX());
}

void
UGSTensor::decrement(IntSequence &v) const
{
  TL_RAISE_IF(v.size() != dimen(),
              "Wrong input/output vector size in UGSTensor::decrement");

  UTensor::decrement(v, tdims.getNVX());
}

/* Return a new instance of folded version. */
std::unique_ptr<FTensor>
UGSTensor::fold() const
{
  return std::make_unique<FGSTensor>(*this);
}

/* Return an offset of a given index. */
int
UGSTensor::getOffset(const IntSequence &v) const
{
  TL_RAISE_IF(v.size() != dimen(),
              "Wrong input vector size in UGSTensor::getOffset");

  return UTensor::getOffset(v, tdims.getNVX());
}

/* Unfold all data. We go through all the columns and for each we
   obtain an index of the first equivalent, and copy the data. */

void
UGSTensor::unfoldData()
{
  for (index in = begin(); in != end(); ++in)
    copyColumn(*getFirstIndexOf(in), *in);
}

/* Here we return the first index which is equivalent in the symmetry
   to the given index. It is a matter of sorting all the symmetry
   partitions of the index. */

Tensor::index
UGSTensor::getFirstIndexOf(const index &in) const
{
  IntSequence v(in.getCoor());
  int last = 0;
  for (int i = 0; i < tdims.getSym().num(); i++)
    {
      IntSequence vtmp(v, last, last+tdims.getSym()[i]);
      vtmp.sort();
      last += tdims.getSym()[i];
    }
  return index(*this, v);
}

/* Here is perfectly same code with the same semantics as in
   FGSTensor::contractAndAdd(). */

void
UGSTensor::contractAndAdd(int i, UGSTensor &out,
                          const URSingleTensor &col) const
{
  TL_RAISE_IF(i < 0 || i >= getSym().num(),
              "Wrong index for UGSTensor::contractAndAdd");
  TL_RAISE_IF(getSym()[i] != col.dimen() || tdims.getNVS()[i] != col.nvar(),
              "Wrong dimensions for UGSTensor::contractAndAdd");

  // set ‘sym_left’ and ‘sym_right’ to symmetries around i
  /* Same code as in FGSTensor::contractAndAdd */
  Symmetry sym_left(getSym());
  Symmetry sym_right(getSym());
  for (int j = 0; j < getSym().num(); j++)
    {
      if (j <= i)
        sym_right[j] = 0;
      if (j >= i)
        sym_left[j] = 0;
    }

  int dleft = TensorDimens(sym_left, tdims.getNVS()).calcUnfoldMaxOffset();
  int dright = TensorDimens(sym_right, tdims.getNVS()).calcUnfoldMaxOffset();
  KronProdAll kp(3);
  kp.setUnit(0, dleft);
  kp.setMat(1, col);
  kp.setUnit(2, dright);
  UGSTensor tmp(out.nrows(), out.getDims());
  kp.mult(*this, tmp);
  out.add(1.0, tmp);
}
