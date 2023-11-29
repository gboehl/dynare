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

#include "rfs_tensor.hh"
#include "kron_prod.hh"
#include "tl_exception.hh"

/* The conversion from unfolded to folded sums up all data from unfolded
   corresponding to one folded index. So we go through all the rows in the
   unfolded tensor ‘ut’, make an index of the folded tensor by sorting the
   coordinates, and add the row. */
FRTensor::FRTensor(const URTensor& ut) :
    FTensor(indor::along_row, IntSequence(ut.dimen(), ut.nvar()),
            FFSTensor::calcMaxOffset(ut.nvar(), ut.dimen()), ut.ncols(), ut.dimen()),
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

std::unique_ptr<UTensor>
FRTensor::unfold() const
{
  return std::make_unique<URTensor>(*this);
}

/* Incrementing is easy. The same as for FFSTensor. */

void
FRTensor::increment(IntSequence& v) const
{
  TL_RAISE_IF(v.size() != dimen(), "Wrong input/output vector size in FRTensor::increment");

  UTensor::increment(v, nv);
  v.monotone();
}

void
FRTensor::decrement(IntSequence& v) const
{
  TL_RAISE_IF(v.size() != dimen(), "Wrong input/output vector size in FRTensor::decrement");

  FTensor::decrement(v, nv);
}

/* Here we convert folded full symmetry tensor to unfolded. We copy all
   columns of folded tensor to unfolded and leave other columns
   (duplicates) zero. In this way, if the unfolded tensor is folded back,
   we should get the same data. */
URTensor::URTensor(const FRTensor& ft) :
    UTensor(indor::along_row, IntSequence(ft.dimen(), ft.nvar()),
            UFSTensor::calcMaxOffset(ft.nvar(), ft.dimen()), ft.ncols(), ft.dimen()),
    nv(ft.nvar())
{
  zeros();
  for (index src = ft.begin(); src != ft.end(); ++src)
    {
      index in(*this, src.getCoor());
      copyRow(ft, *src, *in);
    }
}

std::unique_ptr<FTensor>
URTensor::fold() const
{
  return std::make_unique<FRTensor>(*this);
}

void
URTensor::increment(IntSequence& v) const
{
  TL_RAISE_IF(v.size() != dimen(), "Wrong input/output vector size in URTensor::increment");

  UTensor::increment(v, nv);
}

void
URTensor::decrement(IntSequence& v) const
{
  TL_RAISE_IF(v.size() != dimen(), "Wrong input/output vector size in URTensor::decrement");

  UTensor::decrement(v, nv);
}

int
URTensor::getOffset(const IntSequence& v) const
{
  TL_RAISE_IF(v.size() != dimen(), "Wrong input vector size in URTensor::getOffset");

  return UTensor::getOffset(v, nv);
}

/* Here we construct v₁⊗v₂⊗…⊗vₙ, where v₁,v₂,…,vₙ are stored in a
   std::vector<ConstVector>. */

URSingleTensor::URSingleTensor(const std::vector<ConstVector>& cols) :
    URTensor(1, cols[0].length(), cols.size())
{
  if (dimen() == 1)
    {
      getData() = cols[0];
      return;
    }

  auto last = std::make_unique<Vector>(cols[cols.size() - 1]);
  for (int i = cols.size() - 2; i > 0; i--)
    {
      auto newlast = std::make_unique<Vector>(power(nvar(), cols.size() - i));
      KronProd::kronMult(cols[i], ConstVector(*last), *newlast);
      last = std::move(newlast);
    }
  KronProd::kronMult(cols[0], ConstVector(*last), getData());
}

/* Here we construct v⊗…⊗v, where ‘d’ gives the number of copies of v. */

URSingleTensor::URSingleTensor(const ConstVector& v, int d) : URTensor(1, v.length(), d)
{
  if (d == 1)
    {
      getData() = v;
      return;
    }

  auto last = std::make_unique<Vector>(v);
  for (int i = d - 2; i > 0; i--)
    {
      auto newlast = std::make_unique<Vector>(last->length() * v.length());
      KronProd::kronMult(v, ConstVector(*last), *newlast);
      last = std::move(newlast);
    }
  KronProd::kronMult(v, ConstVector(*last), getData());
}

std::unique_ptr<FTensor>
URSingleTensor::fold() const
{
  return std::make_unique<FRSingleTensor>(*this);
}

/* The conversion from unfolded URSingleTensor to folded FRSingleTensor is
   exactly the same as the conversion from URTensor to FRTensor, except that we
   do not copy rows but elements. */
FRSingleTensor::FRSingleTensor(const URSingleTensor& ut) : FRTensor(1, ut.nvar(), ut.dimen())
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
