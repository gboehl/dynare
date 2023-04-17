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

#include "pyramid_prod2.hh"
#include "rfs_tensor.hh"

/* Here we only call sp.createPackedColumns(c, cols, unit_flag) which
   fills ‘cols’ and ‘unit_flag’ for the given column ‘c’. Then we set
   ‘end_seq’ according to ‘unit_flag’ and columns lengths. */

IrregTensorHeader::IrregTensorHeader(const StackProduct<FGSTensor> &sp,
                                     const IntSequence &c)
  : nv(sp.getAllSize()),
    unit_flag(sp.dimen()),
    cols(sp.createPackedColumns(c, unit_flag)),
    end_seq(sp.dimen())
{
  for (int i = 0; i < sp.dimen(); i++)
    {
      end_seq[i] = cols[i]->length();
      if (unit_flag[i] != -1)
        end_seq[i] = unit_flag[i]+1;
    }
}

/* Here we have to increment the given integer sequence. We do it by
   the following code, whose pattern is valid for all tensor. The only
   difference is how we increment item of coordinates. */

void
IrregTensorHeader::increment(IntSequence &v) const
{
  TL_RAISE_IF(v.size() != dimen(),
              "Wrong size of coordinates in IrregTensorHeader::increment");

  if (v.size() == 0)
    return;
  int i = v.size()-1;

  // increment i-th item in coordinate ‘v’
  /* Here we increment item of coordinates. Whenever we reached end of
     column coming from matrices, and ‘unit_flag’ is not -1, we have to
     jump to that ‘unit_flag’. */
  v[i]++;
  if (unit_flag[i] != -1 && v[i] == cols[i]->length()-1)
    v[i] = unit_flag[i];

  while (i > 0 && v[i] == end_seq[i])
    {
      v[i] = 0;
      i--;
      // increment i-th item in coordinate ‘v’
      /* Same code as above */
      v[i]++;
      if (unit_flag[i] != -1 && v[i] == cols[i]->length()-1)
        v[i] = unit_flag[i];
    }
}

/* It is a product of all column lengths. */

int
IrregTensorHeader::calcMaxOffset() const
{
  int res = 1;
  for (int i = 0; i < dimen(); i++)
    res *= cols[i]->length();
  return res;
}

/* Everything is done in IrregTensorHeader, only we have to Kronecker
   multiply all columns of the header. */

IrregTensor::IrregTensor(const IrregTensorHeader &h)
  : Tensor(indor::along_row, IntSequence(h.dimen(), 0), h.end_seq,
           h.calcMaxOffset(), 1, h.dimen()),
    header(h)
{
  if (header.dimen() == 1)
    {
      getData() = *(header.cols[0]);
      return;
    }

  auto last = std::make_unique<Vector>(*(header.cols[header.dimen()-1]));
  for (int i = header.dimen()-2; i > 0; i--)
    {
      auto newlast = std::make_unique<Vector>(last->length()*header.cols[i]->length());
      KronProd::kronMult(ConstVector(*(header.cols[i])),
                         ConstVector(*last), *newlast);
      last = std::move(newlast);
    }
  KronProd::kronMult(ConstVector(*(header.cols[0])),
                     ConstVector(*last), getData());
}

void
IrregTensor::addTo(FRSingleTensor &out) const
{
  for (index it = begin(); it != end(); ++it)
    {
      IntSequence tmp(it.getCoor());
      tmp.sort();
      Tensor::index ind(out, tmp);
      out.get(*ind, 0) += get(*it, 0);
    }
}
