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
 * along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "pyramid_prod.hh"
#include "permutation.hh"
#include "tl_exception.hh"

/* Here we construct the USubTensor object. We allocate space via the parent
   URTensor. Number of columns is a length of the list of indices ‘lst’, number
   of variables and dimensions are of the tensor $h$, this is given by ‘hdims’.

   We go through all equivalences with number of classes equal to dimension of
   B. For each equivalence we make a permutation ‘per’. Then we fetch all the
   necessary tensors g with symmetries implied by symmetry of B and the
   equivalence. Then we go through the list of indices, permute them by the
   permutation and add the Kronecker product of the selected columns. This is
   done by addKronColumn(). */

USubTensor::USubTensor(const TensorDimens &bdims,
                       const TensorDimens &hdims,
                       const FGSContainer &cont,
                       const std::vector<IntSequence> &lst)
  : URTensor(lst.size(), hdims.getNVX(0), hdims.dimen())
{
  TL_RAISE_IF(!hdims.getNVX().isConstant(),
              "Tensor has not full symmetry in USubTensor()");
  const EquivalenceSet &eset = TLStatic::getEquiv(bdims.dimen());
  zeros();
  for (const auto &it : eset)
    {
      if (it.numClasses() == hdims.dimen())
        {
          Permutation per(it);
          std::vector<const FGSTensor *> ts
            = cont.fetchTensors(bdims.getSym(), it);
          for (int i = 0; i < static_cast<int>(lst.size()); i++)
            {
              IntSequence perindex(lst[i].size());
              per.apply(lst[i], perindex);
              addKronColumn(i, ts, perindex);
            }
        }
    }
}

/* This makes a Kronecker product of appropriate columns from tensors in ‘fs’
   and adds such data to i-th column of this matrix. The appropriate columns
   are defined by ‘pindex’ sequence. A column of a tensor has index created
   from a corresponding part of ‘pindex’. The sizes of these parts are given by
   dimensions of the tensors in ‘ts’.

   Here we break the given index ‘pindex’ according to the dimensions of the
   tensors in ‘ts’, and for each subsequence of the ‘pindex’ we find an index
   of the folded tensor, which involves calling getOffset() for folded tensor,
   which might be costly. We gather all columns to a vector ‘tmpcols’ which are
   Kronecker multiplied in constructor of URSingleTensor. Finally we add data
   of URSingleTensor to the i-th column. */

void
USubTensor::addKronColumn(int i, const std::vector<const FGSTensor *> &ts,
                          const IntSequence &pindex)
{
  std::vector<ConstVector> tmpcols;
  int lastdim = 0;
  for (auto t : ts)
    {
      IntSequence ind(pindex, lastdim, lastdim+t->dimen());
      lastdim += t->dimen();
      index in(*t, ind);
      tmpcols.push_back(t->getCol(*in));
    }

  URSingleTensor kronmult(tmpcols);
  Vector coli{getCol(i)};
  coli.add(1.0, kronmult.getData());
}
