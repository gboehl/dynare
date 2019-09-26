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

#include <memory>

#include "normal_moments.hh"
#include "permutation.hh"
#include "kron_prod.hh"
#include "tl_static.hh"

UNormalMoments::UNormalMoments(int maxdim, const TwoDMatrix &v)
  : TensorContainer<URSingleTensor>(1)
{
  if (maxdim >= 2)
    generateMoments(maxdim, v);
}

/* Here we fill up the container with the tensors for d=2,4,6,… up to the given
   dimension. Each tensor of moments is equal to Fₙ(⊗ⁿv). This has a dimension
   equal to 2n. See the header file for proof and details.

   Here we sequentially construct the Kronecker power ⊗ⁿv and apply Fₙ.
 */
void
UNormalMoments::generateMoments(int maxdim, const TwoDMatrix &v)
{
  TL_RAISE_IF(v.nrows() != v.ncols(),
              "Variance-covariance matrix is not square in UNormalMoments constructor");

  int nv = v.nrows();
  auto mom2 = std::make_unique<URSingleTensor>(nv, 2);
  mom2->getData() = v.getData();
  insert(std::move(mom2));
  auto kronv = std::make_unique<URSingleTensor>(nv, 2);
  kronv->getData() = v.getData();
  for (int d = 4; d <= maxdim; d += 2)
    {
      auto newkronv = std::make_unique<URSingleTensor>(nv, d);
      KronProd::kronMult(ConstVector(v.getData()),
                         ConstVector(kronv->getData()),
                         newkronv->getData());
      kronv = std::move(newkronv);
      auto mom = std::make_unique<URSingleTensor>(nv, d);
      // apply Fₙ to ‘kronv’
      /* Here we go through all equivalences, select only those having 2
         elements in each class, then go through all elements in ‘kronv’ and
         add to permuted location of ‘mom’.

         The permutation must be taken as inverse of the permutation implied by
         the equivalence, since we need a permutation which after application
         to identity of indices yileds indices in the equivalence classes. Note
         how the Equivalence::apply() method works. */
      mom->zeros();
      const EquivalenceSet eset = TLStatic::getEquiv(d);
      for (const auto &cit : eset)
        if (selectEquiv(cit))
          {
            Permutation per(cit);
            per.inverse();
            for (Tensor::index it = kronv->begin(); it != kronv->end(); ++it)
              {
                IntSequence ind(kronv->dimen());
                per.apply(it.getCoor(), ind);
                Tensor::index it2(*mom, ind);
                mom->get(*it2, 0) += kronv->get(*it, 0);
              }
          }
      insert(std::move(mom));
    }
}

/* We return true for an equivalence whose each class has 2 elements. */

bool
UNormalMoments::selectEquiv(const Equivalence &e)
{
  if (2*e.numClasses() != e.getN())
    return false;
  for (const auto &si : e)
    if (si.length() != 2)
      return false;
  return true;
}

/* Here we go through all the unfolded container, fold each tensor and
   insert it. */
FNormalMoments::FNormalMoments(const UNormalMoments &moms)
  : TensorContainer<FRSingleTensor>(1)
{
  for (const auto &mom : moms)
    insert(std::make_unique<FRSingleTensor>(*(mom.second)));
}
