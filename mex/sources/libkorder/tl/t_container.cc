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

#include "t_container.hh"
#include "kron_prod.hh"
#include "ps_tensor.hh"
#include "pyramid_prod.hh"

// UGSContainer conversion from FGSContainer
UGSContainer::UGSContainer(const FGSContainer &c)
  : TensorContainer<UGSTensor>(c.num())
{
  for (const auto &it : c)
    insert(std::make_unique<UGSTensor>(*(it.second)));
}

/* We set ‘l’ to dimension of ‘t’, this is a tensor which multiplies tensors
   from the container from the left. Also we set ‘k’ to a dimension of the
   resulting tensor. We go through all equivalences on ‘k’ element set and
   pickup only those which have ‘l’ classes.

   In each loop, we fetch all necessary tensors for the product to the vector
   ‘ts’. Then we form Kronecker product KronProdAll and feed it with tensors
   from ‘ts’. Then we form unfolded permuted symmetry tensor UPSTensor as
   matrix product of ‘t’ and Kronecker product ‘kp’. Then we add the permuted
   data to ‘out’. This is done by UPSTensor method addTo(). */

void
UGSContainer::multAndAdd(const UGSTensor &t, UGSTensor &out) const
{
  int l = t.dimen();
  int k = out.dimen();
  const EquivalenceSet &eset = TLStatic::getEquiv(k);

  for (const auto &it : eset)
    if (it.numClasses() == l)
      {
        std::vector<const UGSTensor *> ts = fetchTensors(out.getSym(), it);
        KronProdAllOptim kp(l);
        for (int i = 0; i < l; i++)
          kp.setMat(i, *(ts[i]));
        kp.optimizeOrder();
        UPSTensor ups(out.getDims(), it, t, kp);
        ups.addTo(out);
      }
}

// FGSContainer conversion from UGSContainer
FGSContainer::FGSContainer(const UGSContainer &c)
  : TensorContainer<FGSTensor>(c.num())
{
  for (const auto &it : c)
    insert(std::make_unique<FGSTensor>(*(it.second)));
}

// FGSContainer::multAndAdd() folded code
/* Here we perform one step of the Faà Di Bruno operation. We call the
   multAndAdd() for unfolded tensor. */
void
FGSContainer::multAndAdd(const FGSTensor &t, FGSTensor &out) const
{
  UGSTensor ut(t);
  multAndAdd(ut, out);
}

// FGSContainer::multAndAdd() unfolded code
/* This is the same as UGSContainer::multAndAdd() but we do not construct
   UPSTensor from the Kronecker product, but FPSTensor. */
void
FGSContainer::multAndAdd(const UGSTensor &t, FGSTensor &out) const
{
  int l = t.dimen();
  int k = out.dimen();
  const EquivalenceSet &eset = TLStatic::getEquiv(k);

  for (const auto &it : eset)
    if (it.numClasses() == l)
      {
        std::vector<const FGSTensor *> ts
          = fetchTensors(out.getSym(), it);
        KronProdAllOptim kp(l);
        for (int i = 0; i < l; i++)
          kp.setMat(i, *(ts[i]));
        kp.optimizeOrder();
        FPSTensor fps(out.getDims(), it, t, kp);
        fps.addTo(out);
      }
}

/* This fills a given vector with integer sequences corresponding to first
   ‘num’ indices from interval ‘start’ (including) to ‘end’ (excluding). If
   there are not ‘num’ of such indices, the shorter vector is returned. */
Tensor::index
FGSContainer::getIndices(int num, std::vector<IntSequence> &out,
                         const Tensor::index &start,
                         const Tensor::index &end)
{
  out.clear();
  int i = 0;
  Tensor::index run = start;
  while (i < num && run != end)
    {
      out.push_back(run.getCoor());
      i++;
      ++run;
    }
  return run;
}
