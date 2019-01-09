// Copyright 2004, Ondra Kamenik

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

/* Here we fill up the container with the tensors for $d=2,4,6,\ldots$
   up to the given dimension. Each tensor of moments is equal to
   $F_n\left(\otimes^nv\right).$ This has a dimension equal to
   $2n$. See the header file for proof and details.

   Here we sequentially construct the Kronecker power
   $\otimes^nv$, and apply $F_n$. */

void
UNormalMoments::generateMoments(int maxdim, const TwoDMatrix &v)
{
  TL_RAISE_IF(v.nrows() != v.ncols(),
              "Variance-covariance matrix is not square in UNormalMoments constructor");

  int nv = v.nrows();
  URSingleTensor *mom2 = new URSingleTensor(nv, 2);
  mom2->getData() = v.getData();
  insert(mom2);
  URSingleTensor *kronv = new URSingleTensor(nv, 2);
  kronv->getData() = v.getData();
  for (int d = 4; d <= maxdim; d += 2)
    {
      URSingleTensor *newkronv = new URSingleTensor(nv, d);
      KronProd::kronMult(ConstVector(v.getData()),
                         ConstVector(kronv->getData()),
                         newkronv->getData());
      delete kronv;
      kronv = newkronv;
      URSingleTensor *mom = new URSingleTensor(nv, d);
      // apply $F_n$ to |kronv|
      /* Here we go through all equivalences, select only those having 2
         elements in each class, then go through all elements in |kronv| and
         add to permuted location of |mom|.

         The permutation must be taken as inverse of the permutation implied by
         the equivalence, since we need a permutation which after application
         to identity of indices yileds indices in the equivalence classes. Note
         how the |Equivalence::apply| method works. */
      mom->zeros();
      const EquivalenceSet eset = ebundle.get(d);
      for (const auto & cit : eset)
        {
          if (selectEquiv(cit))
            {
              Permutation per(cit);
              per.inverse();
              for (Tensor::index it = kronv->begin(); it != kronv->end(); ++it)
                {
                  IntSequence ind(kronv->dimen());
                  per.apply(it.getCoor(), ind);
                  Tensor::index it2(mom, ind);
                  mom->get(*it2, 0) += kronv->get(*it, 0);
                }
            }
        }
      insert(mom);
    }
  delete kronv;
}

/* We return |true| for an equivalence whose each class has 2 elements. */

bool
UNormalMoments::selectEquiv(const Equivalence &e)
{
  if (2*e.numClasses() != e.getN())
    return false;
  for (const auto & si : e)
    {
      if (si.length() != 2)
        return false;
    }
  return true;
}

/* Here we go through all the unfolded container, fold each tensor and
   insert it. */
FNormalMoments::FNormalMoments(const UNormalMoments &moms)
  : TensorContainer<FRSingleTensor>(1)
{
  for (const auto & mom : moms)
    {
      FRSingleTensor *fm = new FRSingleTensor(*(mom.second));
      insert(fm);
    }
}
