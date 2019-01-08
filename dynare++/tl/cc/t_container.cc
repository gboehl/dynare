// Copyright 2004, Ondra Kamenik

#include "t_container.hh"
#include "kron_prod.hh"
#include "ps_tensor.hh"
#include "pyramid_prod.hh"

const int FGSContainer::num_one_time = 10;

// |UGSContainer| conversion from |FGSContainer|
UGSContainer::UGSContainer(const FGSContainer &c)
  : TensorContainer<UGSTensor>(c.num())
{
  for (FGSContainer::const_iterator it = c.begin();
       it != c.end(); ++it)
    {
      UGSTensor *unfolded = new UGSTensor(*((*it).second));
      insert(unfolded);
    }
}

/* We set |l| to dimension of |t|, this is a tensor which multiplies
   tensors from the container from the left. Also we set |k| to a
   dimension of the resulting tensor. We go through all equivalences on
   |k| element set and pickup only those which have $l$ classes.

   In each loop, we fetch all necessary tensors for the product to the
   vector |ts|. Then we form Kronecker product |KronProdAll| and feed it
   with tensors from |ts|. Then we form unfolded permuted symmetry tensor
   |UPSTensor| as matrix product of |t| and Kronecker product |kp|. Then
   we add the permuted data to |out|. This is done by |UPSTensor| method
   |addTo|. */

void
UGSContainer::multAndAdd(const UGSTensor &t, UGSTensor &out) const
{
  int l = t.dimen();
  int k = out.dimen();
  const EquivalenceSet &eset = ebundle.get(k);

  for (EquivalenceSet::const_iterator it = eset.begin();
       it != eset.end(); ++it)
    {
      if ((*it).numClasses() == l)
        {
          vector<const UGSTensor *> ts
            = fetchTensors(out.getSym(), *it);
          KronProdAllOptim kp(l);
          for (int i = 0; i < l; i++)
            kp.setMat(i, *(ts[i]));
          kp.optimizeOrder();
          UPSTensor ups(out.getDims(), *it, t, kp);
          ups.addTo(out);
        }
    }
}

// |FGSContainer| conversion from |UGSContainer|
FGSContainer::FGSContainer(const UGSContainer &c)
  : TensorContainer<FGSTensor>(c.num())
{
  for (UGSContainer::const_iterator it = c.begin();
       it != c.end(); ++it)
    {
      FGSTensor *folded = new FGSTensor(*((*it).second));
      insert(folded);
    }
}

// |FGSContainer::multAndAdd| folded code
/* Here we perform one step of the Faa Di Bruno operation. We call the
   |multAndAdd| for unfolded tensor. */
void
FGSContainer::multAndAdd(const FGSTensor &t, FGSTensor &out) const
{
  UGSTensor ut(t);
  multAndAdd(ut, out);
}

// |FGSContainer::multAndAdd| unfolded code
/* This is the same as |@<|UGSContainer::multAndAdd| code@>|
   but we do not construct |UPSTensor| from the Kronecker
   product, but |FPSTensor|. */
void
FGSContainer::multAndAdd(const UGSTensor &t, FGSTensor &out) const
{
  int l = t.dimen();
  int k = out.dimen();
  const EquivalenceSet &eset = ebundle.get(k);

  for (EquivalenceSet::const_iterator it = eset.begin();
       it != eset.end(); ++it)
    {
      if ((*it).numClasses() == l)
        {
          vector<const FGSTensor *> ts
            = fetchTensors(out.getSym(), *it);
          KronProdAllOptim kp(l);
          for (int i = 0; i < l; i++)
            kp.setMat(i, *(ts[i]));
          kp.optimizeOrder();
          FPSTensor fps(out.getDims(), *it, t, kp);
          fps.addTo(out);
        }
    }
}

/* This fills a given vector with integer sequences corresponding to
   first |num| indices from interval |start| (including) to |end|
   (excluding). If there are not |num| of such indices, the shorter vector
   is returned. */
Tensor::index
FGSContainer::getIndices(int num, vector<IntSequence> &out,
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
