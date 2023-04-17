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

#include "stack_container.hh"
#include "pyramid_prod2.hh"
#include "ps_tensor.hh"

#include <memory>

// FoldedStackContainer::multAndAdd() sparse code
/* Here we multiply the sparse tensor with the FoldedStackContainer. We have
   four implementations, multAndAddSparse1(), multAndAddSparse2(),
   multAndAddSparse3(), and multAndAddSparse4(). The third is not threaded yet
   and I expect that it is certainly the slowest. The multAndAddSparse4()
   exploits the sparsity, however, it seems to be still worse than
   multAndAddSparse2() even for really sparse matrices. On the other hand, it
   can be more efficient than multAndAddSparse2() for large problems, since it
   does not need that much of memory and can avoid much swapping. Very
   preliminary examination shows that multAndAddSparse2() is the best in terms
   of time. */
void
FoldedStackContainer::multAndAdd(const FSSparseTensor &t,
                                 FGSTensor &out) const
{
  TL_RAISE_IF(t.nvar() != getAllSize(),
              "Wrong number of variables of tensor for FoldedStackContainer::multAndAdd");
  multAndAddSparse2(t, out);
}

// FoldedStackContainer::multAndAdd() dense code
/* Here we perform the Faà Di Bruno step for a given dimension ‘dim’, and for
   the dense fully symmetric tensor which is scattered in the container
   of general symmetric tensors. The implementation is pretty the same as
   UnfoldedStackContainer::multAndAdd() dense code. */
void
FoldedStackContainer::multAndAdd(int dim, const FGSContainer &c, FGSTensor &out) const
{
  TL_RAISE_IF(c.num() != numStacks(),
              "Wrong symmetry length of container for FoldedStackContainer::multAndAdd");

  sthread::detach_thread_group gr;

  for (auto &si : SymmetrySet(dim, c.num()))
    if (c.check(si))
      gr.insert(std::make_unique<WorkerFoldMAADense>(*this, si, c, out));

  gr.run();
}

/* This is analogous to WorkerUnfoldMAADense::operator()() code. */

void
WorkerFoldMAADense::operator()(std::mutex &mut)
{
  Permutation iden(dense_cont.num());
  IntSequence coor(iden.getMap().unfold(sym));
  const FGSTensor &g = dense_cont.get(sym);
  cont.multAndAddStacks(coor, g, out, mut);
}

WorkerFoldMAADense::WorkerFoldMAADense(const FoldedStackContainer &container,
                                       Symmetry s,
                                       const FGSContainer &dcontainer,
                                       FGSTensor &outten)
  : cont(container), sym(std::move(s)), dense_cont(dcontainer), out(outten)
{
}

/* This is analogous to UnfoldedStackContainer::multAndAddSparse1() code. */
void
FoldedStackContainer::multAndAddSparse1(const FSSparseTensor &t,
                                        FGSTensor &out) const
{
  sthread::detach_thread_group gr;
  UFSTensor dummy(0, numStacks(), t.dimen());
  for (Tensor::index ui = dummy.begin(); ui != dummy.end(); ++ui)
    gr.insert(std::make_unique<WorkerFoldMAASparse1>(*this, t, out, ui.getCoor()));

  gr.run();
}

/* This is analogous to WorkerUnfoldMAASparse1::operator()() code.
   The only difference is that instead of UPSTensor as a
   result of multiplication of unfolded tensor and tensors from
   containers, we have FPSTensor with partially folded permuted
   symmetry.

   TODO: make slice vertically narrowed according to the fill of t,
   vertically narrow out accordingly. */

void
WorkerFoldMAASparse1::operator()(std::mutex &mut)
{
  const EquivalenceSet &eset = TLStatic::getEquiv(out.dimen());
  const PermutationSet &pset = TLStatic::getPerm(t.dimen());
  Permutation iden(t.dimen());

  UPSTensor slice(t, cont.getStackSizes(), coor,
                  PerTensorDimens(cont.getStackSizes(), coor));
  for (int iper = 0; iper < pset.getNum(); iper++)
    {
      const Permutation &per = pset.get(iper);
      IntSequence percoor(coor.size());
      per.apply(coor, percoor);
      for (const auto &it : eset)
        {
          if (it.numClasses() == t.dimen())
            {
              StackProduct<FGSTensor> sp(cont, it, out.getSym());
              if (!sp.isZero(percoor))
                {
                  KronProdStack<FGSTensor> kp(sp, percoor);
                  kp.optimizeOrder();
                  const Permutation &oper = kp.getPer();
                  if (Permutation(oper, per) == iden)
                    {
                      FPSTensor fps(out.getDims(), it, slice, kp);
                      {
                        std::unique_lock<std::mutex> lk{mut};
                        fps.addTo(out);
                      }
                    }
                }
            }
        }
    }
}

WorkerFoldMAASparse1::WorkerFoldMAASparse1(const FoldedStackContainer &container,
                                           const FSSparseTensor &ten,
                                           FGSTensor &outten, IntSequence c)
  : cont(container), t(ten), out(outten), coor(std::move(c))
{
}

/* Here is the second implementation of sparse folded multAndAdd(). It
   is pretty similar to implementation of
   UnfoldedStackContainer::multAndAddSparse2() code. We make a
   dense folded slice, and then call folded multAndAddStacks(), which
   multiplies all the combinations compatible with the slice. */

void
FoldedStackContainer::multAndAddSparse2(const FSSparseTensor &t,
                                        FGSTensor &out) const
{
  sthread::detach_thread_group gr;
  FFSTensor dummy_f(0, numStacks(), t.dimen());
  for (Tensor::index fi = dummy_f.begin(); fi != dummy_f.end(); ++fi)
    gr.insert(std::make_unique<WorkerFoldMAASparse2>(*this, t, out, fi.getCoor()));

  gr.run();
}

/* Here we make a sparse slice first and then call multAndAddStacks()
   if the slice is not empty. If the slice is really sparse, we call
   sparse version of multAndAddStacks(). What means “really sparse” is
   given by ‘fill_threshold’. It is not tuned yet, a practice shows that
   it must be a really low number, since sparse multAndAddStacks() is
   much slower than the dense version.

   Further, we take only nonzero rows of the slice, and accordingly of
   the out tensor. We jump over zero initial rows and drop zero tailing
   rows. */

void
WorkerFoldMAASparse2::operator()(std::mutex &mut)
{
  GSSparseTensor slice(t, cont.getStackSizes(), coor,
                       TensorDimens(cont.getStackSizes(), coor));
  if (slice.getNumNonZero())
    {
      if (slice.getUnfoldIndexFillFactor() > FoldedStackContainer::fill_threshold)
        {
          FGSTensor dense_slice(slice);
          int r1 = slice.getFirstNonZeroRow();
          int r2 = slice.getLastNonZeroRow();
          FGSTensor dense_slice1(r1, r2-r1+1, dense_slice);
          FGSTensor out1(r1, r2-r1+1, out);
          cont.multAndAddStacks(coor, dense_slice1, out1, mut);
        }
      else
        cont.multAndAddStacks(coor, slice, out, mut);
    }
}

WorkerFoldMAASparse2::WorkerFoldMAASparse2(const FoldedStackContainer &container,
                                           const FSSparseTensor &ten,
                                           FGSTensor &outten, IntSequence c)
  : cont(container), t(ten), out(outten), coor(std::move(c))
{
}

/* Here is the third implementation of the sparse folded
   multAndAdd(). It is column-wise implementation, and thus is not a good
   candidate for the best performer.

   We go through all columns from the output. For each column we
   calculate folded ‘sumcol’ which is a sum of all appropriate columns
   for all suitable equivalences. So we go through all suitable
   equivalences, for each we construct a StackProduct object and
   construct IrregTensor for a corresponding column of z. The
   IrregTensor is an abstraction for Kronecker multiplication of
   stacked columns of the two containers without zeros. Then the column
   is added to ‘sumcol’. Finally, the ‘sumcol’ is multiplied by the
   sparse tensor. */

void
FoldedStackContainer::multAndAddSparse3(const FSSparseTensor &t,
                                        FGSTensor &out) const
{
  const EquivalenceSet &eset = TLStatic::getEquiv(out.dimen());
  for (Tensor::index run = out.begin(); run != out.end(); ++run)
    {
      Vector outcol{out.getCol(*run)};
      FRSingleTensor sumcol(t.nvar(), t.dimen());
      sumcol.zeros();
      for (const auto &it : eset)
        if (it.numClasses() == t.dimen())
          {
            StackProduct<FGSTensor> sp(*this, it, out.getSym());
            IrregTensorHeader header(sp, run.getCoor());
            IrregTensor irten(header);
            irten.addTo(sumcol);
          }
      t.multColumnAndAdd(sumcol, outcol);
    }
}

/* Here is the fourth implementation of sparse
   FoldedStackContainer::multAndAdd(). It is almost equivalent to
   multAndAddSparse2() with the exception that the FPSTensor as a
   result of a product of a slice and Kronecker product of the stack
   derivatives is calculated in the sparse fashion. For further details, see
   FoldedStackContainer::multAndAddStacks() sparse code and
   FPSTensor| sparse constructor. */

void
FoldedStackContainer::multAndAddSparse4(const FSSparseTensor &t, FGSTensor &out) const
{
  sthread::detach_thread_group gr;
  FFSTensor dummy_f(0, numStacks(), t.dimen());
  for (Tensor::index fi = dummy_f.begin(); fi != dummy_f.end(); ++fi)
    gr.insert(std::make_unique<WorkerFoldMAASparse4>(*this, t, out, fi.getCoor()));

  gr.run();
}

/* The WorkerFoldMAASparse4 is the same as WorkerFoldMAASparse2
   with the exception that we call a sparse version of
   multAndAddStacks(). */

void
WorkerFoldMAASparse4::operator()(std::mutex &mut)
{
  GSSparseTensor slice(t, cont.getStackSizes(), coor,
                       TensorDimens(cont.getStackSizes(), coor));
  if (slice.getNumNonZero())
    cont.multAndAddStacks(coor, slice, out, mut);
}

WorkerFoldMAASparse4::WorkerFoldMAASparse4(const FoldedStackContainer &container,
                                           const FSSparseTensor &ten,
                                           FGSTensor &outten, IntSequence c)
  : cont(container), t(ten), out(outten), coor(std::move(c))
{
}

// FoldedStackContainer::multAndAddStacks() dense code
/* This is almost the same as UnfoldedStackContainer::multAndAddStacks() code.
   The only difference is that we do not construct a UPSTensor from
   KronProdStack, but we construct partially folded permuted symmetry
   FPSTensor. Note that the tensor ‘g’ must be unfolded in order to be able to
   multiply with unfolded rows of Kronecker product. However, columns of such a
   product are partially folded giving a rise to the FPSTensor. */
void
FoldedStackContainer::multAndAddStacks(const IntSequence &coor,
                                       const FGSTensor &g,
                                       FGSTensor &out, std::mutex &mut) const
{
  const EquivalenceSet &eset = TLStatic::getEquiv(out.dimen());

  UGSTensor ug(g);
  UFSTensor dummy_u(0, numStacks(), g.dimen());
  for (Tensor::index ui = dummy_u.begin(); ui != dummy_u.end(); ++ui)
    {
      IntSequence tmp(ui.getCoor());
      tmp.sort();
      if (tmp == coor)
        {
          Permutation sort_per(ui.getCoor());
          sort_per.inverse();
          for (const auto &it : eset)
            if (it.numClasses() == g.dimen())
              {
                StackProduct<FGSTensor> sp(*this, it, sort_per, out.getSym());
                if (!sp.isZero(coor))
                  {
                    KronProdStack<FGSTensor> kp(sp, coor);
                    if (ug.getSym().isFull())
                      kp.optimizeOrder();
                    FPSTensor fps(out.getDims(), it, sort_per, ug, kp);
                    {
                      std::unique_lock<std::mutex> lk{mut};
                      fps.addTo(out);
                    }
                  }
              }
        }
    }
}

// FoldedStackContainer::multAndAddStacks() sparse code
/* This is almost the same as FoldedStackContainer::multAndAddStacks() dense code. The only
   difference is that the Kronecker product of the stacks is multiplied
   with sparse slice GSSparseTensor (not dense slice FGSTensor). The
   multiplication is done in FPSTensor sparse constructor. */
void
FoldedStackContainer::multAndAddStacks(const IntSequence &coor,
                                       const GSSparseTensor &g,
                                       FGSTensor &out, std::mutex &mut) const
{
  const EquivalenceSet &eset = TLStatic::getEquiv(out.dimen());
  UFSTensor dummy_u(0, numStacks(), g.dimen());
  for (Tensor::index ui = dummy_u.begin(); ui != dummy_u.end(); ++ui)
    {
      IntSequence tmp(ui.getCoor());
      tmp.sort();
      if (tmp == coor)
        {
          Permutation sort_per(ui.getCoor());
          sort_per.inverse();
          for (const auto &it : eset)
            if (it.numClasses() == g.dimen())
              {
                StackProduct<FGSTensor> sp(*this, it, sort_per, out.getSym());
                if (!sp.isZero(coor))
                  {
                    KronProdStack<FGSTensor> kp(sp, coor);
                    FPSTensor fps(out.getDims(), it, sort_per, g, kp);
                    {
                      std::unique_lock<std::mutex> lk{mut};
                      fps.addTo(out);
                    }
                  }
              }
        }
    }
}

// UnfoldedStackContainer::multAndAdd() sparse code
/*  Here we simply call either multAndAddSparse1() or
    multAndAddSparse2(). The first one allows for optimization of
    Kronecker products, so it seems to be more efficient. */
void
UnfoldedStackContainer::multAndAdd(const FSSparseTensor &t,
                                   UGSTensor &out) const
{
  TL_RAISE_IF(t.nvar() != getAllSize(),
              "Wrong number of variables of tensor for UnfoldedStackContainer::multAndAdd");
  multAndAddSparse2(t, out);
}

// UnfoldedStackContainer::multAndAdd() dense code
/* Here we implement the formula for stacks for fully symmetric tensor
   scattered in a number of general symmetry tensors contained in a given
   container. The implementations is pretty the same as in
   multAndAddSparse2() but we do not do the slices of sparse tensor, but
   only a lookup to the container.

   This means that we do not iterate through a dummy folded tensor to
   obtain folded coordinates of stacks, rather we iterate through all
   symmetries contained in the container and the coordinates of stacks
   are obtained as unfolded identity sequence via the symmetry. The
   reason of doing this is that we are unable to calculate symmetry from
   stack coordinates as easily as stack coordinates from the symmetry. */
void
UnfoldedStackContainer::multAndAdd(int dim, const UGSContainer &c,
                                   UGSTensor &out) const
{
  TL_RAISE_IF(c.num() != numStacks(),
              "Wrong symmetry length of container for UnfoldedStackContainer::multAndAdd");

  sthread::detach_thread_group gr;
  for (auto &si : SymmetrySet(dim, c.num()))
    if (c.check(si))
      gr.insert(std::make_unique<WorkerUnfoldMAADense>(*this, si, c, out));

  gr.run();
}

void
WorkerUnfoldMAADense::operator()(std::mutex &mut)
{
  Permutation iden(dense_cont.num());
  IntSequence coor(iden.getMap().unfold(sym));
  const UGSTensor &g = dense_cont.get(sym);
  cont.multAndAddStacks(coor, g, out, mut);
}

WorkerUnfoldMAADense::WorkerUnfoldMAADense(const UnfoldedStackContainer &container,
                                           Symmetry s,
                                           const UGSContainer &dcontainer,
                                           UGSTensor &outten)
  : cont(container), sym(std::move(s)), dense_cont(dcontainer), out(outten)
{
}

/* Here we implement the formula for unfolded tensors. If, for instance,
   a coordinate z of a tensor [f_z²] is partitioned as
   z=(a, b), then we perform the following:

               ⎛a_c(x)⎞ ⎛a_c(y)⎞
    [f_z²] · ∑ ⎢      ⎥⊗⎢      ⎥ = [f_aa] · ∑ a_c(x)⊗a_c(y) + [f_ab] · ∑ a_c(x)⊗b_c(y)
             ᶜ ⎝b_c(x)⎠ ⎝b_c(y)⎠            ᶜ                          ᶜ

                                 + [f_ba] · ∑ b_c(x)⊗a_c(y) + [f_bb] · ∑ b_c(x)⊗b_c(y)

   This is exactly what happens here. The code is clear. It goes through
   all combinations of stacks, and each thread is responsible for
   operation for the slice corresponding to the combination of the stacks. */

void
UnfoldedStackContainer::multAndAddSparse1(const FSSparseTensor &t,
                                          UGSTensor &out) const
{
  sthread::detach_thread_group gr;
  UFSTensor dummy(0, numStacks(), t.dimen());
  for (Tensor::index ui = dummy.begin(); ui != dummy.end(); ++ui)
    gr.insert(std::make_unique<WorkerUnfoldMAASparse1>(*this, t, out, ui.getCoor()));

  gr.run();
}

/* This does a step of UnfoldedStackContainer::multAndAddSparse1() for
   a given coordinates. First it makes the slice of the given stack coordinates.
   Then it multiplies everything what should be multiplied with the slice.
   That is it goes through all equivalences, creates StackProduct, then
   KronProdStack, which is added to ‘out’. So far everything is clear.

   However, we want to use optimized KronProdAllOptim to minimize a number of
   flops and memory needed in the Kronecker product. So we go through all
   permutations ‘per’, permute the coordinates to get ‘percoor’, go through all
   equivalences, and make KronProdStack and optimize it. The result of
   optimization is a permutation ‘oper’. Now, we multiply the Kronecker product
   with the slice, only if the slice has the same ordering of coordinates as
   the Kronecker product KronProdStack. However, it is not perfectly true.
   Since we go through *all* permutations ‘per’, there might be two different
   permutations leading to the same ordering in KronProdStack and thus the same
   ordering in the optimized KronProdStack. The two cases would be counted
   twice, which is wrong. That is why we do not condition on
   coor∘oper∘per = coor, but we condition on oper∘per = id. In this way, we
   rule out permutations ‘per’ leading to the same ordering of stacks when
   applied on ‘coor’.

   TODO: vertically narrow slice and out according to the fill in t. */

void
WorkerUnfoldMAASparse1::operator()(std::mutex &mut)
{
  const EquivalenceSet &eset = TLStatic::getEquiv(out.dimen());
  const PermutationSet &pset = TLStatic::getPerm(t.dimen());
  Permutation iden(t.dimen());

  UPSTensor slice(t, cont.getStackSizes(), coor,
                  PerTensorDimens(cont.getStackSizes(), coor));
  for (int iper = 0; iper < pset.getNum(); iper++)
    {
      const Permutation &per = pset.get(iper);
      IntSequence percoor(coor.size());
      per.apply(coor, percoor);
      for (const auto &it : eset)
        if (it.numClasses() == t.dimen())
          {
            StackProduct<UGSTensor> sp(cont, it, out.getSym());
            if (!sp.isZero(percoor))
              {
                KronProdStack<UGSTensor> kp(sp, percoor);
                kp.optimizeOrder();
                const Permutation &oper = kp.getPer();
                if (Permutation(oper, per) == iden)
                  {
                    UPSTensor ups(out.getDims(), it, slice, kp);
                    {
                      std::unique_lock<std::mutex> lk{mut};
                      ups.addTo(out);
                    }
                  }
              }
          }
    }
}

WorkerUnfoldMAASparse1::WorkerUnfoldMAASparse1(const UnfoldedStackContainer &container,
                                               const FSSparseTensor &ten,
                                               UGSTensor &outten, IntSequence c)
  : cont(container), t(ten), out(outten), coor(std::move(c))
{
}

/* In here we implement the formula by a bit different way. We use the
   fact, using notation of UnfoldedStackContainer::multAndAddSparse2(),
   that
                                      ⎛               ⎞
   [f_ba] · ∑ b_c(x)⊗a_c(y) = [f_ba] ·⎢∑ a_c(y)⊗b_c(x)⎥· P
            ᶜ                         ⎝ᶜ              ⎠

   where P is a suitable permutation of columns. The permutation
   corresponds to (in this example) a swap of a and b. An advantage
   of this approach is that we do not need UPSTensor for [f_ba], and
   thus we decrease the number of needed slices.

   So we go through all folded indices of stack coordinates, then for
   each such index ‘fi’ we make a slice and call multAndAddStacks(). This
   goes through all corresponding unfolded indices to perform the
   formula. Each unsorted (unfold) index implies a sorting permutation
   ‘sort_per’ which must be used to permute stacks in StackProduct, and
   permute equivalence classes when UPSTensor is formed. In this way
   the column permutation P from the formula is factored to the
   permutation of UPSTensor. */

void
UnfoldedStackContainer::multAndAddSparse2(const FSSparseTensor &t,
                                          UGSTensor &out) const
{
  sthread::detach_thread_group gr;
  FFSTensor dummy_f(0, numStacks(), t.dimen());
  for (Tensor::index fi = dummy_f.begin(); fi != dummy_f.end(); ++fi)
    gr.insert(std::make_unique<WorkerUnfoldMAASparse2>(*this, t, out, fi.getCoor()));

  gr.run();
}

/* This does a step of UnfoldedStackContainer::multAndAddSparse2() for a given
   coordinates.

   TODO: implement multAndAddStacks() for sparse slice as
   FoldedStackContainer::multAndAddStacks() sparse code and do this method as
   WorkerFoldMAASparse2::operator()(). */

void
WorkerUnfoldMAASparse2::operator()(std::mutex &mut)
{
  GSSparseTensor slice(t, cont.getStackSizes(), coor,
                       TensorDimens(cont.getStackSizes(), coor));
  if (slice.getNumNonZero())
    {
      FGSTensor fslice(slice);
      UGSTensor dense_slice(fslice);
      int r1 = slice.getFirstNonZeroRow();
      int r2 = slice.getLastNonZeroRow();
      UGSTensor dense_slice1(r1, r2-r1+1, dense_slice);
      UGSTensor out1(r1, r2-r1+1, out);

      cont.multAndAddStacks(coor, dense_slice1, out1, mut);
    }
}

WorkerUnfoldMAASparse2::WorkerUnfoldMAASparse2(const UnfoldedStackContainer &container,
                                               const FSSparseTensor &ten,
                                               UGSTensor &outten, IntSequence c)
  : cont(container), t(ten), out(outten), coor(std::move(c))
{
}

/* For a given unfolded coordinates of stacks ‘fi’, and appropriate
   tensor g, whose symmetry is a symmetry of ‘fi’, the method
   contributes to ‘out’ all tensors in unfolded stack formula involving
   stacks chosen by ‘fi’.

   We go through all ‘ui’ coordinates which yield ‘fi’ after sorting. We
   construct a permutation ‘sort_per’ which sorts ‘ui’ to ‘fi’. We go
   through all appropriate equivalences, and construct StackProduct
   from equivalence classes permuted by ‘sort_per’, then UPSTensor with
   implied permutation of columns by the permuted equivalence by
   ‘sort_per’. The UPSTensor is then added to ‘out’.

   We cannot use here the optimized KronProdStack, since the symmetry
   of ‘UGSTensor& g’ prescribes the ordering of the stacks. However, if
   ‘g’ is fully symmetric, we can do the optimization harmlessly. */

void
UnfoldedStackContainer::multAndAddStacks(const IntSequence &fi,
                                         const UGSTensor &g,
                                         UGSTensor &out, std::mutex &mut) const
{
  const EquivalenceSet &eset = TLStatic::getEquiv(out.dimen());

  UFSTensor dummy_u(0, numStacks(), g.dimen());
  for (Tensor::index ui = dummy_u.begin(); ui != dummy_u.end(); ++ui)
    {
      IntSequence tmp(ui.getCoor());
      tmp.sort();
      if (tmp == fi)
        {
          Permutation sort_per(ui.getCoor());
          sort_per.inverse();
          for (const auto &it : eset)
            if (it.numClasses() == g.dimen())
              {
                StackProduct<UGSTensor> sp(*this, it, sort_per, out.getSym());
                if (!sp.isZero(fi))
                  {
                    KronProdStack<UGSTensor> kp(sp, fi);
                    if (g.getSym().isFull())
                      kp.optimizeOrder();
                    UPSTensor ups(out.getDims(), it, sort_per, g, kp);
                    {
                      std::unique_lock<std::mutex> lk{mut};
                      ups.addTo(out);
                    }
                  }
              }
        }
    }
}
