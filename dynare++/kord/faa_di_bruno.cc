/*
 * Copyright © 2005 Ondra Kamenik
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

#include "faa_di_bruno.hh"
#include "fine_container.hh"

#include <cmath>

// FaaDiBruno::calculate() folded sparse code
/* We take an opportunity to refine the stack container to avoid allocation of
   more memory than available. */
void
FaaDiBruno::calculate(const StackContainer<FGSTensor> &cont,
                      const TensorContainer<FSSparseTensor> &f,
                      FGSTensor &out)
{
  out.zeros();
  for (int l = 1; l <= out.dimen(); l++)
    {
      auto [max, mem_mb, p_size_mb] = estimRefinement(out.getDims(), out.nrows(), l);
      FoldedFineContainer fine_cont(cont, max);
      fine_cont.multAndAdd(l, f, out);
      JournalRecord recc(journal);
      recc << "dim=" << l << " avmem=" << mem_mb << " tmpmem=" << p_size_mb << " max=" << max
           << " stacks=" << cont.numStacks() << u8"→" << fine_cont.numStacks() << endrec;
    }
}

// FaaDiBruno::calculate() folded dense code
/* Here we just simply evaluate multAndAdd() for the dense container. There is
   no opportunity for tuning. */
void
FaaDiBruno::calculate(const FoldedStackContainer &cont, const FGSContainer &g,
                      FGSTensor &out)
{
  out.zeros();
  for (int l = 1; l <= out.dimen(); l++)
    {
      long int mem = SystemResources::availableMemory();
      cont.multAndAdd(l, g, out);
      JournalRecord rec(journal);
      int mem_mb = mem/1024/1024;
      rec << "dim=" << l << " avmem=" << mem_mb << endrec;
    }
}

// FaaDiBruno::calculate() unfolded sparse code
/* This is the same as FaaDiBruno::calculate() folded sparse code. The only
   difference is that we construct unfolded fine container. */
void
FaaDiBruno::calculate(const StackContainer<UGSTensor> &cont,
                      const TensorContainer<FSSparseTensor> &f,
                      UGSTensor &out)
{
  out.zeros();
  for (int l = 1; l <= out.dimen(); l++)
    {
      auto [max, mem_mb, p_size_mb] = estimRefinement(out.getDims(), out.nrows(), l);
      UnfoldedFineContainer fine_cont(cont, max);
      fine_cont.multAndAdd(l, f, out);
      JournalRecord recc(journal);
      recc << "dim=" << l << " avmem=" << mem_mb << " tmpmem=" << p_size_mb << " max=" << max
           << " stacks=" << cont.numStacks() << u8"→" << fine_cont.numStacks() << endrec;
    }
}

// FaaDiBruno::calculate() unfolded dense code
/* Again, no tuning opportunity here. */
void
FaaDiBruno::calculate(const UnfoldedStackContainer &cont, const UGSContainer &g,
                      UGSTensor &out)
{
  out.zeros();
  for (int l = 1; l <= out.dimen(); l++)
    {
      long int mem = SystemResources::availableMemory();
      cont.multAndAdd(l, g, out);
      JournalRecord rec(journal);
      int mem_mb = mem/1024/1024;
      rec << "dim=" << l << " avmem=" << mem_mb << endrec;
    }
}

/* This function returns a number of maximum rows used for refinement of the
   stacked container. We want to set the maximum so that the expected memory
   consumption for the number of paralel threads would be less than available
   memory. On the other hand we do not want to be too pesimistic since a very
   fine refinement can be very slow.

   Besides memory needed for a dense unfolded slice of a tensor from ‘f’, each
   thread needs ‘magic_mult*per_size’ bytes of memory. In the worst case,
   ‘magic_mult’ will be equal to 2, this means memory ‘per_size’ for target
   temporary (permuted symmetry) tensor plus one copy for intermediate result.
   However, this shows to be too pesimistic, so we set ‘magic_mult’ to 1.5. The
   memory for permuted symmetry temporary tensor ‘per_size’ is estimated as a
   weigthed average of unfolded memory of the ‘out’ tensor and unfolded memory
   of a symetric tensor with the largest coordinate size. Some experiments
   showed that the best combination of the two is to take 100% if the latter,
   so we set ‘lambda’ to zero.

   The ‘max’ number of rows in the refined ‘cont’ must be such that each
   slice fits to remaining memory. Number of columns of the slice are
   never greater maxˡ. (This is not true, since stacks corresponding to
   unit/zero matrices cannot be further refined). We get an equation:

    nthreads·maxˡ·8·r = mem − magic_mult·nthreads·per_size·8·r

   where ‘mem’ is available memory in bytes, ‘nthreads’ is a number of threads,
   r is a number of rows, and 8 is ‘sizeof(double)’.

   If the right hand side is less than zero, we set ‘max’ to 10, just to let it
   do something. */

std::tuple<int, int, int>
FaaDiBruno::estimRefinement(const TensorDimens &tdims, int nr, int l)
{
  int nthreads = sthread::detach_thread_group::max_parallel_threads;
  long per_size1 = tdims.calcUnfoldMaxOffset();
  long per_size2 = static_cast<long>(std::pow(tdims.getNVS().getMax(), l));
  double lambda = 0.0;
  long per_size = sizeof(double)*nr
    *static_cast<long>(lambda*per_size1+(1-lambda)*per_size2);
  long mem = SystemResources::availableMemory();
  int max = 0;
  double num_cols = static_cast<double>(mem-magic_mult*nthreads*per_size)
    /nthreads/sizeof(double)/nr;
  if (num_cols > 0)
    {
      double maxd = std::pow(num_cols, 1.0/l);
      max = static_cast<int>(std::floor(maxd));
    }
  if (max == 0)
    {
      max = 10;
      JournalRecord rec(journal);
      rec << "dim=" << l << " run out of memory, imposing max=" << max;
      if (nthreads > 1)
        rec << " (decrease number of threads)";
      rec << endrec;
    }
  int avmem_mb = mem/1024/1024;
  int tmpmem_mb = nthreads*per_size/1024/1024;
  return { max, avmem_mb, tmpmem_mb };
}
