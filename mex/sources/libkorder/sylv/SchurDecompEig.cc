/*
 * Copyright © 2004-2011 Ondra Kamenik
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

#include "SchurDecompEig.hh"
#include "SylvException.hh"

#include <dynlapack.h>

#include <memory>

/* Bubble diagonal 1×1 or 2×2 block from position ‘from’ to position
   ‘to’. If an eigenvalue cannot be swapped with its neighbour, the
   neighbour is bubbled also in front. The method returns a new
   position ‘to’, where the original block pointed by ‘to’ happens to
   appear at the end. ‘from’ must be greater than ‘to’.
*/
SchurDecompEig::diag_iter
SchurDecompEig::bubbleEigen(diag_iter from, diag_iter to)
{
  diag_iter run = from;
  while (run != to)
    {
      diag_iter runm = run;
      if (!tryToSwap(run, runm) && runm == to)
        ++to;
      else
        {
          /* Bubble all eigenvalues from runm(incl.) to run(excl.),
             this includes either bubbling generated eigenvalues due
             to split, or an eigenvalue which couldn't be swapped */
          while (runm != run)
            {
              to = bubbleEigen(runm, to);
              ++runm;
            }
        }
    }
  return to;
}

/* This tries to swap two neighbouring eigenvalues, ‘it’ and ‘--it’,
   and returns ‘itadd’. If the blocks can be swapped, new eigenvalues
   can emerge due to possible 2×2 block splits. ‘it’ then points to
   the last eigenvalue coming from block pointed by ‘it’ at the
   begining, and ‘itadd’ points to the first. On swap failure, ‘it’ is
   not changed, and ‘itadd’ points to previous eignevalue (which must
   be moved backwards before). In either case, it is necessary to
   resolve eigenvalues from ‘itadd’ to ‘it’, before the ‘it’ can be
   resolved.
   The success is signaled by returned true.
*/
bool
SchurDecompEig::tryToSwap(diag_iter &it, diag_iter &itadd)
{
  itadd = it;
  --itadd;

  lapack_int n = getDim(), ldt = getT().getLD(), ldq = getQ().getLD();
  lapack_int ifst = it->getIndex() + 1;
  lapack_int ilst = itadd->getIndex() + 1;
  auto work = std::make_unique<double[]>(n);
  lapack_int info;
  dtrexc("V", &n, getT().base(), &ldt, getQ().base(), &ldq, &ifst, &ilst, work.get(),
         &info);
  if (info < 0)
    throw SYLV_MES_EXCEPTION("Wrong argument to dtrexc.");

  if (info == 0)
    {
      // swap successful
      getT().swapDiagLogically(itadd);
      // check for 2×2 block splits
      getT().checkDiagConsistency(it);
      getT().checkDiagConsistency(itadd);
      // and go back by ‘it’ in NEW eigenvalue set
      --it;
      return true;
    }
  return false;
}

void
SchurDecompEig::orderEigen()
{
  diag_iter run = getT().diag_begin();
  diag_iter runp = run;
  ++runp;
  double last_size = 0.0;
  while (runp != getT().diag_end())
    {
      diag_iter least = getT().findNextLargerBlock(run, getT().diag_end(),
                                                   last_size);
      last_size = least->getSize();
      if (run == least)
        ++run;
      else
        run = bubbleEigen(least, run);
      runp = run;
      ++runp;
    }
}
