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

#include "SimilarityDecomp.hh"
#include "SchurDecomp.hh"
#include "SchurDecompEig.hh"
#include "SylvException.hh"

#include <dynlapack.h>

#include <cmath>

SimilarityDecomp::SimilarityDecomp(const ConstVector &d, int d_size, double log10norm)
{
  SchurDecomp sd(SqSylvMatrix(Vector{d}, d_size));
  q = std::make_unique<SqSylvMatrix>(sd.getQ());
  b = std::make_unique<BlockDiagonal>(sd.getT());
  invq = std::make_unique<SqSylvMatrix>(d_size);
  invq->setUnit();
  invq->multLeftTrans(sd.getQ());
  double norm = pow(10.0, log10norm);
  diagonalize(norm);
}

void
SimilarityDecomp::getXDim(diag_iter start, diag_iter end,
                          int &rows, int &cols) const
{
  int si = start->getIndex();
  int ei = end->getIndex();
  cols = b->nrows() - ei;
  rows = ei - si;
}

/* Find solution of X for diagonal block given by start(incl.) and
   end(excl.). If the solution cannot be found, or it is greater than
   norm, X is not changed and flase is returned.
*/
bool
SimilarityDecomp::solveX(diag_iter start, diag_iter end,
                         GeneralMatrix &X, double norm) const
{
  int si = start->getIndex();
  int ei = end->getIndex();

  SqSylvMatrix A(const_cast<const BlockDiagonal &>(*b), si, si, X.nrows());
  SqSylvMatrix B(const_cast<const BlockDiagonal &>(*b), ei, ei, X.ncols());
  GeneralMatrix C(const_cast<const BlockDiagonal &>(*b), si, ei, X.nrows(), X.ncols());

  lapack_int isgn = -1;
  lapack_int m = A.nrows();
  lapack_int n = B.nrows();
  lapack_int lda = A.getLD(), ldb = B.getLD();
  double scale;
  lapack_int info;
  dtrsyl("N", "N", &isgn, &m, &n, A.base(), &lda, B.base(), &ldb,
         C.base(), &m, &scale, &info);
  if (info < -1)
    throw SYLV_MES_EXCEPTION("Wrong parameter to LAPACK dtrsyl.");

  if (info == 1 || scale < 1)
    return false;
  if (C.getData().getMax() > norm)
    return false;

  X = C;
  return true;
}

/*                         ⎛I −X⎞     ⎛I X⎞
  Multiply Q and invQ with ⎝0  I⎠ and ⎝0 I⎠ respectively. This also sets X=−X. */
void
SimilarityDecomp::updateTransform(diag_iter start, diag_iter end,
                                  GeneralMatrix &X)
{
  int si = start->getIndex();
  int ei = end->getIndex();

  SqSylvMatrix iX(q->nrows());
  iX.setUnit();
  iX.place(X, si, ei);
  invq->GeneralMatrix::multLeft(iX);

  iX.setUnit();
  X.mult(-1.0);
  iX.place(X, si, ei);
  q->multRight(iX);
}

void
SimilarityDecomp::bringGuiltyBlock(diag_iter start, diag_iter &end)
{
  double av = b->getAverageDiagSize(start, end);
  diag_iter guilty = b->findClosestDiagBlock(end, b->diag_end(), av);
  SchurDecompEig sd(*b); // works on b including diagonal structure
  end = sd.bubbleEigen(guilty, end); // iterators are valid
  ++end;
  q->multRight(sd.getQ());
  invq->multLeftTrans(sd.getQ());
}

void
SimilarityDecomp::diagonalize(double norm)
{
  diag_iter start = b->diag_begin();
  diag_iter end = start;
  ++end;

  while (end != b->diag_end())
    {
      int xrows;
      int xcols;
      getXDim(start, end, xrows, xcols);
      GeneralMatrix X(xrows, xcols);
      if (solveX(start, end, X, norm))
        {
          updateTransform(start, end, X);
          b->setZeroBlockEdge(end);
          start = end;
          ++end;
        }
      else
        bringGuiltyBlock(start, end); // moves with end
    }
}

void
SimilarityDecomp::check(SylvParams &pars, const GeneralMatrix &m) const
{
  // M − Q·B·Q⁻¹
  SqSylvMatrix c(getQ() * getB());
  c.multRight(getInvQ());
  c.add(-1.0, m);
  pars.f_err1 = c.getNorm1();
  pars.f_errI = c.getNormInf();

  // I − Q·Q⁻¹
  c.setUnit();
  c.mult(-1);
  c.multAndAdd(getQ(), getInvQ());
  pars.viv_err1 = c.getNorm1();
  pars.viv_errI = c.getNormInf();

  // I − Q⁻¹·Q
  c.setUnit();
  c.mult(-1);
  c.multAndAdd(getInvQ(), getQ());
  pars.ivv_err1 = c.getNorm1();
  pars.ivv_errI = c.getNormInf();
}

void
SimilarityDecomp::infoToPars(SylvParams &pars) const
{
  pars.f_blocks = getB().getNumBlocks();
  pars.f_largest = getB().getLargestBlock();
  pars.f_zeros = getB().getNumZeros();
  pars.f_offdiag = getB().getNumOffdiagonal();
}
