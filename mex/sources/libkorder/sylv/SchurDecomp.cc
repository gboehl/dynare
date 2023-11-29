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

#include "SchurDecomp.hh"

#include <memory>

#include <dynlapack.h>

SchurDecomp::SchurDecomp(const SqSylvMatrix& m) : q(m.nrows())
{
  lapack_int rows = m.nrows();
  SqSylvMatrix auxt(m);
  lapack_int lda = auxt.getLD(), ldvs = q.getLD();
  lapack_int sdim;
  auto wr = std::make_unique<double[]>(rows);
  auto wi = std::make_unique<double[]>(rows);
  lapack_int lwork = 6 * rows;
  auto work = std::make_unique<double[]>(lwork);
  lapack_int info;
  dgees("V", "N", nullptr, &rows, auxt.base(), &lda, &sdim, wr.get(), wi.get(), q.base(), &ldvs,
        work.get(), &lwork, nullptr, &info);
  t_storage = std::make_unique<QuasiTriangular>(auxt.getData(), rows);
  t = t_storage.get();
}

SchurDecomp::SchurDecomp(const QuasiTriangular& tr) :
    q(tr.nrows()), t_storage {std::make_unique<QuasiTriangular>(tr)}, t {t_storage.get()}
{
  q.setUnit();
}

SchurDecomp::SchurDecomp(QuasiTriangular& tr) : q(tr.nrows()), t {&tr}
{
  q.setUnit();
}

int
SchurDecomp::getDim() const
{
  return t->nrows();
}

SchurDecompZero::SchurDecompZero(const GeneralMatrix& m) :
    SchurDecomp(SqSylvMatrix(m, m.nrows() - m.ncols(), 0, m.ncols())),
    ru(m, 0, 0, m.nrows() - m.ncols(), m.ncols())
{
  ru.multRight(getQ());
}

int
SchurDecompZero::getDim() const
{
  return getT().nrows() + ru.nrows();
}
