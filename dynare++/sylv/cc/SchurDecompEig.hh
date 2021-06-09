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

// contains algorithms for eigenvalue reordering

#ifndef SCHUR_DECOMP_EIG_H
#define SCHUR_DECOMP_EIG_H

#include "SchurDecomp.hh"
#include "QuasiTriangular.hh"

class SchurDecompEig : public SchurDecomp
{
public:
  using diag_iter = QuasiTriangular::diag_iter;
  SchurDecompEig(const SqSylvMatrix &m) : SchurDecomp(m)
  {
  }
  SchurDecompEig(const QuasiTriangular &tr) : SchurDecomp(tr)
  {
  };
  SchurDecompEig(QuasiTriangular &tr) : SchurDecomp(tr)
  {
  }
  diag_iter bubbleEigen(diag_iter from, diag_iter to);
  void orderEigen();
protected:
  bool tryToSwap(diag_iter &it, diag_iter &itadd);
};

#endif /* SCHUR_DECOMP_EIG_H */
