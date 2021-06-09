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

#include "QuasiTriangularZero.hh"
#include "SchurDecomp.hh"
#include "SylvMatrix.hh"
#include "SylvException.hh"

#include <iostream>

QuasiTriangularZero::QuasiTriangularZero(int num_zeros, const ConstVector &d,
                                         int d_size)
  : QuasiTriangular(SqSylvMatrix(GeneralMatrix(Vector{d}, num_zeros+d_size, d_size),
                                 num_zeros, 0, d_size).getData(),
                    d_size),
    nz(num_zeros),
    ru(GeneralMatrix(Vector{d}, num_zeros+d_size, d_size), 0, 0, num_zeros, d_size)
{
}

QuasiTriangularZero::QuasiTriangularZero(double r,
                                         const QuasiTriangularZero &t)
  : QuasiTriangular(r, t),
    nz(t.nz),
    ru(t.ru)
{
  ru.mult(r);
}

QuasiTriangularZero::QuasiTriangularZero(double r,
                                         const QuasiTriangularZero &t,
                                         double r2,
                                         const QuasiTriangularZero &t2)
  : QuasiTriangular(r, t, r2, t2),
    nz(t.nz),
    ru(t.ru)
{
  ru.mult(r);
  ru.add(r2, t2.ru);
}

QuasiTriangularZero::QuasiTriangularZero(const std::string &dummy, const QuasiTriangularZero &t)
  : QuasiTriangular(dummy, t),
    nz(t.nz),
    ru(t.ru)
{
  ru.multRight(t);
}

QuasiTriangularZero::QuasiTriangularZero(const SchurDecompZero &decomp)
  : QuasiTriangular(decomp.getT().getData(),
                    decomp.getT().nrows()),
    nz(decomp.getZeroCols()),
    ru(decomp.getRU())
{
}

QuasiTriangularZero::QuasiTriangularZero(const QuasiTriangular &t)
  : QuasiTriangular(t),
    nz(0), ru(0, t.getDiagonal().getSize())
{
}

void
QuasiTriangularZero::solvePre(Vector &x, double &eig_min)
{
  Vector xu(x, 0, nz);
  Vector xl(x, nz, x.length()-nz);
  QuasiTriangular::solvePre(xl, eig_min);
  ru.multsVec(xu, xl);
  if (nz > 0)
    eig_min = (eig_min > 1.0) ? 1.0 : eig_min;
}

void
QuasiTriangularZero::solvePreTrans(Vector &x, double &eig_min)
{
  Vector xu(x, 0, nz);
  Vector xl(x, nz, x.length()-nz);
  ru.multsVecTrans(xl, xu);
  QuasiTriangular::solvePreTrans(xl, eig_min);
  if (nz > 0)
    eig_min = (eig_min > 1.0) ? 1.0 : eig_min;
}

void
QuasiTriangularZero::multVec(Vector &x, const ConstVector &b) const
{
  x.zeros();
  multaVec(x, b);
}

void
QuasiTriangularZero::multVecTrans(Vector &x, const ConstVector &b) const
{
  x.zeros();
  multaVecTrans(x, b);
}

void
QuasiTriangularZero::multaVec(Vector &x, const ConstVector &b) const
{
  ConstVector bl(b, nz, b.length()-nz);
  Vector xu(x, 0, nz);
  Vector xl(x, nz, x.length()-nz);
  xu.zeros();
  ru.multaVec(xu, bl);
  QuasiTriangular::multVec(xl, bl);
}

void
QuasiTriangularZero::multaVecTrans(Vector &x, const ConstVector &b) const
{
  ConstVector bu(b, 0, b.length());
  ConstVector bl(b, nz, b.length()-nz);
  Vector xu(x, 0, nz);
  Vector xl(x, nz, x.length()-nz);
  xu.zeros();
  QuasiTriangular::multVecTrans(xl, bl);
  ru.multaVecTrans(xl, bu);
}

void
QuasiTriangularZero::multLeftOther(GeneralMatrix &a) const
{
  GeneralMatrix a1(a, 0, 0, nz, a.ncols());
  GeneralMatrix a2(a, nz, 0, a.nrows()-nz, a.ncols());
  a1.mult(ru, a2);
  QuasiTriangular::multLeftOther(a2);
}

void
QuasiTriangularZero::print() const
{
  std::cout << "super=" << std::endl;
  QuasiTriangular::print();
  std::cout << "nz=" << nz << std::endl
            << "ru=" << std::endl;
  ru.print();
}

void
QuasiTriangularZero::multKron(KronVector &x) const
{
  throw SYLV_MES_EXCEPTION("Attempt to run QuasiTriangularZero::multKron.");
}

void
QuasiTriangularZero::multKronTrans(KronVector &x) const
{
  throw SYLV_MES_EXCEPTION("Attempt to run QuasiTriangularZero::multKronTrans.");
}
