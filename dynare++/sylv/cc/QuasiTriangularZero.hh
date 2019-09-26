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
 * along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef QUASI_TRIANGULAR_ZERO_H
#define QUASI_TRIANGULAR_ZERO_H

#include "QuasiTriangular.hh"
#include "GeneralMatrix.hh"

#include <memory>

/*
   Stores a (square) quasi-triangular matrix whose first columns are zero:
    ⎛0 R⎞
    ⎝0 M⎠
   where M (square quasi-triangular) is stored in the super-class, and R (rectangular)
   is stored in ‘ru’.
*/

class QuasiTriangularZero : public QuasiTriangular
{
  int nz; // number of zero columns
  GeneralMatrix ru; // data in right upper part (of size nz×d_size)
public:
  QuasiTriangularZero(int num_zeros, const ConstVector &d, int d_size);
  // Initializes with r·t
  QuasiTriangularZero(double r, const QuasiTriangularZero &t);
  // Initializes with r·t+r₂·t₂
  QuasiTriangularZero(double r, const QuasiTriangularZero &t,
                      double r2, const QuasiTriangularZero &t2);
  // Initializes with t²
  QuasiTriangularZero(const std::string &dummy, const QuasiTriangularZero &t);
  explicit QuasiTriangularZero(const QuasiTriangular &t);
  explicit QuasiTriangularZero(const SchurDecompZero &decomp);
  ~QuasiTriangularZero() override = default;
  void solvePre(Vector &x, double &eig_min) override;
  void solvePreTrans(Vector &x, double &eig_min) override;
  void multVec(Vector &x, const ConstVector &b) const override;
  void multVecTrans(Vector &x, const ConstVector &b) const override;
  void multaVec(Vector &x, const ConstVector &b) const override;
  void multaVecTrans(Vector &x, const ConstVector &b) const override;
  void multKron(KronVector &x) const override;
  void multKronTrans(KronVector &x) const override;
  void multLeftOther(GeneralMatrix &a) const override;

  std::unique_ptr<QuasiTriangular>
  clone() const override
  {
    return std::make_unique<QuasiTriangularZero>(*this);
  }
  std::unique_ptr<QuasiTriangular>
  square() const override
  {
    return std::make_unique<QuasiTriangularZero>("square", *this);
  }
  std::unique_ptr<QuasiTriangular>
  scale(double r) const override
  {
    return std::make_unique<QuasiTriangularZero>(r, *this);
  }
  std::unique_ptr<QuasiTriangular>
  linearlyCombine(double r, double r2, const QuasiTriangular &t2) const override
  {
    return std::make_unique<QuasiTriangularZero>(r, *this, r2, dynamic_cast<const QuasiTriangularZero &>(t2));
  }
  void print() const override;
};

#endif /* QUASI_TRIANGULAR_ZERO_H */
