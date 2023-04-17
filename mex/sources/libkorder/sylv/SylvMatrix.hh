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

#ifndef SYLV_MATRIX_H
#define SYLV_MATRIX_H

#include "GeneralMatrix.hh"
#include "KronVector.hh"

#include <utility>

class SqSylvMatrix;

class SylvMatrix : public GeneralMatrix
{
public:
  SylvMatrix(int m, int n)
    : GeneralMatrix(m, n)
  {
  }
  SylvMatrix(Vector d, int m, int n)
    : GeneralMatrix(std::move(d), m, n)
  {
  }
  SylvMatrix(const GeneralMatrix &m)
    : GeneralMatrix(m)
  {
  }
  SylvMatrix(const GeneralMatrix &m, int i, int j, int nrows, int ncols)
    : GeneralMatrix(m, i, j, nrows, ncols)
  {
  }
  SylvMatrix(GeneralMatrix &m, int i, int j, int nrows, int ncols)
    : GeneralMatrix(m, i, j, nrows, ncols)
  {
  }
  SylvMatrix(const SylvMatrix &m) = default;
  SylvMatrix(SylvMatrix &&m) = default;
  SylvMatrix &operator=(const SylvMatrix &m) = default;
  SylvMatrix &operator=(SylvMatrix &&m) = default;

  /*        ⎛I 0⎞
     this = ⎝0 m⎠·this */
  void multLeftI(const SqSylvMatrix &m);
  /*        ⎛I  0⎞
     this = ⎝0 mᵀ⎠·this */
  void multLeftITrans(const SqSylvMatrix &m);
  // this = (0 a)·b, so that (0 a) is square
  void multLeft(int zero_cols, const GeneralMatrix &a, const GeneralMatrix &b);
  // this = this·(m⊗m⊗…⊗m)
  void multRightKron(const SqSylvMatrix &m, int order);
  // this = this·(mᵀ⊗mᵀ⊗…⊗mᵀ)
  void multRightKronTrans(const SqSylvMatrix &m, int order);
  /* this = P·this, x = P·x, where P is gauss transformation setting
     a given element to zero */
  void eliminateLeft(int row, int col, Vector &x);
  /* this = this·P, x = Pᵀ·x, where P is gauss transformation setting
     a given element to zero */
  void eliminateRight(int row, int col, Vector &x);
};

class SqSylvMatrix : public SylvMatrix
{
public:
  SqSylvMatrix(int m) : SylvMatrix(m, m)
  {
  }
  SqSylvMatrix(Vector d, int m) : SylvMatrix(std::move(d), m, m)
  {
  }
  SqSylvMatrix(const SylvMatrix &m) : SylvMatrix(m)
  {
  }
  SqSylvMatrix(const SqSylvMatrix &m) = default;
  SqSylvMatrix(SqSylvMatrix &&m) = default;
  SqSylvMatrix(const GeneralMatrix &m, int i, int j, int nrows)
    : SylvMatrix(m, i, j, nrows, nrows)
  {
  }
  SqSylvMatrix(GeneralMatrix &m, int i, int j, int nrows)
    : SylvMatrix(m, i, j, nrows, nrows)
  {
  }
  SqSylvMatrix &operator=(const SqSylvMatrix &m) = default;
  SqSylvMatrix &operator=(SqSylvMatrix &&m) = default;
  // x = (this⊗this⊗…⊗this)·d
  void multVecKron(KronVector &x, const ConstKronVector &d) const;
  // x = (thisᵀ⊗thisᵀ⊗…⊗thisᵀ)·d
  void multVecKronTrans(KronVector &x, const ConstKronVector &d) const;
  // a = this⁻¹·a, b=this⁻¹·b */
  void multInvLeft2(GeneralMatrix &a, GeneralMatrix &b,
                    double &rcond1, double &rcondinf) const;
  // this = I
  void setUnit();
};

#endif /* SYLV_MATRIX_H */
