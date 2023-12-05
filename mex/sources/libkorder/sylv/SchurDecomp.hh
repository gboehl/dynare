/*
 * Copyright © 2004-2011 Ondra Kamenik
 * Copyright © 2019-2023 Dynare Team
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

#ifndef SCHUR_DECOMP_H
#define SCHUR_DECOMP_H

#include "QuasiTriangular.hh"
#include "SylvMatrix.hh"

#include <memory>

class QuasiTriangular;
class SchurDecomp
{
  SqSylvMatrix q;
  // Stores t if is owned
  std::unique_ptr<QuasiTriangular> t_storage;
  QuasiTriangular* t;

public:
  SchurDecomp(const SqSylvMatrix& m);
  SchurDecomp(const QuasiTriangular& tr);
  SchurDecomp(QuasiTriangular& tr);
  [[nodiscard]] const SqSylvMatrix&
  getQ() const
  {
    return q;
  }
  [[nodiscard]] const QuasiTriangular&
  getT() const
  {
    return *t;
  }
  SqSylvMatrix&
  getQ()
  {
    return q;
  }
  QuasiTriangular&
  getT()
  {
    return *t;
  }
  [[nodiscard]] virtual int getDim() const;
  virtual ~SchurDecomp() = default;
};

class SchurDecompZero : public SchurDecomp
{
  GeneralMatrix ru; // right upper matrix
public:
  SchurDecompZero(const GeneralMatrix& m);
  [[nodiscard]] ConstGeneralMatrix
  getRU() const
  {
    return ru;
  }
  [[nodiscard]] int getDim() const override;
  [[nodiscard]] int
  getZeroCols() const
  {
    return ru.nrows();
  }
};

#endif /* SCHUR_DECOMP_H */
