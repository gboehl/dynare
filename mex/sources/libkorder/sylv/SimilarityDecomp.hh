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

#ifndef SIMILARITY_DECOMP_H
#define SIMILARITY_DECOMP_H

#include "BlockDiagonal.hh"
#include "SylvMatrix.hh"
#include "SylvParams.hh"

#include <memory>

class SimilarityDecomp
{
  std::unique_ptr<SqSylvMatrix> q;
  std::unique_ptr<BlockDiagonal> b;
  std::unique_ptr<SqSylvMatrix> invq;
  using diag_iter = BlockDiagonal::diag_iter;

public:
  SimilarityDecomp(const ConstVector& d, int d_size, double log10norm = 3.0);
  virtual ~SimilarityDecomp() = default;
  [[nodiscard]] const SqSylvMatrix&
  getQ() const
  {
    return *q;
  }
  [[nodiscard]] const SqSylvMatrix&
  getInvQ() const
  {
    return *invq;
  }
  [[nodiscard]] const BlockDiagonal&
  getB() const
  {
    return *b;
  }
  void check(SylvParams& pars, const GeneralMatrix& m) const;
  void infoToPars(SylvParams& pars) const;

protected:
  void getXDim(diag_iter start, diag_iter end, int& rows, int& cols) const;
  bool solveX(diag_iter start, diag_iter end, GeneralMatrix& X, double norm) const;
  void updateTransform(diag_iter start, diag_iter end, GeneralMatrix& X);
  void bringGuiltyBlock(diag_iter start, diag_iter& end);
  void diagonalize(double norm);
};

#endif /* SIMILARITY_DECOMP_H */
