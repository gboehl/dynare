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

#ifndef BLOCK_DIAGONAL_H
#define BLOCK_DIAGONAL_H

#include <memory>
#include <vector>

#include "QuasiTriangular.hh"

class BlockDiagonal : public QuasiTriangular
{
  std::vector<int> row_len, col_len;
public:
  BlockDiagonal(ConstVector d, int d_size);
  BlockDiagonal(const BlockDiagonal &b) = default;
  explicit BlockDiagonal(const QuasiTriangular &t);
  BlockDiagonal &operator=(const QuasiTriangular &t)
  {
    GeneralMatrix::operator=(t);
    return *this;
  }
  BlockDiagonal &operator=(const BlockDiagonal &b) = default;
  ~BlockDiagonal() override = default;
  void setZeroBlockEdge(diag_iter edge);
  int getNumZeros() const;
  int getNumBlocks() const;
  int getLargestBlock() const;
  void printInfo() const;

  void multKron(KronVector &x) const override;
  void multKronTrans(KronVector &x) const override;

  const_col_iter col_begin(const DiagonalBlock &b) const override;
  col_iter col_begin(const DiagonalBlock &b) override;
  const_row_iter row_end(const DiagonalBlock &b) const override;
  row_iter row_end(const DiagonalBlock &b) override;
  std::unique_ptr<QuasiTriangular>
  clone() const override
  {
    return std::make_unique<BlockDiagonal>(*this);
  }
private:
  void setZerosToRU(diag_iter edge);
  const_diag_iter findBlockStart(const_diag_iter from) const;
  static void savePartOfX(int si, int ei, const KronVector &x, Vector &work);
  void multKronBlock(const_diag_iter start, const_diag_iter end,
                     KronVector &x, Vector &work) const;
  void multKronBlockTrans(const_diag_iter start, const_diag_iter end,
                          KronVector &x, Vector &work) const;
};

#endif /* BLOCK_DIAGONAL_H */
