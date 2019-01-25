/* $Header: /var/lib/cvs/dynare_cpp/sylv/cc/BlockDiagonal.h,v 1.1.1.1 2004/06/04 13:00:20 kamenik Exp $ */

/* Tag $Name:  $ */

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
  BlockDiagonal(int p, const BlockDiagonal &b);
  BlockDiagonal(const BlockDiagonal &b) = default;
  BlockDiagonal(const QuasiTriangular &t);
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
