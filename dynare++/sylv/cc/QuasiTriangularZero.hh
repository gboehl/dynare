/* $Header: /var/lib/cvs/dynare_cpp/sylv/cc/QuasiTriangularZero.h,v 1.1.1.1 2004/06/04 13:00:44 kamenik Exp $ */

/* Tag $Name:  $ */

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
