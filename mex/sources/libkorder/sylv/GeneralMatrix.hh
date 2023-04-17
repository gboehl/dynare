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

#ifndef GENERAL_MATRIX_H
#define GENERAL_MATRIX_H

#include "Vector.hh"
#include "SylvException.hh"

#include <algorithm>
#include <memory>
#include <utility>
#include <string>

template<class T>
class TransposedMatrix
{
  friend class GeneralMatrix;
  template<class T2>
  friend GeneralMatrix operator*(const ConstGeneralMatrix &a, const TransposedMatrix<T2> &b);
  template<class T2>
  friend GeneralMatrix operator*(const TransposedMatrix<T2> &a, const ConstGeneralMatrix &b);
  template<class T1, class T2>
  friend GeneralMatrix operator*(const TransposedMatrix<T1> &a, const TransposedMatrix<T2> &b);
private:
  T &orig;
public:
  TransposedMatrix(T &orig_arg) : orig{orig_arg}
  {
  };
};

// Syntactic sugar for representing a transposed matrix
template<class T>
TransposedMatrix<T>
transpose(T &m)
{
  return TransposedMatrix<T>(m);
}

class GeneralMatrix;

class ConstGeneralMatrix
{
  friend class GeneralMatrix;
  friend GeneralMatrix operator*(const ConstGeneralMatrix &a, const ConstGeneralMatrix &b);
  template<class T>
  friend GeneralMatrix operator*(const ConstGeneralMatrix &a, const TransposedMatrix<T> &b);
  template<class T>
  friend GeneralMatrix operator*(const TransposedMatrix<T> &a, const ConstGeneralMatrix &b);
  template<class T1, class T2>
  friend GeneralMatrix operator*(const TransposedMatrix<T1> &a, const TransposedMatrix<T2> &b);
protected:
  ConstVector data; // Has unit-stride
  int rows;
  int cols;
  int ld;
public:
  ConstGeneralMatrix(ConstVector d, int m, int n)
    : data(std::move(d)), rows(m), cols(n), ld(m)
  {
    if (data.skip() > 1)
      throw SYLV_MES_EXCEPTION("Vector must have unit-stride");
    if (data.length() < m*n)
      throw SYLV_MES_EXCEPTION("Vector is too small");
  }
  ConstGeneralMatrix(const ConstGeneralMatrix &m) = default;
  ConstGeneralMatrix(ConstGeneralMatrix &&m) = default;
  // Implicit conversion from GeneralMatrix is ok, since it's cheap
  ConstGeneralMatrix(const GeneralMatrix &m);
  // Create submatrix (with data sharing)
  ConstGeneralMatrix(const GeneralMatrix &m, int i, int j, int nrows, int ncols);
  // Create submatrix (with data sharing)
  ConstGeneralMatrix(const ConstGeneralMatrix &m, int i, int j, int nrows, int ncols);
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
  explicit ConstGeneralMatrix(const mxArray *p)
    : data(p), rows{static_cast<int>(mxGetM(p))}, cols{static_cast<int>(mxGetN(p))}, ld{rows}
  {
  }
#endif
  virtual ~ConstGeneralMatrix() = default;

  ConstGeneralMatrix &operator=(const ConstGeneralMatrix &v) = delete;
  ConstGeneralMatrix &operator=(ConstGeneralMatrix &&v) = delete;

  const double &
  get(int i, int j) const
  {
    return data[j*ld+i];
  }
  int
  nrows() const
  {
    return rows;
  }
  int
  ncols() const
  {
    return cols;
  }
  int
  getLD() const
  {
    return ld;
  }
  const double *
  base() const
  {
    return data.base();
  }
  const ConstVector &
  getData() const
  {
    return data;
  }
  ConstVector getRow(int row) const;
  ConstVector getCol(int col) const;

  double getNormInf() const;
  double getNorm1() const;
  /* x = scalar(a)*x + scalar(b)*this*d */
  void multVec(double a, Vector &x, double b, const ConstVector &d) const;
  /* x = scalar(a)*x + scalar(b)*this'*d */
  void multVecTrans(double a, Vector &x, double b, const ConstVector &d) const;
  /* x = x + this*d */
  void
  multaVec(Vector &x, const ConstVector &d) const
  {
    multVec(1.0, x, 1.0, d);
  }
  /* x = x + this'*d */
  void
  multaVecTrans(Vector &x, const ConstVector &d) const
  {
    multVecTrans(1.0, x, 1.0, d);
  }
  /* x = x - this*d */
  void
  multsVec(Vector &x, const ConstVector &d) const
  {
    multVec(1.0, x, -1.0, d);
  }
  /* x = x - this'*d */
  void
  multsVecTrans(Vector &x, const ConstVector &d) const
  {
    multVecTrans(1.0, x, -1.0, d);
  }
  /* m = inv(this)*m */
  void multInvLeft(GeneralMatrix &m) const;
  /* m = inv(this')*m */
  void multInvLeftTrans(GeneralMatrix &m) const;
  /* d = inv(this)*d */
  void multInvLeft(Vector &d) const;
  /* d = inv(this')*d */
  void multInvLeftTrans(Vector &d) const;

  bool isFinite() const;
  /** Returns true of the matrix is exactly zero. */
  bool isZero() const;

  virtual void print() const;
protected:
  void multInvLeft(const std::string &trans, int mrows, int mcols, int mld, double *d) const;
};

class GeneralMatrix
{
  friend class ConstGeneralMatrix;
  friend GeneralMatrix operator*(const ConstGeneralMatrix &a, const ConstGeneralMatrix &b);
  template<class T>
  friend GeneralMatrix operator*(const ConstGeneralMatrix &a, const TransposedMatrix<T> &b);
  template<class T>
  friend GeneralMatrix operator*(const TransposedMatrix<T> &a, const ConstGeneralMatrix &b);
  template<class T1, class T2>
  friend GeneralMatrix operator*(const TransposedMatrix<T1> &a, const TransposedMatrix<T2> &b);
protected:
  Vector data; // Has unit-stride
  int rows;
  int cols;
  int ld;
public:
  GeneralMatrix(int m, int n)
    : data(m*n), rows(m), cols(n), ld(m)
  {
  }
  GeneralMatrix(Vector d, int m, int n)
    : data(std::move(d)), rows(m), cols(n), ld(m)
  {
    if (data.skip() > 1)
      throw SYLV_MES_EXCEPTION("Vector must have unit-stride");
    if (data.length() < m*n)
      throw SYLV_MES_EXCEPTION("Vector is too small");
  }

  /* The copies will have ld==rows, for memory efficiency (hence we do not use
     the default copy constructor) */
  GeneralMatrix(const GeneralMatrix &m);
  // We don't want implict conversion from ConstGeneralMatrix, since it's expensive
  explicit GeneralMatrix(const ConstGeneralMatrix &m);

  GeneralMatrix(GeneralMatrix &&m) = default;

  template<class T>
  explicit GeneralMatrix(const TransposedMatrix<T> &m)
    : data(m.orig.rows*m.orig.cols), rows(m.orig.cols), cols(m.orig.rows), ld(rows)
  {
    for (int i = 0; i < rows; i++)
      for (int j = 0; j < cols; j++)
        get(i, j) = m.orig.get(j, i);
  }

  // Create submatrix (with data copy)
  GeneralMatrix(const GeneralMatrix &m, int i, int j, int nrows, int ncols);
  // Create submatrix (with data sharing)
  GeneralMatrix(GeneralMatrix &m, int i, int j, int nrows, int ncols);

#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
  explicit GeneralMatrix(mxArray *p)
    : data(p), rows{static_cast<int>(mxGetM(p))}, cols{static_cast<int>(mxGetN(p))}, ld{rows}
  {
  }
#endif

  virtual ~GeneralMatrix() = default;
  GeneralMatrix &operator=(const GeneralMatrix &m) = default;
  GeneralMatrix &operator=(GeneralMatrix &&m) = default;
  GeneralMatrix &operator=(const ConstGeneralMatrix &m);

  const double &
  get(int i, int j) const
  {
    return data[j*ld+i];
  }
  double &
  get(int i, int j)
  {
    return data[j*ld+i];
  }
  int
  nrows() const
  {
    return rows;
  }
  int
  ncols() const
  {
    return cols;
  }
  int
  getLD() const
  {
    return ld;
  }
  double *
  base()
  {
    return data.base();
  }
  const double *
  base() const
  {
    return data.base();
  }
  Vector &
  getData()
  {
    return data;
  }
  ConstVector
  getData() const
  {
    return data;
  }
  Vector getRow(int row);
  Vector getCol(int col);
  ConstVector getRow(int row) const;
  ConstVector getCol(int col) const;

  double
  getNormInf() const
  {
    return ConstGeneralMatrix(*this).getNormInf();
  }
  double
  getNorm1() const
  {
    return ConstGeneralMatrix(*this).getNorm1();
  }

  /* place matrix m to the position (i,j) */
  void place(const ConstGeneralMatrix &m, int i, int j);
  void
  place(const GeneralMatrix &m, int i, int j)
  {
    place(ConstGeneralMatrix(m), i, j);
  }

  // this = a·b
  void mult(const ConstGeneralMatrix &a, const ConstGeneralMatrix &b);
  void
  mult(const GeneralMatrix &a, const GeneralMatrix &b)
  {
    mult(ConstGeneralMatrix(a), ConstGeneralMatrix(b));
  }

  // this = this + scalar·a·b
  void multAndAdd(const ConstGeneralMatrix &a, const ConstGeneralMatrix &b,
                  double mult = 1.0);
  void
  multAndAdd(const GeneralMatrix &a, const GeneralMatrix &b,
             double mult = 1.0)
  {
    multAndAdd(ConstGeneralMatrix(a), ConstGeneralMatrix(b), mult);
  }

  // this = this + scalar·a·bᵀ
  void multAndAdd(const ConstGeneralMatrix &a, const ConstGeneralMatrix &b,
                  const std::string &dum, double mult = 1.0);
  void
  multAndAdd(const GeneralMatrix &a, const GeneralMatrix &b,
             const std::string &dum, double mult = 1.0)
  {
    multAndAdd(ConstGeneralMatrix(a), ConstGeneralMatrix(b), dum, mult);
  }

  // this = this + scalar·aᵀ·b
  void multAndAdd(const ConstGeneralMatrix &a, const std::string &dum, const ConstGeneralMatrix &b,
                  double mult = 1.0);
  void
  multAndAdd(const GeneralMatrix &a, const std::string &dum, const GeneralMatrix &b,
             double mult = 1.0)
  {
    multAndAdd(ConstGeneralMatrix(a), dum, ConstGeneralMatrix(b), mult);
  }

  // this = this + scalar·aᵀ·bᵀ
  void multAndAdd(const ConstGeneralMatrix &a, const std::string &dum1,
                  const ConstGeneralMatrix &b, const std::string &dum2, double mult = 1.0);
  void
  multAndAdd(const GeneralMatrix &a, const std::string &dum1,
             const GeneralMatrix &b, const std::string &dum2, double mult = 1.0)
  {
    multAndAdd(ConstGeneralMatrix(a), dum1, ConstGeneralMatrix(b), dum2, mult);
  }

  // this = this + scalar·a·aᵀ
  void addOuter(const ConstVector &a, double mult = 1.0);
  void
  addOuter(const Vector &a, double mult = 1.0)
  {
    addOuter(ConstVector(a), mult);
  }

  // this = this·m
  void multRight(const ConstGeneralMatrix &m);
  void
  multRight(const GeneralMatrix &m)
  {
    multRight(ConstGeneralMatrix(m));
  }

  // this = m·this
  void multLeft(const ConstGeneralMatrix &m);
  void
  multLeft(const GeneralMatrix &m)
  {
    multLeft(ConstGeneralMatrix(m));
  }

  // this = this·mᵀ
  void multRightTrans(const ConstGeneralMatrix &m);
  void
  multRightTrans(const GeneralMatrix &m)
  {
    multRightTrans(ConstGeneralMatrix(m));
  }

  // this = mᵀ·this
  void multLeftTrans(const ConstGeneralMatrix &m);
  void
  multLeftTrans(const GeneralMatrix &m)
  {
    multLeftTrans(ConstGeneralMatrix(m));
  }

  // x = scalar(a)·x + scalar(b)·this·d
  void
  multVec(double a, Vector &x, double b, const ConstVector &d) const
  {
    ConstGeneralMatrix(*this).multVec(a, x, b, d);
  }

  // x = scalar(a)·x + scalar(b)·thisᵀ·d
  void
  multVecTrans(double a, Vector &x, double b, const ConstVector &d) const
  {
    ConstGeneralMatrix(*this).multVecTrans(a, x, b, d);
  }

  // x = x + this·d
  void
  multaVec(Vector &x, const ConstVector &d) const
  {
    ConstGeneralMatrix(*this).multaVec(x, d);
  }

  // x = x + thisᵀ·d */
  void
  multaVecTrans(Vector &x, const ConstVector &d) const
  {
    ConstGeneralMatrix(*this).multaVecTrans(x, d);
  }

  // x = x - this·d
  void
  multsVec(Vector &x, const ConstVector &d) const
  {
    ConstGeneralMatrix(*this).multsVec(x, d);
  }

  // x = x - thisᵀ·d
  void
  multsVecTrans(Vector &x, const ConstVector &d) const
  {
    ConstGeneralMatrix(*this).multsVecTrans(x, d);
  }

  // this = zero
  void zeros();

  // this = unit (on main diagonal)
  void unit();

  // this = NaN
  void nans();

  // this = ∞
  void infs();

  // this = scalar·this
  void mult(double a);

  // this = this + scalar·m
  void add(double a, const ConstGeneralMatrix &m);
  void
  add(double a, const GeneralMatrix &m)
  {
    add(a, ConstGeneralMatrix(m));
  }

  // this = this + scalar·mᵀ
  void add(double a, const ConstGeneralMatrix &m, const std::string &dum);
  void
  add(double a, const GeneralMatrix &m, const std::string &dum)
  {
    add(a, ConstGeneralMatrix(m), dum);
  }

  bool
  isFinite() const
  {
    return (ConstGeneralMatrix(*this)).isFinite();
  }

  bool
  isZero() const
  {
    return (ConstGeneralMatrix(*this)).isZero();
  }

  virtual void
  print() const
  {
    ConstGeneralMatrix(*this).print();
  }
private:
  void copy(const ConstGeneralMatrix &m, int ioff = 0, int joff = 0);
  void
  copy(const GeneralMatrix &m, int ioff = 0, int joff = 0)
  {
    copy(ConstGeneralMatrix(m), ioff, joff);
  }

  void gemm(const std::string &transa, const ConstGeneralMatrix &a,
            const std::string &transb, const ConstGeneralMatrix &b,
            double alpha, double beta);
  void
  gemm(const std::string &transa, const GeneralMatrix &a,
       const std::string &transb, const GeneralMatrix &b,
       double alpha, double beta)
  {
    gemm(transa, ConstGeneralMatrix(a), transb, ConstGeneralMatrix(b),
         alpha, beta);
  }

  /* this = this * op(m) (without whole copy of this) */
  void gemm_partial_right(const std::string &trans, const ConstGeneralMatrix &m,
                          double alpha, double beta);
  void
  gemm_partial_right(const std::string &trans, const GeneralMatrix &m,
                     double alpha, double beta)
  {
    gemm_partial_right(trans, ConstGeneralMatrix(m), alpha, beta);
  }

  // this = op(m)·this (without whole copy of ‘this’)
  void gemm_partial_left(const std::string &trans, const ConstGeneralMatrix &m,
                         double alpha, double beta);
  void
  gemm_partial_left(const std::string &trans, const GeneralMatrix &m,
                    double alpha, double beta)
  {
    gemm_partial_left(trans, ConstGeneralMatrix(m), alpha, beta);
  }

  /* number of rows/columns for copy used in gemm_partial_* */
  static constexpr int md_length = 23;
};

// Computes a·b
inline GeneralMatrix
operator*(const ConstGeneralMatrix &a, const ConstGeneralMatrix &b)
{
  GeneralMatrix m(a.rows, b.cols);
  m.gemm("N", a, "N", b, 1.0, 0.0);
  return m;
}

// Computes a·bᵀ
template<class T>
GeneralMatrix
operator*(const ConstGeneralMatrix &a, const TransposedMatrix<T> &b)
{
  GeneralMatrix m(a.rows, b.orig.rows);
  m.gemm("N", a, "T", b.orig, 1.0, 0.0);
  return m;
}

// Computes aᵀ·b
template<class T>
GeneralMatrix
operator*(const TransposedMatrix<T> &a, const ConstGeneralMatrix &b)
{
  GeneralMatrix m(a.orig.cols, b.cols);
  m.gemm("T", a.orig, "N", b, 1.0, 0.0);
  return m;
}

// Computes aᵀ·bᵀ
template<class T1, class T2>
GeneralMatrix
operator*(const TransposedMatrix<T1> &a, const TransposedMatrix<T2> &b)
{
  GeneralMatrix m(a.orig.cols, b.orig.rows);
  m.gemm("T", a.orig, "T", b.orig, 1.0, 0.0);
  return m;
}

class SVDDecomp
{
protected:
  // Minimum of number of rows and columns of the decomposed matrix
  const int minmn;
  // Singular values
  Vector sigma;
  // Orthogonal matrix U
  GeneralMatrix U;
  // Orthogonal matrix Vᵀ
  GeneralMatrix VT;
  // Convered flag
  bool conv;
public:
  SVDDecomp(const GeneralMatrix &A)
    : minmn(std::min<int>(A.nrows(), A.ncols())),
      sigma(minmn),
      U(A.nrows(), A.nrows()),
      VT(A.ncols(), A.ncols()),
      conv(false)
  {
    construct(A);
  }
  const GeneralMatrix &
  getU() const
  {
    return U;
  }
  const GeneralMatrix &
  getVT() const
  {
    return VT;
  }
  void solve(const ConstGeneralMatrix &B, GeneralMatrix &X) const;
  void
  solve(const ConstVector &b, Vector &x) const
  {
    GeneralMatrix xmat(x, x.length(), 1);
    solve(ConstGeneralMatrix(b, b.length(), 1), xmat);
  }
private:
  void construct(const GeneralMatrix &A);
};

#endif /* GENERAL_MATRIX_H */
