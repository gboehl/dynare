/* $Header: /var/lib/cvs/dynare_cpp/sylv/cc/GeneralMatrix.h,v 1.3 2004/11/24 20:41:59 kamenik Exp $ */

/* Tag $Name:  $ */

#ifndef GENERAL_MATRIX_H
#define GENERAL_MATRIX_H

#include "Vector.hh"
#include "SylvException.hh"

#include <algorithm>
#include <memory>
#include <utility>
#include <string>

class GeneralMatrix;

class ConstGeneralMatrix
{
  friend class GeneralMatrix;
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
  // Implicit conversion from ConstGeneralMatrix is ok, since it's cheap
  ConstGeneralMatrix(const GeneralMatrix &m);
  ConstGeneralMatrix(const GeneralMatrix &m, int i, int j, int nrows, int ncols);
  ConstGeneralMatrix(const ConstGeneralMatrix &m, int i, int j, int nrows, int ncols);
  virtual ~ConstGeneralMatrix() = default;

  ConstGeneralMatrix &operator=(const ConstGeneralMatrix &v) = delete;
  ConstGeneralMatrix &operator=(ConstGeneralMatrix &&v) = delete;

  const double &
  get(int i, int j) const
  {
    return data[j*ld+i];
  }
  int
  numRows() const
  {
    return rows;
  }
  int
  numCols() const
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
  GeneralMatrix(const GeneralMatrix &m, const std::string &dummy); // transpose
  GeneralMatrix(const ConstGeneralMatrix &m, const std::string &dummy); // transpose
  GeneralMatrix(const GeneralMatrix &m, int i, int j, int nrows, int ncols);
  GeneralMatrix(GeneralMatrix &m, int i, int j, int nrows, int ncols);
  /* this = a*b */
  GeneralMatrix(const ConstGeneralMatrix &a, const ConstGeneralMatrix &b);
  /* this = a*b' */
  GeneralMatrix(const ConstGeneralMatrix &a, const ConstGeneralMatrix &b, const std::string &dum);
  /* this = a'*b */
  GeneralMatrix(const ConstGeneralMatrix &a, const std::string &dum, const ConstGeneralMatrix &b);
  /* this = a'*b */
  GeneralMatrix(const ConstGeneralMatrix &a, const std::string &dum1,
                const ConstGeneralMatrix &b, const std::string &dum2);

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
  numRows() const
  {
    return rows;
  }
  int
  numCols() const
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

  /* this = a*b */
  void mult(const ConstGeneralMatrix &a, const ConstGeneralMatrix &b);
  void
  mult(const GeneralMatrix &a, const GeneralMatrix &b)
  {
    mult(ConstGeneralMatrix(a), ConstGeneralMatrix(b));
  }

  /* this = this + scalar*a*b */
  void multAndAdd(const ConstGeneralMatrix &a, const ConstGeneralMatrix &b,
                  double mult = 1.0);
  void
  multAndAdd(const GeneralMatrix &a, const GeneralMatrix &b,
             double mult = 1.0)
  {
    multAndAdd(ConstGeneralMatrix(a), ConstGeneralMatrix(b), mult);
  }

  /* this = this + scalar*a*b' */
  void multAndAdd(const ConstGeneralMatrix &a, const ConstGeneralMatrix &b,
                  const std::string &dum, double mult = 1.0);
  void
  multAndAdd(const GeneralMatrix &a, const GeneralMatrix &b,
             const std::string &dum, double mult = 1.0)
  {
    multAndAdd(ConstGeneralMatrix(a), ConstGeneralMatrix(b), dum, mult);
  }

  /* this = this + scalar*a'*b */
  void multAndAdd(const ConstGeneralMatrix &a, const std::string &dum, const ConstGeneralMatrix &b,
                  double mult = 1.0);
  void
  multAndAdd(const GeneralMatrix &a, const std::string &dum, const GeneralMatrix &b,
             double mult = 1.0)
  {
    multAndAdd(ConstGeneralMatrix(a), dum, ConstGeneralMatrix(b), mult);
  }

  /* this = this + scalar*a'*b' */
  void multAndAdd(const ConstGeneralMatrix &a, const std::string &dum1,
                  const ConstGeneralMatrix &b, const std::string &dum2, double mult = 1.0);
  void
  multAndAdd(const GeneralMatrix &a, const std::string &dum1,
             const GeneralMatrix &b, const std::string &dum2, double mult = 1.0)
  {
    multAndAdd(ConstGeneralMatrix(a), dum1, ConstGeneralMatrix(b), dum2, mult);
  }

  /* this = this + scalar*a*a' */
  void addOuter(const ConstVector &a, double mult = 1.0);
  void
  addOuter(const Vector &a, double mult = 1.0)
  {
    addOuter(ConstVector(a), mult);
  }

  /* this = this * m */
  void multRight(const ConstGeneralMatrix &m);
  void
  multRight(const GeneralMatrix &m)
  {
    multRight(ConstGeneralMatrix(m));
  }

  /* this = m * this */
  void multLeft(const ConstGeneralMatrix &m);
  void
  multLeft(const GeneralMatrix &m)
  {
    multLeft(ConstGeneralMatrix(m));
  }

  /* this = this * m' */
  void multRightTrans(const ConstGeneralMatrix &m);
  void
  multRightTrans(const GeneralMatrix &m)
  {
    multRightTrans(ConstGeneralMatrix(m));
  }

  /* this = m' * this */
  void multLeftTrans(const ConstGeneralMatrix &m);
  void
  multLeftTrans(const GeneralMatrix &m)
  {
    multLeftTrans(ConstGeneralMatrix(m));
  }

  /* x = scalar(a)*x + scalar(b)*this*d */
  void
  multVec(double a, Vector &x, double b, const ConstVector &d) const
  {
    ConstGeneralMatrix(*this).multVec(a, x, b, d);
  }

  /* x = scalar(a)*x + scalar(b)*this'*d */
  void
  multVecTrans(double a, Vector &x, double b, const ConstVector &d) const
  {
    ConstGeneralMatrix(*this).multVecTrans(a, x, b, d);
  }

  /* x = x + this*d */
  void
  multaVec(Vector &x, const ConstVector &d) const
  {
    ConstGeneralMatrix(*this).multaVec(x, d);
  }

  /* x = x + this'*d */
  void
  multaVecTrans(Vector &x, const ConstVector &d) const
  {
    ConstGeneralMatrix(*this).multaVecTrans(x, d);
  }

  /* x = x - this*d */
  void
  multsVec(Vector &x, const ConstVector &d) const
  {
    ConstGeneralMatrix(*this).multsVec(x, d);
  }

  /* x = x - this'*d */
  void
  multsVecTrans(Vector &x, const ConstVector &d) const
  {
    ConstGeneralMatrix(*this).multsVecTrans(x, d);
  }

  /* this = zero */
  void zeros();

  /** this = unit (on main diagonal) */
  void unit();

  /* this = NaN */
  void nans();

  /* this = Inf */
  void infs();

  /* this = scalar*this */
  void mult(double a);

  /* this = this + scalar*m */
  void add(double a, const ConstGeneralMatrix &m);
  void
  add(double a, const GeneralMatrix &m)
  {
    add(a, ConstGeneralMatrix(m));
  }

  /* this = this + scalar*m' */
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

  /* this = op(m) *this (without whole copy of this) */
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

class SVDDecomp
{
protected:
  /** Minimum of number of rows and columns of the decomposed
   * matrix. */
  const int minmn;
  /** Singular values. */
  Vector sigma;
  /** Orthogonal matrix U. */
  GeneralMatrix U;
  /** Orthogonal matrix V^T. */
  GeneralMatrix VT;
  /** Convered flag. */
  bool conv;
public:
  SVDDecomp(const GeneralMatrix &A)
    : minmn(std::min<int>(A.numRows(), A.numCols())),
      sigma(minmn),
      U(A.numRows(), A.numRows()),
      VT(A.numCols(), A.numCols()),
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
