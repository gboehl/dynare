/* $Header: /var/lib/cvs/dynare_cpp/sylv/cc/QuasiTriangular.h,v 1.1.1.1 2004/06/04 13:00:44 kamenik Exp $ */

/* Tag $Name:  $ */

#ifndef QUASI_TRIANGULAR_H
#define QUASI_TRIANGULAR_H

#include "Vector.hh"
#include "KronVector.hh"
#include "SylvMatrix.hh"

#include <list>
#include <memory>

class DiagonalBlock;
class Diagonal;
class DiagPair
{
private:
  double *a1;
  double *a2;
public:
  DiagPair() = default;
  DiagPair(double *aa1, double *aa2) : a1{aa1}, a2{aa2}
  {
  }
  DiagPair(const DiagPair &p) = default;
  DiagPair &operator=(const DiagPair &p) = default;
  DiagPair &
  operator=(double v)
  {
    *a1 = v;
    *a2 = v;
    return *this;
  }
  const double &
  operator*() const
  {
    return *a1;
  }
  /** here we must not define double& operator*(), since it wouldn't
      rewrite both values, we use operator= for this */
  friend class Diagonal;
  friend class DiagonalBlock;
};

// Stores a diagonal block: either a scalar, or a 2x2 block
/* alpha points to the diagonal element(s); beta1 and beta2 point to the
   off-diagonal elements of the 2x2 block */
class DiagonalBlock
{
private:
  int jbar;
  bool real;
  DiagPair alpha;
  double *beta1;
  double *beta2;

public:
  DiagonalBlock() = default;
  DiagonalBlock(int jb, bool r, double *a1, double *a2,
                double *b1, double *b2)
    : jbar{jb}, real{r}, alpha{a1, a2}, beta1{b1}, beta2{b2}
  {
  }
  // construct complex block
  DiagonalBlock(int jb, double *a1, double *a2)
    : jbar{jb}, real{false}, alpha{a1, a2}, beta1{a2-1}, beta2{a1+1}
  {
  }
  // construct real block
  DiagonalBlock(int jb, double *a1)
    : jbar{jb}, real{true}, alpha{a1, a1}, beta1{nullptr}, beta2{nullptr}
  {
  }
  DiagonalBlock(const DiagonalBlock &b) = default;
  DiagonalBlock &operator=(const DiagonalBlock &b) = default;
  int
  getIndex() const
  {
    return jbar;
  }
  bool
  isReal() const
  {
    return real;
  }
  const DiagPair &
  getAlpha() const
  {
    return alpha;
  }
  DiagPair &
  getAlpha()
  {
    return alpha;
  }
  double &
  getBeta1() const
  {
    return *beta1;
  }
  double &
  getBeta2() const
  {
    return *beta2;
  }
  double getDeterminant() const;
  double getSBeta() const;
  double getSize() const;
  void setReal();
  // for debugging
  void checkBlock(const double *d, int d_size);
  friend class Diagonal;
};

class Diagonal
{
public:
  using const_diag_iter = std::list<DiagonalBlock>::const_iterator;
  using diag_iter = std::list<DiagonalBlock>::iterator;
private:
  int num_all{0};
  std::list<DiagonalBlock> blocks;
  int num_real{0};
public:
  Diagonal() = default;
  Diagonal(double *data, int d_size);
  Diagonal(double *data, const Diagonal &d);
  Diagonal(const Diagonal &d) = default;
  Diagonal &operator=(const Diagonal &d) = default;
  virtual ~Diagonal() = default;

  int
  getNumComplex() const
  {
    return num_all - num_real;
  }
  int
  getNumReal() const
  {
    return num_real;
  }
  int
  getSize() const
  {
    return getNumReal() + 2*getNumComplex();
  }
  int
  getNumBlocks() const
  {
    return num_all;
  }
  void getEigenValues(Vector &eig) const;
  void swapLogically(diag_iter it);
  void checkConsistency(diag_iter it);
  double getAverageSize(diag_iter start, diag_iter end);
  diag_iter findClosestBlock(diag_iter start, diag_iter end, double a);
  diag_iter findNextLargerBlock(diag_iter start, diag_iter end, double a);
  void print() const;

  diag_iter
  begin()
  {
    return blocks.begin();
  }
  const_diag_iter
  begin() const
  {
    return blocks.begin();
  }
  diag_iter
  end()
  {
    return blocks.end();
  }
  const_diag_iter
  end() const
  {
    return blocks.end();
  }

  /* redefine pointers as data start at p */
  void changeBase(double *p);
private:
  constexpr static double EPS = 1.0e-300;
  static int getNumComplex(const double *data, int d_size);
  static bool isZero(double p);
};

template <class _TRef, class _TPtr>
struct _matrix_iter
{
  using _Self = _matrix_iter<_TRef, _TPtr>;
  int d_size;
  bool real;
  _TPtr ptr;
public:
  _matrix_iter(_TPtr base, int ds, bool r)
  {
    ptr = base;
    d_size = ds;
    real = r;
  }
  virtual ~_matrix_iter() = default;
  bool
  operator==(const _Self &it) const
  {
    return ptr == it.ptr;
  }
  bool
  operator!=(const _Self &it) const
  {
    return ptr != it.ptr;
  }
  _TRef
  operator*() const
  {
    return *ptr;
  }
  _TRef
  a() const
  {
    return *ptr;
  }
  virtual _Self &operator++() = 0;
};

template <class _TRef, class _TPtr>
class _column_iter : public _matrix_iter<_TRef, _TPtr>
{
  using _Tparent = _matrix_iter<_TRef, _TPtr>;
  using _Self = _column_iter<_TRef, _TPtr>;
  int row;
public:
  _column_iter(_TPtr base, int ds, bool r, int rw)
    : _matrix_iter<_TRef, _TPtr>(base, ds, r), row(rw)
  {
  };
  _Self &
  operator++() override
  {
    _Tparent::ptr++;
    row++;
    return *this;
  }
  _TRef
  b() const
  {
    if (_Tparent::real)
      return *(_Tparent::ptr);
    else
      return *(_Tparent::ptr+_Tparent::d_size);
  }
  int
  getRow() const
  {
    return row;
  }
};

template <class _TRef, class _TPtr>
class _row_iter : public _matrix_iter<_TRef, _TPtr>
{
  using _Tparent = _matrix_iter<_TRef, _TPtr>;
  using _Self = _row_iter<_TRef, _TPtr>;
  int col;
public:
  _row_iter(_TPtr base, int ds, bool r, int cl)
    : _matrix_iter<_TRef, _TPtr>(base, ds, r), col(cl)
  {
  };
  _Self &
  operator++() override
  {
    _Tparent::ptr += _Tparent::d_size;
    col++;
    return *this;
  }
  virtual _TRef
  b() const
  {
    if (_Tparent::real)
      return *(_Tparent::ptr);
    else
      return *(_Tparent::ptr+1);
  }
  int
  getCol() const
  {
    return col;
  }
};

class SchurDecomp;
class SchurDecompZero;

/* Represents an upper quasi-triangular matrix.
   All the elements are stored in the SqSylvMatrix super-class.
   Additionally, a list of the diagonal blocks (1x1 or 2x2), is stored in the
   "diagonal" member, in order to optimize some operations (where the matrix is
   seen as an upper-triangular matrix, plus sub-diagonal elements of the 2x2
   diagonal blocks) */
class QuasiTriangular : public SqSylvMatrix
{
public:
  using const_col_iter = _column_iter<const double &, const double *>;
  using col_iter = _column_iter<double &, double *>;
  using const_row_iter = _row_iter<const double &, const double *>;
  using row_iter = _row_iter<double &, double *>;
  using const_diag_iter = Diagonal::const_diag_iter;
  using diag_iter = Diagonal::diag_iter;
protected:
  Diagonal diagonal;
public:
  QuasiTriangular(const ConstVector &d, int d_size);
  QuasiTriangular(double r, const QuasiTriangular &t);
  QuasiTriangular(double r, const QuasiTriangular &t,
                  double rr, const QuasiTriangular &tt);
  QuasiTriangular(int p, const QuasiTriangular &t);
  QuasiTriangular(const SchurDecomp &decomp);
  QuasiTriangular(const SchurDecompZero &decomp);
  QuasiTriangular(const QuasiTriangular &t);
  
  ~QuasiTriangular() override = default;
  const Diagonal &
  getDiagonal() const
  {
    return diagonal;
  }
  int getNumOffdiagonal() const;
  void swapDiagLogically(diag_iter it);
  void checkDiagConsistency(diag_iter it);
  double getAverageDiagSize(diag_iter start, diag_iter end);
  diag_iter findClosestDiagBlock(diag_iter start, diag_iter end, double a);
  diag_iter findNextLargerBlock(diag_iter start, diag_iter end, double a);

  /* (I+T)y = x, y-->x  */
  virtual void solvePre(Vector &x, double &eig_min);
  /* (I+T')y = x, y-->x */
  virtual void solvePreTrans(Vector &x, double &eig_min);
  /* (I+T)x = b */
  virtual void solve(Vector &x, const ConstVector &b, double &eig_min);
  /* (I+T')x = b */
  virtual void solveTrans(Vector &x, const ConstVector &b, double &eig_min);
  /* x = Tb */
  virtual void multVec(Vector &x, const ConstVector &b) const;
  /* x = T'b */
  virtual void multVecTrans(Vector &x, const ConstVector &b) const;
  /* x = x + Tb */
  virtual void multaVec(Vector &x, const ConstVector &b) const;
  /* x = x + T'b */
  virtual void multaVecTrans(Vector &x, const ConstVector &b) const;
  /* x = (T\otimes I)x */
  virtual void multKron(KronVector &x) const;
  /* x = (T'\otimes I)x */
  virtual void multKronTrans(KronVector &x) const;
  /* A = T*A */
  virtual void multLeftOther(GeneralMatrix &a) const;
  /* A = T'*A */
  virtual void multLeftOtherTrans(GeneralMatrix &a) const;

  const_diag_iter
  diag_begin() const
  {
    return diagonal.begin();
  }
  diag_iter
  diag_begin()
  {
    return diagonal.begin();
  }
  const_diag_iter
  diag_end() const
  {
    return diagonal.end();
  }
  diag_iter
  diag_end()
  {
    return diagonal.end();
  }

  /* iterators for off diagonal elements */
  virtual const_col_iter col_begin(const DiagonalBlock &b) const;
  virtual col_iter col_begin(const DiagonalBlock &b);
  virtual const_row_iter row_begin(const DiagonalBlock &b) const;
  virtual row_iter row_begin(const DiagonalBlock &b);
  virtual const_col_iter col_end(const DiagonalBlock &b) const;
  virtual col_iter col_end(const DiagonalBlock &b);
  virtual const_row_iter row_end(const DiagonalBlock &b) const;
  virtual row_iter row_end(const DiagonalBlock &b);

  /* clone */
  virtual std::unique_ptr<QuasiTriangular>
  clone() const
  {
    return std::make_unique<QuasiTriangular>(*this);
  }
  virtual std::unique_ptr<QuasiTriangular>
  clone(int p, const QuasiTriangular &t) const
  {
    return std::make_unique<QuasiTriangular>(p, t);
  }
  virtual std::unique_ptr<QuasiTriangular>
  clone(double r) const
  {
    return std::make_unique<QuasiTriangular>(r, *this);
  }
  virtual std::unique_ptr<QuasiTriangular>
  clone(double r, double rr, const QuasiTriangular &tt) const
  {
    return std::make_unique<QuasiTriangular>(r, *this, rr, tt);
  }
protected:
  void setMatrix(double r, const QuasiTriangular &t);
  void addMatrix(double r, const QuasiTriangular &t);
private:
  void addUnit();
  /* x = x + (T\otimes I)b */
  void multaKron(KronVector &x, const ConstKronVector &b) const;
  /* x = x + (T'\otimes I)b */
  void multaKronTrans(KronVector &x, const ConstKronVector &b) const;
  /* implementation via iterators, useful for large matrices */
  void setMatrixViaIter(double r, const QuasiTriangular &t);
  void addMatrixViaIter(double r, const QuasiTriangular &t);
  /* hide noneffective implementations of parents */
  void multsVec(Vector &x, const ConstVector &d) const;
  void multsVecTrans(Vector &x, const ConstVector &d) const;
};

#endif /* QUASI_TRIANGULAR_H */
