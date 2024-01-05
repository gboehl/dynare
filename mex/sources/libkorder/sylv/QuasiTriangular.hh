/*
 * Copyright © 2004-2011 Ondra Kamenik
 * Copyright © 2019-2024 Dynare Team
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

#ifndef QUASI_TRIANGULAR_HH
#define QUASI_TRIANGULAR_HH

#include "KronVector.hh"
#include "SylvMatrix.hh"
#include "Vector.hh"

#include <list>
#include <memory>

class DiagonalBlock;
class Diagonal;
class DiagPair
{
private:
  double* a1;
  double* a2;

public:
  DiagPair() = default;
  DiagPair(double* aa1, double* aa2) : a1 {aa1}, a2 {aa2}
  {
  }
  DiagPair(const DiagPair& p) = default;
  DiagPair& operator=(const DiagPair& p) = default;
  DiagPair&
  operator=(double v)
  {
    *a1 = v;
    *a2 = v;
    return *this;
  }
  const double&
  operator*() const
  {
    return *a1;
  }
  /* Here we must not define double& operator*(), since it wouldn't
     rewrite both values, we use operator=() for this */
  friend class Diagonal;
  friend class DiagonalBlock;
};

/* Stores a diagonal block of a quasi-triangular real matrix:
   – either a 1×1 block, i.e. a real scalar, stored in α₁
                               ⎛α₁ β₁⎞
   – or a 2×2 block, stored as ⎝β₂ α₂⎠
*/
class DiagonalBlock
{
private:
  int jbar; // Index of block in the diagonal
  bool real;
  DiagPair alpha;
  double* beta1;
  double* beta2;

public:
  DiagonalBlock() = default;
  DiagonalBlock(int jb, bool r, double* a1, double* a2, double* b1, double* b2) :
      jbar {jb}, real {r}, alpha {a1, a2}, beta1 {b1}, beta2 {b2}
  {
  }
  // Construct a complex 2×2 block
  /* β₁ and β₂ will be deduced from pointers to α₁ and α₂ */
  DiagonalBlock(int jb, double* a1, double* a2) :
      jbar {jb}, real {false}, alpha {a1, a2}, beta1 {a2 - 1}, beta2 {a1 + 1}
  {
  }
  // Construct a real 1×1 block
  DiagonalBlock(int jb, double* a1) :
      jbar {jb}, real {true}, alpha {a1, a1}, beta1 {nullptr}, beta2 {nullptr}
  {
  }
  DiagonalBlock(const DiagonalBlock& b) = default;
  DiagonalBlock& operator=(const DiagonalBlock& b) = default;
  [[nodiscard]] int
  getIndex() const
  {
    return jbar;
  }
  [[nodiscard]] bool
  isReal() const
  {
    return real;
  }
  [[nodiscard]] const DiagPair&
  getAlpha() const
  {
    return alpha;
  }
  DiagPair&
  getAlpha()
  {
    return alpha;
  }
  [[nodiscard]] double&
  getBeta1() const
  {
    return *beta1;
  }
  [[nodiscard]] double&
  getBeta2() const
  {
    return *beta2;
  }
  // Returns determinant of this block (assuming it is 2×2)
  [[nodiscard]] double getDeterminant() const;
  // Returns −β₁β₂
  [[nodiscard]] double getSBeta() const;
  // Returns the modulus of the eigenvalue(s) contained in this block
  [[nodiscard]] double getSize() const;
  // Transforms this block into a real one
  void setReal();
  // Verifies that the block information is consistent with the matrix d (for debugging)
  void checkBlock(const double* d, int d_size);
  friend class Diagonal;
};

// Stores the diagonal blocks of a quasi-triangular real matrix
class Diagonal
{
public:
  using const_diag_iter = std::list<DiagonalBlock>::const_iterator;
  using diag_iter = std::list<DiagonalBlock>::iterator;

private:
  int num_all {0}; // Total number of blocks
  std::list<DiagonalBlock> blocks;
  int num_real {0}; // Number of 1×1 (real) blocks
public:
  Diagonal() = default;
  // Construct the diagonal blocks of (quasi-triangular) matrix ‘data’
  Diagonal(double* data, int d_size);
  /* Construct the diagonal blocks of (quasi-triangular) matrix ‘data’,
     assuming it has the same shape as ‘d’ */
  Diagonal(double* data, const Diagonal& d);
  Diagonal(const Diagonal& d) = default;
  Diagonal& operator=(const Diagonal& d) = default;
  virtual ~Diagonal() = default;

  // Returns number of 2×2 blocks on the diagonal
  [[nodiscard]] int
  getNumComplex() const
  {
    return num_all - num_real;
  }
  // Returns number of 1×1 blocks on the diagonal
  [[nodiscard]] int
  getNumReal() const
  {
    return num_real;
  }
  // Returns number of scalar elements on the diagonal
  [[nodiscard]] int
  getSize() const
  {
    return getNumReal() + 2 * getNumComplex();
  }
  // Returns total number of blocks on the diagonal
  [[nodiscard]] int
  getNumBlocks() const
  {
    return num_all;
  }
  void getEigenValues(Vector& eig) const;
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
  [[nodiscard]] const_diag_iter
  begin() const
  {
    return blocks.begin();
  }
  diag_iter
  end()
  {
    return blocks.end();
  }
  [[nodiscard]] const_diag_iter
  end() const
  {
    return blocks.end();
  }

  /* redefine pointers as data start at p */
  void changeBase(double* p);

private:
  constexpr static double EPS = 1.0e-300;
  /* Computes number of 2×2 diagonal blocks on the quasi-triangular matrix
     represented by data (of size d_size×d_size) */
  static int getNumComplex(const double* data, int d_size);
  // Checks whether |p|<EPS
  static bool isZero(double p);
};

template<class _TRef, class _TPtr>
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
  operator==(const _Self& it) const
  {
    return ptr == it.ptr;
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
  virtual _Self& operator++() = 0;
};

template<class _TRef, class _TPtr>
class _column_iter : public _matrix_iter<_TRef, _TPtr>
{
  using _Tparent = _matrix_iter<_TRef, _TPtr>;
  using _Self = _column_iter<_TRef, _TPtr>;
  int row;

public:
  _column_iter(_TPtr base, int ds, bool r, int rw) :
      _matrix_iter<_TRef, _TPtr>(base, ds, r), row(rw) {};
  _Self&
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
      return *(_Tparent::ptr + _Tparent::d_size);
  }
  [[nodiscard]] int
  getRow() const
  {
    return row;
  }
};

template<class _TRef, class _TPtr>
class _row_iter : public _matrix_iter<_TRef, _TPtr>
{
  using _Tparent = _matrix_iter<_TRef, _TPtr>;
  using _Self = _row_iter<_TRef, _TPtr>;
  int col;

public:
  _row_iter(_TPtr base, int ds, bool r, int cl) :
      _matrix_iter<_TRef, _TPtr>(base, ds, r), col(cl) {};
  _Self&
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
      return *(_Tparent::ptr + 1);
  }
  [[nodiscard]] int
  getCol() const
  {
    return col;
  }
};

class SchurDecomp;
class SchurDecompZero;

/* Represents an upper quasi-triangular matrix.
   All the elements are stored in the SqSylvMatrix super-class.
   Additionally, a list of the diagonal blocks (1×1 or 2×2), is stored in the
   “diagonal” member, in order to optimize some operations (where the matrix is
   seen as an upper-triangular matrix, plus sub-diagonal elements of the 2×2
   diagonal blocks) */
class QuasiTriangular : public SqSylvMatrix
{
public:
  using const_col_iter = _column_iter<const double&, const double*>;
  using col_iter = _column_iter<double&, double*>;
  using const_row_iter = _row_iter<const double&, const double*>;
  using row_iter = _row_iter<double&, double*>;
  using const_diag_iter = Diagonal::const_diag_iter;
  using diag_iter = Diagonal::diag_iter;

protected:
  Diagonal diagonal;

public:
  QuasiTriangular(const ConstVector& d, int d_size);
  // Initializes with r·t
  QuasiTriangular(double r, const QuasiTriangular& t);
  // Initializes with r·t+r₂·t₂
  QuasiTriangular(double r, const QuasiTriangular& t, double r2, const QuasiTriangular& t2);
  // Initializes with t²
  QuasiTriangular(const std::string& dummy, const QuasiTriangular& t);
  explicit QuasiTriangular(const SchurDecomp& decomp);
  explicit QuasiTriangular(const SchurDecompZero& decomp);
  QuasiTriangular(const QuasiTriangular& t);

  ~QuasiTriangular() override = default;
  [[nodiscard]] const Diagonal&
  getDiagonal() const
  {
    return diagonal;
  }
  [[nodiscard]] int getNumOffdiagonal() const;
  void swapDiagLogically(diag_iter it);
  void checkDiagConsistency(diag_iter it);
  double getAverageDiagSize(diag_iter start, diag_iter end);
  diag_iter findClosestDiagBlock(diag_iter start, diag_iter end, double a);
  diag_iter findNextLargerBlock(diag_iter start, diag_iter end, double a);

  /* (I+this)·y = x, y→x  */
  virtual void solvePre(Vector& x, double& eig_min);
  /* (I+thisᵀ)·y = x, y→x */
  virtual void solvePreTrans(Vector& x, double& eig_min);
  /* (I+this)·x = b */
  virtual void solve(Vector& x, const ConstVector& b, double& eig_min);
  /* (I+thisᵀ)·x = b */
  virtual void solveTrans(Vector& x, const ConstVector& b, double& eig_min);
  /* x = this·b */
  virtual void multVec(Vector& x, const ConstVector& b) const;
  /* x = thisᵀ·b */
  virtual void multVecTrans(Vector& x, const ConstVector& b) const;
  /* x = x + this·b */
  virtual void multaVec(Vector& x, const ConstVector& b) const;
  /* x = x + thisᵀ·b */
  virtual void multaVecTrans(Vector& x, const ConstVector& b) const;
  /* x = (this⊗I)·x */
  virtual void multKron(KronVector& x) const;
  /* x = (thisᵀ⊗I)·x */
  virtual void multKronTrans(KronVector& x) const;
  /* A = this·A */
  virtual void multLeftOther(GeneralMatrix& a) const;
  /* A = thisᵀ·A */
  virtual void multLeftOtherTrans(GeneralMatrix& a) const;

  [[nodiscard]] const_diag_iter
  diag_begin() const
  {
    return diagonal.begin();
  }
  diag_iter
  diag_begin()
  {
    return diagonal.begin();
  }
  [[nodiscard]] const_diag_iter
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
  [[nodiscard]] virtual const_col_iter col_begin(const DiagonalBlock& b) const;
  virtual col_iter col_begin(const DiagonalBlock& b);
  [[nodiscard]] virtual const_row_iter row_begin(const DiagonalBlock& b) const;
  virtual row_iter row_begin(const DiagonalBlock& b);
  [[nodiscard]] virtual const_col_iter col_end(const DiagonalBlock& b) const;
  virtual col_iter col_end(const DiagonalBlock& b);
  [[nodiscard]] virtual const_row_iter row_end(const DiagonalBlock& b) const;
  virtual row_iter row_end(const DiagonalBlock& b);

  [[nodiscard]] virtual std::unique_ptr<QuasiTriangular>
  clone() const
  {
    return std::make_unique<QuasiTriangular>(*this);
  }
  // Returns this²
  [[nodiscard]] virtual std::unique_ptr<QuasiTriangular>
  square() const
  {
    return std::make_unique<QuasiTriangular>("square", *this);
  }
  // Returns r·this
  [[nodiscard]] virtual std::unique_ptr<QuasiTriangular>
  scale(double r) const
  {
    return std::make_unique<QuasiTriangular>(r, *this);
  }
  // Returns r·this + r₂·t₂
  [[nodiscard]] virtual std::unique_ptr<QuasiTriangular>
  linearlyCombine(double r, double r2, const QuasiTriangular& t2) const
  {
    return std::make_unique<QuasiTriangular>(r, *this, r2, t2);
  }

protected:
  // this = r·t
  void setMatrix(double r, const QuasiTriangular& t);
  // this = this + r·t
  void addMatrix(double r, const QuasiTriangular& t);

private:
  // this = this + I
  void addUnit();
  /* x = x + (this⊗I)·b */
  void multaKron(KronVector& x, const ConstKronVector& b) const;
  /* x = x + (thisᵀ⊗I)·b */
  void multaKronTrans(KronVector& x, const ConstKronVector& b) const;
  /* hide noneffective implementations of parents */
  void multsVec(Vector& x, const ConstVector& d) const;
  void multsVecTrans(Vector& x, const ConstVector& d) const;
};

#endif
