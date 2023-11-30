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

#include "GeneralMatrix.hh"
#include "SylvException.hh"

#include <dynblas.h>
#include <dynlapack.h>

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>

GeneralMatrix::GeneralMatrix(const GeneralMatrix& m) :
    data(m.rows * m.cols), rows(m.rows), cols(m.cols), ld(m.rows)
{
  copy(m);
}

GeneralMatrix::GeneralMatrix(const ConstGeneralMatrix& m) :
    data(m.rows * m.cols), rows(m.rows), cols(m.cols), ld(m.rows)
{
  copy(m);
}

GeneralMatrix::GeneralMatrix(const GeneralMatrix& m, int i, int j, int nrows, int ncols) :
    data(nrows * ncols), rows(nrows), cols(ncols), ld(nrows)
{
  copy(m, i, j);
}

GeneralMatrix::GeneralMatrix(GeneralMatrix& m, int i, int j, int nrows, int ncols) :
    data(m.data, m.ld * j + i, m.ld * (ncols - 1) + nrows), rows(nrows), cols(ncols), ld(m.ld)
{
}

GeneralMatrix&
GeneralMatrix::operator=(const ConstGeneralMatrix& m)
{
  data = m.data;
  rows = m.rows;
  cols = m.cols;
  ld = m.ld;
  return *this;
}

Vector
GeneralMatrix::getRow(int row)
{
  return Vector {data, row, ld, cols};
}

Vector
GeneralMatrix::getCol(int col)
{
  return Vector {data, ld * col, rows};
}

ConstVector
GeneralMatrix::getRow(int row) const
{
  return ConstVector {data, row, ld, cols};
}

ConstVector
GeneralMatrix::getCol(int col) const
{
  return ConstVector {data, ld * col, rows};
}

void
GeneralMatrix::place(const ConstGeneralMatrix& m, int i, int j)
{
  if (i + m.nrows() > nrows() || j + m.ncols() > ncols())
    throw SYLV_MES_EXCEPTION("Bad submatrix placement, matrix dimensions exceeded.");

  GeneralMatrix tmpsub(*this, i, j, m.nrows(), m.ncols());
  tmpsub.copy(m);
}

void
GeneralMatrix::mult(const ConstGeneralMatrix& a, const ConstGeneralMatrix& b)
{
  gemm("N", a, "N", b, 1.0, 0.0);
}

void
GeneralMatrix::multAndAdd(const ConstGeneralMatrix& a, const ConstGeneralMatrix& b, double mult)
{
  gemm("N", a, "N", b, mult, 1.0);
}

void
GeneralMatrix::multAndAdd(const ConstGeneralMatrix& a, const ConstGeneralMatrix& b,
                          [[maybe_unused]] const std::string& dum, double mult)
{
  gemm("N", a, "T", b, mult, 1.0);
}

void
GeneralMatrix::multAndAdd(const ConstGeneralMatrix& a, [[maybe_unused]] const std::string& dum,
                          const ConstGeneralMatrix& b, double mult)
{
  gemm("T", a, "N", b, mult, 1.0);
}

void
GeneralMatrix::multAndAdd(const ConstGeneralMatrix& a, [[maybe_unused]] const std::string& dum1,
                          const ConstGeneralMatrix& b, [[maybe_unused]] const std::string& dum2,
                          double mult)
{
  gemm("T", a, "T", b, mult, 1.0);
}

void
GeneralMatrix::addOuter(const ConstVector& a, double mult)
{
  if (nrows() != ncols())
    throw SYLV_MES_EXCEPTION("Matrix is not square in GeneralMatrix::addOuter.");
  if (nrows() != a.length())
    throw SYLV_MES_EXCEPTION("Wrong length of a vector in GeneralMatrix::addOuter.");

  // since BLAS dsyr (symmetric rank 1 update) assumes symmetricity, we do this manually
  for (int i = 0; i < nrows(); i++)
    for (int j = i; j < nrows(); j++)
      {
        double s = mult * a[i] * a[j];
        get(i, j) = get(i, j) + s;
        if (i != j)
          get(j, i) = get(j, i) + s;
      }
}

void
GeneralMatrix::multRight(const ConstGeneralMatrix& m)
{
  gemm_partial_right("N", m, 1.0, 0.0);
}

void
GeneralMatrix::multLeft(const ConstGeneralMatrix& m)
{
  gemm_partial_left("N", m, 1.0, 0.0);
}

void
GeneralMatrix::multRightTrans(const ConstGeneralMatrix& m)
{
  gemm_partial_right("T", m, 1.0, 0.0);
}

void
GeneralMatrix::multLeftTrans(const ConstGeneralMatrix& m)
{
  gemm_partial_left("T", m, 1.0, 0.0);
}

// here we must be careful for ld
void
GeneralMatrix::zeros()
{
  if (ld == rows)
    data.zeros();
  else
    for (int i = 0; i < rows; i++)
      for (int j = 0; j < cols; j++)
        get(i, j) = 0;
}

void
GeneralMatrix::unit()
{
  for (int i = 0; i < rows; i++)
    for (int j = 0; j < cols; j++)
      if (i == j)
        get(i, j) = 1.0;
      else
        get(i, j) = 0.0;
}

void
GeneralMatrix::nans()
{
  for (int i = 0; i < rows; i++)
    for (int j = 0; j < cols; j++)
      get(i, j) = std::numeric_limits<double>::quiet_NaN();
}

void
GeneralMatrix::infs()
{
  for (int i = 0; i < rows; i++)
    for (int j = 0; j < cols; j++)
      get(i, j) = std::numeric_limits<double>::infinity();
}

// here we must be careful for ld
void
GeneralMatrix::mult(double a)
{
  if (ld == rows)
    data.mult(a);
  else
    for (int i = 0; i < rows; i++)
      for (int j = 0; j < cols; j++)
        get(i, j) *= a;
}

// here we must be careful for ld
void
GeneralMatrix::add(double a, const ConstGeneralMatrix& m)
{
  if (m.nrows() != rows || m.ncols() != cols)
    throw SYLV_MES_EXCEPTION("Matrix has different size in GeneralMatrix::add.");

  if (ld == rows && m.ld == m.rows)
    data.add(a, m.data);
  else
    for (int i = 0; i < rows; i++)
      for (int j = 0; j < cols; j++)
        get(i, j) += a * m.get(i, j);
}

void
GeneralMatrix::add(double a, const ConstGeneralMatrix& m, [[maybe_unused]] const std::string& dum)
{
  if (m.nrows() != cols || m.ncols() != rows)
    throw SYLV_MES_EXCEPTION("Matrix has different size in GeneralMatrix::add.");

  for (int i = 0; i < rows; i++)
    for (int j = 0; j < cols; j++)
      get(i, j) += a * m.get(j, i);
}

void
GeneralMatrix::copy(const ConstGeneralMatrix& m, int ioff, int joff)
{
  for (int i = 0; i < rows; i++)
    for (int j = 0; j < cols; j++)
      get(i, j) = m.get(i + ioff, j + joff);
}

void
GeneralMatrix::gemm(const std::string& transa, const ConstGeneralMatrix& a,
                    const std::string& transb, const ConstGeneralMatrix& b, double alpha,
                    double beta)
{
  int opa_rows = a.nrows();
  int opa_cols = a.ncols();
  if (transa == "T")
    {
      opa_rows = a.ncols();
      opa_cols = a.nrows();
    }
  int opb_rows = b.nrows();
  int opb_cols = b.ncols();
  if (transb == "T")
    {
      opb_rows = b.ncols();
      opb_cols = b.nrows();
    }

  if (opa_rows != nrows() || opb_cols != ncols() || opa_cols != opb_rows)
    throw SYLV_MES_EXCEPTION("Wrong dimensions for matrix multiplication.");

  blas_int m = opa_rows;
  blas_int n = opb_cols;
  blas_int k = opa_cols;
  blas_int lda = a.ld;
  blas_int ldb = b.ld;
  blas_int ldc = ld;
  if (lda > 0 && ldb > 0 && ldc > 0)
    dgemm(transa.c_str(), transb.c_str(), &m, &n, &k, &alpha, a.data.base(), &lda, b.data.base(),
          &ldb, &beta, data.base(), &ldc);
  else if (nrows() * ncols() > 0)
    {
      if (beta == 0.0)
        zeros();
      else
        mult(beta);
    }
}

void
GeneralMatrix::gemm_partial_left(const std::string& trans, const ConstGeneralMatrix& m,
                                 double alpha, double beta)
{
  int icol;
  for (icol = 0; icol + md_length < cols; icol += md_length)
    {
      GeneralMatrix incopy(const_cast<const GeneralMatrix&>(*this), 0, icol, rows, md_length);
      GeneralMatrix inplace(*this, 0, icol, rows, md_length);
      inplace.gemm(trans, m, "N", ConstGeneralMatrix(incopy), alpha, beta);
    }
  if (cols > icol)
    {
      GeneralMatrix incopy(const_cast<const GeneralMatrix&>(*this), 0, icol, rows, cols - icol);
      GeneralMatrix inplace(*this, 0, icol, rows, cols - icol);
      inplace.gemm(trans, m, "N", ConstGeneralMatrix(incopy), alpha, beta);
    }
}

void
GeneralMatrix::gemm_partial_right(const std::string& trans, const ConstGeneralMatrix& m,
                                  double alpha, double beta)
{
  int irow;
  for (irow = 0; irow + md_length < rows; irow += md_length)
    {
      GeneralMatrix incopy(const_cast<const GeneralMatrix&>(*this), irow, 0, md_length, cols);
      GeneralMatrix inplace(*this, irow, 0, md_length, cols);
      inplace.gemm("N", ConstGeneralMatrix(incopy), trans, m, alpha, beta);
    }
  if (rows > irow)
    {
      GeneralMatrix incopy(const_cast<const GeneralMatrix&>(*this), irow, 0, rows - irow, cols);
      GeneralMatrix inplace(*this, irow, 0, rows - irow, cols);
      inplace.gemm("N", ConstGeneralMatrix(incopy), trans, m, alpha, beta);
    }
}

ConstGeneralMatrix::ConstGeneralMatrix(const GeneralMatrix& m, int i, int j, int nrows, int ncols) :
    data(m.getData(), j * m.getLD() + i, (ncols - 1) * m.getLD() + nrows),
    rows(nrows),
    cols(ncols),
    ld(m.getLD())
{
  // FIXME: check that the submatrix is fully in the matrix
}

ConstGeneralMatrix::ConstGeneralMatrix(const ConstGeneralMatrix& m, int i, int j, int nrows,
                                       int ncols) :
    data(m.getData(), j * m.getLD() + i, (ncols - 1) * m.getLD() + nrows),
    rows(nrows),
    cols(ncols),
    ld(m.getLD())
{
  // FIXME: check that the submatrix is fully in the matrix
}

ConstGeneralMatrix::ConstGeneralMatrix(const GeneralMatrix& m) :
    data(m.data), rows(m.rows), cols(m.cols), ld(m.ld)
{
}

ConstVector
ConstGeneralMatrix::getRow(int row) const
{
  return ConstVector {data, row, ld, cols};
}

ConstVector
ConstGeneralMatrix::getCol(int col) const
{
  return ConstVector {data, ld * col, rows};
}

double
ConstGeneralMatrix::getNormInf() const
{
  double norm = 0.0;
  for (int i = 0; i < nrows(); i++)
    {
      double normi = getRow(i).getNorm1();
      norm = std::max(normi, norm);
    }
  return norm;
}

double
ConstGeneralMatrix::getNorm1() const
{
  double norm = 0.0;
  for (int j = 0; j < ncols(); j++)
    {
      double normj = getCol(j).getNorm1();
      norm = std::max(normj, norm);
    }
  return norm;
}

void
ConstGeneralMatrix::multVec(double a, Vector& x, double b, const ConstVector& d) const
{
  if (x.length() != rows || cols != d.length())
    throw SYLV_MES_EXCEPTION("Wrong dimensions for vector multiply.");
  if (rows > 0)
    {
      blas_int mm = rows;
      blas_int nn = cols;
      double alpha = b;
      blas_int lda = ld;
      blas_int incx = d.skip();
      double beta = a;
      blas_int incy = x.skip();
      dgemv("N", &mm, &nn, &alpha, data.base(), &lda, d.base(), &incx, &beta, x.base(), &incy);
    }
}

void
ConstGeneralMatrix::multVecTrans(double a, Vector& x, double b, const ConstVector& d) const
{
  if (x.length() != cols || rows != d.length())
    throw SYLV_MES_EXCEPTION("Wrong dimensions for vector multiply.");
  if (rows > 0)
    {
      blas_int mm = rows;
      blas_int nn = cols;
      double alpha = b;
      blas_int lda = ld;
      blas_int incx = d.skip();
      double beta = a;
      blas_int incy = x.skip();
      dgemv("T", &mm, &nn, &alpha, data.base(), &lda, d.base(), &incx, &beta, x.base(), &incy);
    }
}

/* m = inv(this)*m */
void
ConstGeneralMatrix::multInvLeft(const std::string& trans, int mrows, int mcols, int mld,
                                double* d) const
{
  if (rows != cols)
    throw SYLV_MES_EXCEPTION("The matrix is not square for inversion.");
  if (cols != mrows)
    throw SYLV_MES_EXCEPTION("Wrong dimensions for matrix inverse mutliply.");

  if (rows > 0)
    {
      GeneralMatrix inv(*this);
      auto ipiv = std::make_unique<lapack_int[]>(rows);
      lapack_int info;
      lapack_int rows2 = rows, mcols2 = mcols, mld2 = mld, lda = inv.ld;
      dgetrf(&rows2, &rows2, inv.getData().base(), &lda, ipiv.get(), &info);
      dgetrs(trans.c_str(), &rows2, &mcols2, inv.base(), &lda, ipiv.get(), d, &mld2, &info);
    }
}

/* m = inv(this)*m */
void
ConstGeneralMatrix::multInvLeft(GeneralMatrix& m) const
{
  multInvLeft("N", m.nrows(), m.ncols(), m.getLD(), m.getData().base());
}

/* m = inv(this')*m */
void
ConstGeneralMatrix::multInvLeftTrans(GeneralMatrix& m) const
{
  multInvLeft("T", m.nrows(), m.ncols(), m.getLD(), m.getData().base());
}

/* d = inv(this)*d */
void
ConstGeneralMatrix::multInvLeft(Vector& d) const
{
  if (d.skip() != 1)
    throw SYLV_MES_EXCEPTION("Skip≠1 not implemented in ConstGeneralMatrix::multInvLeft(Vector&)");

  multInvLeft("N", d.length(), 1, d.length(), d.base());
}

/* d = inv(this')*d */
void
ConstGeneralMatrix::multInvLeftTrans(Vector& d) const
{
  if (d.skip() != 1)
    throw SYLV_MES_EXCEPTION("Skip≠1 not implemented in ConstGeneralMatrix::multInvLeft(Vector&)");

  multInvLeft("T", d.length(), 1, d.length(), d.base());
}

bool
ConstGeneralMatrix::isFinite() const
{
  for (int i = 0; i < nrows(); i++)
    for (int j = 0; j < ncols(); j++)
      if (!std::isfinite(get(i, j)))
        return false;
  return true;
}

bool
ConstGeneralMatrix::isZero() const
{
  for (int i = 0; i < nrows(); i++)
    for (int j = 0; j < ncols(); j++)
      if (get(i, j) != 0.0)
        return false;
  return true;
}

void
ConstGeneralMatrix::print() const
{
  auto ff = std::cout.flags();
  std::cout << "rows=" << rows << ", cols=" << cols << std::endl;
  for (int i = 0; i < rows; i++)
    {
      std::cout << "row " << i << ':' << std::endl << std::setprecision(3);
      for (int j = 0; j < cols; j++)
        std::cout << std::setw(6) << get(i, j) << ' ';
      std::cout << std::endl;
    }
  std::cout.flags(ff);
}

void
SVDDecomp::construct(const GeneralMatrix& A)
{
  // quick exit if empty matrix
  if (minmn == 0)
    {
      U.unit();
      VT.unit();
      conv = true;
      return;
    }

  // make copy of the matrix
  GeneralMatrix AA(A);

  lapack_int m = AA.nrows();
  lapack_int n = AA.ncols();
  double* a = AA.base();
  lapack_int lda = AA.getLD();
  double* s = sigma.base();
  double* u = U.base();
  lapack_int ldu = U.getLD();
  double* vt = VT.base();
  lapack_int ldvt = VT.getLD();
  double tmpwork;
  lapack_int lwork = -1;
  lapack_int info;

  auto iwork = std::make_unique<lapack_int[]>(8 * minmn);
  // query for optimal lwork
  dgesdd("A", &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, &tmpwork, &lwork, iwork.get(), &info);
  lwork = static_cast<lapack_int>(tmpwork);
  Vector work(lwork);
  // do the decomposition
  dgesdd("A", &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work.base(), &lwork, iwork.get(), &info);
  if (info < 0)
    throw SYLV_MES_EXCEPTION("Internal error in SVDDecomp constructor");
  if (info == 0)
    conv = true;
}

void
SVDDecomp::solve(const ConstGeneralMatrix& B, GeneralMatrix& X) const
{
  if (B.nrows() != U.nrows())
    throw SYLV_MES_EXCEPTION("Incompatible number of rows ");

  /* Reciprocal condition number for determination of zeros in the
     end of sigma */
  constexpr double rcond = 1e-13;

  // determine nz=number of zeros in the end of sigma
  int nz = 0;
  while (nz < minmn && sigma[minmn - 1 - nz] < rcond * sigma[0])
    nz++;
  // take relevant B for sigma inversion
  int m = U.nrows();
  int n = VT.ncols();
  GeneralMatrix Bprime(transpose(U) * B, m - minmn, 0, minmn - nz, B.ncols());
  // solve sigma
  for (int i = 0; i < minmn - nz; i++)
    Bprime.getRow(i).mult(1.0 / sigma[i]);
  // solve VT
  X.zeros();
  //- copy Bprime to right place of X
  for (int i = 0; i < minmn - nz; i++)
    X.getRow(n - minmn + i) = Bprime.getRow(i);
  //- multiply with VT
  X.multLeftTrans(VT);
}
