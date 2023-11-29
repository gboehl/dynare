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

#include "SylvMatrix.hh"
#include "SylvException.hh"
#include "int_power.hh"

#include <dynblas.h>
#include <dynlapack.h>

#include <cmath>
#include <memory>

void
SylvMatrix::multLeftI(const SqSylvMatrix& m)
{
  int off = rows - m.nrows();
  if (off < 0)
    throw SYLV_MES_EXCEPTION("Wrong matrix dimensions for multLeftI.");

  GeneralMatrix subtmp(*this, off, 0, m.nrows(), cols);
  subtmp.multLeft(m);
}

void
SylvMatrix::multLeftITrans(const SqSylvMatrix& m)
{
  int off = rows - m.nrows();
  if (off < 0)
    throw SYLV_MES_EXCEPTION("Wrong matrix dimensions for multLeftITrans.");

  GeneralMatrix subtmp(*this, off, 0, m.nrows(), cols);
  subtmp.multLeftTrans(m);
}

void
SylvMatrix::multLeft(int zero_cols, const GeneralMatrix& a, const GeneralMatrix& b)
{
  int off = a.nrows() - a.ncols();
  if (off < 0 || a.nrows() != rows || off != zero_cols || rows != b.nrows() || cols != b.ncols())
    throw SYLV_MES_EXCEPTION("Wrong matrix dimensions for multLeft.");

  /* Here we cannot call SylvMatrix::gemm() since it would require
     another copy of (usually big) b (we are not able to do inplace
     submatrix of const GeneralMatrix) */
  if (a.getLD() > 0 && ld > 0)
    {
      blas_int mm = a.nrows();
      blas_int nn = cols;
      blas_int kk = a.ncols();
      double alpha = 1.0;
      blas_int lda = a.getLD();
      blas_int ldb = ld;
      double beta = 0.0;
      blas_int ldc = ld;
      dgemm("N", "N", &mm, &nn, &kk, &alpha, a.getData().base(), &lda, b.getData().base() + off,
            &ldb, &beta, data.base(), &ldc);
    }
}

void
SylvMatrix::multRightKron(const SqSylvMatrix& m, int order)
{
  if (power(m.nrows(), order) != cols)
    throw SYLV_MES_EXCEPTION("Wrong number of cols for right kron multiply.");

  KronVector auxrow(m.nrows(), m.nrows(), order - 1);
  for (int i = 0; i < rows; i++)
    {
      Vector rowi {getRow(i)};
      KronVector rowikron(rowi, m.nrows(), m.nrows(), order - 1);
      auxrow = rowi; // copy data
      m.multVecKronTrans(rowikron, auxrow);
    }
}

void
SylvMatrix::multRightKronTrans(const SqSylvMatrix& m, int order)
{
  if (power(m.nrows(), order) != cols)
    throw SYLV_MES_EXCEPTION("Wrong number of cols for right kron multiply.");

  KronVector auxrow(m.nrows(), m.nrows(), order - 1);
  for (int i = 0; i < rows; i++)
    {
      Vector rowi {getRow(i)};
      KronVector rowikron(rowi, m.nrows(), m.nrows(), order - 1);
      auxrow = rowi; // copy data
      m.multVecKron(rowikron, auxrow);
    }
}

void
SylvMatrix::eliminateLeft(int row, int col, Vector& x)
{
  double d = get(col, col);
  double e = get(row, col);
  if (std::abs(d) > std::abs(e))
    {
      get(row, col) = 0.0;
      double mult = e / d;
      for (int i = col + 1; i < ncols(); i++)
        get(row, i) = get(row, i) - mult * get(col, i);
      x[row] = x[row] - mult * x[col];
    }
  else if (std::abs(e) > std::abs(d))
    {
      get(row, col) = 0.0;
      get(col, col) = e;
      double mult = d / e;
      for (int i = col + 1; i < ncols(); i++)
        {
          double tx = get(col, i);
          double ty = get(row, i);
          get(col, i) = ty;
          get(row, i) = tx - mult * ty;
        }
      double tx = x[col];
      double ty = x[row];
      x[col] = ty;
      x[row] = tx - mult * ty;
    }
}

void
SylvMatrix::eliminateRight(int row, int col, Vector& x)
{
  double d = get(row, row);
  double e = get(row, col);

  if (std::abs(d) > std::abs(e))
    {
      get(row, col) = 0.0;
      double mult = e / d;
      for (int i = 0; i < row; i++)
        get(i, col) = get(i, col) - mult * get(i, row);
      x[col] = x[col] - mult * x[row];
    }
  else if (std::abs(e) > std::abs(d))
    {
      get(row, col) = 0.0;
      get(row, row) = e;
      double mult = d / e;
      for (int i = 0; i < row; i++)
        {
          double tx = get(i, row);
          double ty = get(i, col);
          get(i, row) = ty;
          get(i, col) = tx - mult * ty;
        }
      double tx = x[row];
      double ty = x[col];
      x[row] = ty;
      x[col] = tx - mult * ty;
    }
}

void
SqSylvMatrix::multVecKron(KronVector& x, const ConstKronVector& d) const
{
  x.zeros();
  if (d.getDepth() == 0)
    multaVec(x, d);
  else
    {
      KronVector aux(x.getM(), x.getN(), x.getDepth());
      for (int i = 0; i < x.getM(); i++)
        {
          KronVector auxi(aux, i);
          ConstKronVector di(d, i);
          multVecKron(auxi, di);
        }
      for (int i = 0; i < rows; i++)
        {
          KronVector xi(x, i);
          for (int j = 0; j < cols; j++)
            {
              KronVector auxj(aux, j);
              xi.add(get(i, j), auxj);
            }
        }
    }
}

void
SqSylvMatrix::multVecKronTrans(KronVector& x, const ConstKronVector& d) const
{
  x.zeros();
  if (d.getDepth() == 0)
    multaVecTrans(x, d);
  else
    {
      KronVector aux(x.getM(), x.getN(), x.getDepth());
      for (int i = 0; i < x.getM(); i++)
        {
          KronVector auxi(aux, i);
          ConstKronVector di(d, i);
          multVecKronTrans(auxi, di);
        }
      for (int i = 0; i < rows; i++)
        {
          KronVector xi(x, i);
          for (int j = 0; j < cols; j++)
            {
              KronVector auxj(aux, j);
              xi.add(get(j, i), auxj);
            }
        }
    }
}

void
SqSylvMatrix::multInvLeft2(GeneralMatrix& a, GeneralMatrix& b, double& rcond1,
                           double& rcondinf) const
{
  if (rows != a.nrows() || rows != b.nrows())
    throw SYLV_MES_EXCEPTION("Wrong dimensions for multInvLeft2.");

  // PLU factorization
  Vector inv(data);
  auto ipiv = std::make_unique<lapack_int[]>(rows);
  lapack_int info;
  lapack_int rows2 = rows, lda = ld;
  dgetrf(&rows2, &rows2, inv.base(), &lda, ipiv.get(), &info);
  // solve a
  lapack_int acols = a.ncols();
  double* abase = a.base();
  dgetrs("N", &rows2, &acols, inv.base(), &lda, ipiv.get(), abase, &rows2, &info);
  // solve b
  lapack_int bcols = b.ncols();
  double* bbase = b.base();
  dgetrs("N", &rows2, &bcols, inv.base(), &lda, ipiv.get(), bbase, &rows2, &info);

  // condition numbers
  auto work = std::make_unique<double[]>(4 * rows);
  auto iwork = std::make_unique<lapack_int[]>(rows);
  double norm1 = getNorm1();
  dgecon("1", &rows2, inv.base(), &lda, &norm1, &rcond1, work.get(), iwork.get(), &info);
  double norminf = getNormInf();
  dgecon("I", &rows2, inv.base(), &lda, &norminf, &rcondinf, work.get(), iwork.get(), &info);
}

void
SqSylvMatrix::setUnit()
{
  for (int i = 0; i < rows; i++)
    for (int j = 0; j < cols; j++)
      if (i == j)
        get(i, j) = 1.0;
      else
        get(i, j) = 0.0;
}
