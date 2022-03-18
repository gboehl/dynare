/*
 * Copyright © 2004-2011 Ondra Kamenik
 * Copyright © 2019-2022 Dynare Team
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

#include "Vector.hh"
#include "GeneralMatrix.hh"
#include "SylvException.hh"

#include <dynblas.h>

#include <cmath>
#include <algorithm>
#include <limits>
#include <iostream>
#include <iomanip>

Vector::Vector(const Vector &v)
  : len(v.len), data{new double[len]}
{
  copy(v.data, v.s);
}

Vector::Vector(const ConstVector &v)
  : len(v.len), data{new double[len]}
{
  copy(v.data, v.s);
}

Vector &
Vector::operator=(const Vector &v)
{
  if (this == &v)
    return *this;

  if (v.len != len)
    throw SYLV_MES_EXCEPTION("Attempt to assign vectors with different lengths.");

  if (s == v.s
      && (data <= v.data && v.data < data+len*s
          || v.data <= data && data < v.data+v.len*v.s)
      && (data-v.data) % s == 0)
    throw SYLV_MES_EXCEPTION("Attempt to assign overlapping vectors.");

  copy(v.data, v.s);
  return *this;
}

Vector &
Vector::operator=(Vector &&v)
{
  if (v.len != len)
    throw SYLV_MES_EXCEPTION("Attempt to assign vectors with different lengths.");
  copy(v.data, v.s);
  return *this;
}

Vector &
Vector::operator=(const ConstVector &v)
{
  if (v.len != len)
    throw SYLV_MES_EXCEPTION("Attempt to assign vectors with different lengths.");
  if (s == v.s
      && (data <= v.data && v.data < data+len*s
          || v.data <= data && data < v.data+v.len*v.s)
      && (data-v.data) % s == 0)
    throw SYLV_MES_EXCEPTION("Attempt to assign overlapping vectors.");

  copy(v.data, v.s);
  return *this;
}

void
Vector::copy(const double *d, int inc)
{
  blas_int n = len;
  blas_int incy = s;
  blas_int inc2 = inc;
  dcopy(&n, d, &inc2, data, &incy);
}

Vector::Vector(Vector &v, int off_arg, int l)
  : len(l), s(v.s), data{v.data+off_arg*v.s}, destroy{false}
{
  if (off_arg < 0 || off_arg + len > v.len)
    throw SYLV_MES_EXCEPTION("Subvector not contained in supvector.");
}

Vector::Vector(const Vector &v, int off_arg, int l)
  : len(l), data{new double[len]}
{
  if (off_arg < 0 || off_arg + len > v.len)
    throw SYLV_MES_EXCEPTION("Subvector not contained in supvector.");
  copy(v.data+off_arg*v.s, v.s);
}

Vector::Vector(Vector &v, int off_arg, int skip, int l)
  : len(l), s(v.s*skip), data{v.data+off_arg*v.s}, destroy{false}
{
}

Vector::Vector(const Vector &v, int off_arg, int skip, int l)
  : len(l), data{new double[len]}
{
  copy(v.data+off_arg*v.s, v.s*skip);
}

#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
Vector::Vector(mxArray *p)
  : len{static_cast<int>(mxGetNumberOfElements(p))},
    data{mxGetPr(p)}, destroy{false}
{
  if (!mxIsDouble(p) || mxIsComplex(p) || mxIsSparse(p))
    throw SYLV_MES_EXCEPTION("This is not a dense array of real doubles.");
}
#endif

bool
Vector::operator==(const Vector &y) const
{
  return ConstVector(*this) == y;
}

bool
Vector::operator!=(const Vector &y) const
{
  return ConstVector(*this) != y;
}

bool
Vector::operator<(const Vector &y) const
{
  return ConstVector(*this) < y;
}

bool
Vector::operator<=(const Vector &y) const
{
  return ConstVector(*this) <= y;
}

bool
Vector::operator>(const Vector &y) const
{
  return ConstVector(*this) > y;
}

bool
Vector::operator>=(const Vector &y) const
{
  return ConstVector(*this) >= y;
}

void
Vector::zeros()
{
  if (s == 1)
    std::fill_n(data, len, 0.0);
  else
    for (int i = 0; i < len; i++)
      operator[](i) = 0.0;
}

void
Vector::nans()
{
  for (int i = 0; i < len; i++)
    operator[](i) = std::numeric_limits<double>::quiet_NaN();
}

void
Vector::infs()
{
  for (int i = 0; i < len; i++)
    operator[](i) = std::numeric_limits<double>::infinity();
}

void
Vector::rotatePair(double alpha, double beta1, double beta2, int i)
{
  double tmp = alpha*operator[](i) - beta1*operator[](i+1);
  operator[](i+1) = alpha*operator[](i+1) - beta2*operator[](i);
  operator[](i) = tmp;
}

void
Vector::add(double r, const Vector &v)
{
  add(r, ConstVector(v));
}

void
Vector::add(double r, const ConstVector &v)
{
  blas_int n = len;
  blas_int incx = v.s;
  blas_int incy = s;
  daxpy(&n, &r, v.data, &incx, data, &incy);
}

void
Vector::addComplex(const std::complex<double> &z, const Vector &v)
{
  addComplex(z, ConstVector(v));
}

void
Vector::addComplex(const std::complex<double> &z, const ConstVector &v)
{
  blas_int n = len/2;
  blas_int incx = v.s;
  blas_int incy = s;
  zaxpy(&n, reinterpret_cast<const double(&)[2]>(z), v.data, &incx, data, &incy);
}

void
Vector::mult(double r)
{
  blas_int n = len;
  blas_int incx = s;
  dscal(&n, &r, data, &incx);
}

void
Vector::mult2(double alpha, double beta1, double beta2,
              Vector &x1, Vector &x2,
              const Vector &b1, const Vector &b2)
{
  x1.zeros();
  x2.zeros();
  mult2a(alpha, beta1, beta2, x1, x2, b1, b2);
}

void
Vector::mult2a(double alpha, double beta1, double beta2,
               Vector &x1, Vector &x2,
               const Vector &b1, const Vector &b2)
{
  x1.add(alpha, b1);
  x1.add(-beta1, b2);
  x2.add(alpha, b2);
  x2.add(-beta2, b1);
}

double
Vector::getNorm() const
{
  ConstVector v(*this);
  return v.getNorm();
}

double
Vector::getMax() const
{
  ConstVector v(*this);
  return v.getMax();
}

double
Vector::getNorm1() const
{
  ConstVector v(*this);
  return v.getNorm1();
}

double
Vector::dot(const Vector &y) const
{
  return ConstVector(*this).dot(ConstVector(y));
}

bool
Vector::isFinite() const
{
  return (ConstVector(*this)).isFinite();
}

void
Vector::print() const
{
  auto ff = std::cout.flags();
  std::cout << std::setprecision(4);
  for (int i = 0; i < len; i++)
    std::cout << i << '\t' << std::setw(8) << operator[](i) << std::endl;
  std::cout.flags(ff);
}

ConstVector::ConstVector(const Vector &v)
  : len{v.len}, s{v.s}, data{v.data}
{
}

ConstVector::ConstVector(const ConstVector &v, int off_arg, int l)
  : len{l}, s{v.s}, data{v.data+off_arg*v.s}
{
  if (off_arg < 0 || off_arg + len > v.len)
    throw SYLV_MES_EXCEPTION("Subvector not contained in supvector.");
}

ConstVector::ConstVector(const ConstVector &v, int off_arg, int skip, int l)
  : len(l), s{v.s*skip}, data{v.data+off_arg*v.s}
{
}

ConstVector::ConstVector(const double *d, int skip, int l)
  : len{l}, s{skip}, data{d}
{
}

#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
ConstVector::ConstVector(const mxArray *p)
  : len{static_cast<int>(mxGetNumberOfElements(p))},
    data{mxGetPr(p)}
{
  if (!mxIsDouble(p))
    throw SYLV_MES_EXCEPTION("This is not a MATLAB array of doubles.");
}
#endif

bool
ConstVector::operator==(const ConstVector &y) const
{
  if (len != y.len)
    return false;
  if (len == 0)
    return true;
  int i = 0;
  while (i < len && operator[](i) == y[i])
    i++;
  return i == len;
}

bool
ConstVector::operator<(const ConstVector &y) const
{
  int i = std::min(len, y.len);
  int ii = 0;
  while (ii < i && operator[](ii) == y[ii])
    ii++;
  if (ii < i)
    return operator[](ii) < y[ii];
  else
    return len < y.len;
}

double
ConstVector::getNorm() const
{
  double s = 0;
  for (int i = 0; i < len; i++)
    s += operator[](i)*operator[](i);
  return sqrt(s);
}

double
ConstVector::getMax() const
{
  double r = 0;
  for (int i = 0; i < len; i++)
    r = std::max(r, std::abs(operator[](i)));
  return r;
}

double
ConstVector::getNorm1() const
{
  double norm = 0.0;
  for (int i = 0; i < len; i++)
    norm += std::abs(operator[](i));
  return norm;
}

double
ConstVector::dot(const ConstVector &y) const
{
  if (len != y.len)
    throw SYLV_MES_EXCEPTION("Vector has different length in ConstVector::dot.");
  blas_int n = len;
  blas_int incx = s;
  blas_int incy = y.s;
  return ddot(&n, data, &incx, y.data, &incy);
}

bool
ConstVector::isFinite() const
{
  int i = 0;
  while (i < len && std::isfinite(operator[](i)))
    i++;
  return i == len;
}

void
ConstVector::print() const
{
  auto ff = std::cout.flags();
  std::cout << std::setprecision(4);
  for (int i = 0; i < len; i++)
    std::cout << i << '\t' << std::setw(8) << operator[](i) << std::endl;
  std::cout.flags(ff);
}
