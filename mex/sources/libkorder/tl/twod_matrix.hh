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

// Matrix interface.

/* Here we make an interface to 2-dimensional matrix defined in the Sylvester
   module. That abstraction provides an interface to BLAS. The main purpose of
   this file is to only make its subclass in order to keep the tensor library
   and Sylvester module independent. So here is mainly renaming of methods.

   Similarly as in the Sylvester module we declare two classes TwoDMatrix and
   ConstTwoDMatrix. The only purpose of the latter is to allow submatrix
   construction from const reference arguments. */

#ifndef TWOD_MATRIX_H
#define TWOD_MATRIX_H

#include "GeneralMatrix.hh"

#include <matio.h>

#include <string>
#include <utility>

class TwoDMatrix;

class ConstTwoDMatrix : public ConstGeneralMatrix
{
public:
  ConstTwoDMatrix(int m, int n, ConstVector d)
    : ConstGeneralMatrix(std::move(d), m, n)
  {
  }
  ConstTwoDMatrix(const ConstTwoDMatrix &m) = default;
  ConstTwoDMatrix(ConstTwoDMatrix &&m) = default;

  // Implicit conversion from TwoDMatrix is ok, since it's cheap
  ConstTwoDMatrix(const TwoDMatrix &m);

  // Constructors creating a submatrix of consecutive columns
  ConstTwoDMatrix(const TwoDMatrix &m, int first_col, int num);
  ConstTwoDMatrix(const ConstTwoDMatrix &m, int first_col, int num);

  // Constructors creating a submatrix of consecutive rows
  ConstTwoDMatrix(int first_row, int num, const TwoDMatrix &m);
  ConstTwoDMatrix(int first_row, int num, const ConstTwoDMatrix &m);

  ConstTwoDMatrix(const ConstTwoDMatrix &m, int first_row, int first_col, int rows, int cols)
    : ConstGeneralMatrix(m, first_row, first_col, rows, cols)
  {
  }
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
  explicit ConstTwoDMatrix(const mxArray *p) : ConstGeneralMatrix(p)
  {
  }
#endif
  ~ConstTwoDMatrix() override = default;

  ConstTwoDMatrix &operator=(const ConstTwoDMatrix &v) = delete;
  ConstTwoDMatrix &operator=(ConstTwoDMatrix &&v) = delete;

  void writeMat(mat_t *fd, const std::string &vname) const;
};

class TwoDMatrix : public GeneralMatrix
{
public:
  TwoDMatrix(const TwoDMatrix &m) = default;
  TwoDMatrix(TwoDMatrix &&m) = default;
  TwoDMatrix(int r, int c)
    : GeneralMatrix(r, c)
  {
  }
  TwoDMatrix(int r, int c, Vector d)
    : GeneralMatrix(std::move(d), r, c)
  {
  }
  TwoDMatrix(const GeneralMatrix &m)
    : GeneralMatrix(m)
  {
  }
  // We don't want implict conversion from ConstGeneralMatrix, since it's expensive
  explicit TwoDMatrix(const ConstGeneralMatrix &m)
    : GeneralMatrix(m)
  {
  }
  template<class T>
  explicit TwoDMatrix(const TransposedMatrix<T> &m)
    : GeneralMatrix(m)
  {
  }
  // Select only some columns (with data copy)
  TwoDMatrix(const TwoDMatrix &m, int first_col, int num)
    : GeneralMatrix(m, 0, first_col, m.nrows(), num)
  {
  }
  // Select only some columns (with data sharing)
  TwoDMatrix(TwoDMatrix &m, int first_col, int num)
    : GeneralMatrix(m, 0, first_col, m.nrows(), num)
  {
  }
  // Select only some rows (with data copy)
  TwoDMatrix(int first_row, int num, const TwoDMatrix &m)
    : GeneralMatrix(m, first_row, 0, num, m.ncols())
  {
  }
  // Select only some rows (with data sharing)
  TwoDMatrix(int first_row, int num, TwoDMatrix &m)
    : GeneralMatrix(m, first_row, 0, num, m.ncols())
  {
  }
  // Select a submatrix (with data sharing)
  TwoDMatrix(TwoDMatrix &m, int first_row, int first_col, int rows, int cols)
    : GeneralMatrix(m, first_row, first_col, rows, cols)
  {
  }
  // Select a submatrix (with data copy)
  TwoDMatrix(const TwoDMatrix &m, int first_row, int first_col, int rows, int cols)
    : GeneralMatrix(m, first_row, first_col, rows, cols)
  {
  }

#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
  explicit TwoDMatrix(mxArray *p) : GeneralMatrix(p)
  {
  }
#endif
  ~TwoDMatrix() override = default;

  TwoDMatrix &operator=(const TwoDMatrix &m) = default;
  TwoDMatrix &operator=(TwoDMatrix &&m) = default;
  TwoDMatrix &operator=(const ConstTwoDMatrix &m);

  // TwoDMatrix row methods declarations
  void copyRow(int from, int to);
  void copyRow(const ConstTwoDMatrix &m, int from, int to);
  void
  copyRow(const TwoDMatrix &m, int from, int to)
  {
    copyRow(ConstTwoDMatrix(m), from, to);
  }
  void
  addRow(const ConstTwoDMatrix &m, int from, int to)
  {
    addRow(1.0, m, from, to);
  }
  void
  addRow(const TwoDMatrix &m, int from, int to)
  {
    addRow(1.0, ConstTwoDMatrix(m), from, to);
  }
  void addRow(double d, const ConstTwoDMatrix &m, int from, int to);
  void
  addRow(double d, const TwoDMatrix &m, int from, int to)
  {
    addRow(d, ConstTwoDMatrix(m), from, to);
  }

  // TwoDMatrix column methods declarations
  void copyColumn(int from, int to);
  void copyColumn(const ConstTwoDMatrix &m, int from, int to);
  void
  copyColumn(const TwoDMatrix &m, int from, int to)
  {
    copyColumn(ConstTwoDMatrix(m), from, to);
  }
  void
  addColumn(const ConstTwoDMatrix &m, int from, int to)
  {
    addColumn(1.0, ConstTwoDMatrix(m), from, to);
  }
  void
  addColumn(const TwoDMatrix &m, int from, int to)
  {
    addColumn(1.0, ConstTwoDMatrix(m), from, to);
  }
  void addColumn(double d, const ConstTwoDMatrix &m, int from, int to);
  void
  addColumn(double d, const TwoDMatrix &m, int from, int to)
  {
    addColumn(d, ConstTwoDMatrix(m), from, to);
  }

  // Saves the matrix to a text file
  void save(const std::string &fname) const;

  void
  writeMat(mat_t *fd, const std::string &vname) const
  {
    ConstTwoDMatrix(*this).writeMat(fd, vname);
  }
};

#endif
