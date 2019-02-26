// Copyright (C) 2004-2011, Ondra Kamenik

// Matrix interface.

/* Here we make an interface to 2-dimensional matrix defined in the
   Sylvester module. That abstraction provides an interface to BLAS. The
   main purpose of this file is to only make its subclass in order to
   keep the tensor library and Sylvester module independent. So here is
   mainly renaming of methods.

   Similarly as in the Sylvester module we declare two classes
   |TwoDMatrix| and |ConstTwoDMatrix|. The only purpose of the latter is
   to allow submatrix construction from const reference arguments. */

#ifndef TWOD_MATRIX_H
#define TWOD_MATRIX_H

#include "GeneralMatrix.hh"

#include <matio.h>

#include <string>
#include <utility>

class TwoDMatrix;

/* We make two obvious constructors, and then a constructor making
   submatrix of subsequent columns. We also rename
   |GeneralMatrix::numRows()| and |GeneralMatrix::numCols()|. */

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
  ConstTwoDMatrix(const TwoDMatrix &m, int first_col, int num);
  ConstTwoDMatrix(const ConstTwoDMatrix &m, int first_col, int num);
  ConstTwoDMatrix(int first_row, int num, const TwoDMatrix &m);
  ConstTwoDMatrix(int first_row, int num, const ConstTwoDMatrix &m);
  ConstTwoDMatrix(const ConstTwoDMatrix &m, int first_row, int first_col, int rows, int cols)
    : ConstGeneralMatrix(m, first_row, first_col, rows, cols)
  {
  }
  ~ConstTwoDMatrix() override = default;

  ConstTwoDMatrix &operator=(const ConstTwoDMatrix &v) = delete;
  ConstTwoDMatrix &operator=(ConstTwoDMatrix &&v) = delete;

  int
  nrows() const
  {
    return numRows();
  }
  int
  ncols() const
  {
    return numCols();
  }
  void writeMat(mat_t *fd, const std::string &vname) const;
};

/* Here we do the same as for |ConstTwoDMatrix| plus define
   methods for copying and adding rows and columns.

   Also we have |save| method which dumps the matrix to a file with a
   given name. The file can be read by Scilab {\tt fscanfMat} function. */

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
    : GeneralMatrix(m, 0, first_col, m.numRows(), num)
  {
  }
  // Select only some columns (with data sharing)
  TwoDMatrix(TwoDMatrix &m, int first_col, int num)
    : GeneralMatrix(m, 0, first_col, m.numRows(), num)
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
  TwoDMatrix(const ConstTwoDMatrix &a, const ConstTwoDMatrix &b)
    : GeneralMatrix(a, b)
  {
  }
  ~TwoDMatrix() override = default;

  TwoDMatrix &operator=(const TwoDMatrix &m) = default;
  TwoDMatrix &operator=(TwoDMatrix &&m) = default;
  TwoDMatrix &operator=(const ConstTwoDMatrix &m);

  int
  nrows() const
  {
    return numRows();
  }
  int
  ncols() const
  {
    return numCols();
  }

  // |TwoDMatrix| row methods declarations
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

  // |TwoDMatrix| column methods declarations
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

  void save(const std::string &fname) const;
  void
  writeMat(mat_t *fd, const std::string &vname) const
  {
    ConstTwoDMatrix(*this).writeMat(fd, vname);
  }
};

#endif
