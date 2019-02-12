// Copyright 2004, Ondra Kamenik

#include "twod_matrix.hh"
#include "tl_exception.hh"

#include <vector>
#include <fstream>
#include <iomanip>
#include <limits>

ConstTwoDMatrix::ConstTwoDMatrix(const TwoDMatrix &m)
  : ConstGeneralMatrix(m)
{
}

ConstTwoDMatrix::ConstTwoDMatrix(const TwoDMatrix &m, int first_col, int num)
  : ConstGeneralMatrix(m, 0, first_col, m.nrows(), num)
{
}

ConstTwoDMatrix::ConstTwoDMatrix(const ConstTwoDMatrix &m, int first_col, int num)
  : ConstGeneralMatrix(m, 0, first_col, m.nrows(), num)
{
}

ConstTwoDMatrix::ConstTwoDMatrix(int first_row, int num, const TwoDMatrix &m)
  : ConstGeneralMatrix(m, first_row, 0, num, m.ncols())
{
}

ConstTwoDMatrix::ConstTwoDMatrix(int first_row, int num, const ConstTwoDMatrix &m)
  : ConstGeneralMatrix(m, first_row, 0, num, m.ncols())
{
}

void
ConstTwoDMatrix::writeMat(mat_t *fd, const std::string &vname) const
{
  size_t dims[2];
  dims[0] = nrows();
  dims[1] = ncols();
  std::vector<double> data(nrows()*ncols());

  for (int j = 0; j < ncols(); j++)
    for (int i = 0; i < nrows(); i++)
      data[j*nrows()+i] = get(i, j);

  matvar_t *v = Mat_VarCreate(vname.c_str(), MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, data.data(), 0);

  Mat_VarWrite(fd, v, MAT_COMPRESSION_NONE);

  Mat_VarFree(v);
}

TwoDMatrix &
TwoDMatrix::operator=(const ConstTwoDMatrix &m)
{
  GeneralMatrix::operator=(m);
  return *this;
}

void
TwoDMatrix::copyRow(int from, int to)
{
  if (from != to)
    copyRow(ConstTwoDMatrix(*this), from, to);
}

void
TwoDMatrix::copyRow(const ConstTwoDMatrix &m, int from, int to)
{
  getRow(to) = m.getRow(from);
}

void
TwoDMatrix::addRow(double d, const ConstTwoDMatrix &m, int from, int to)
{
  getRow(to).add(d, m.getRow(from));
}

void
TwoDMatrix::copyColumn(int from, int to)
{
  if (from != to)
    copyColumn(ConstTwoDMatrix(*this), from, to);
}

void
TwoDMatrix::copyColumn(const ConstTwoDMatrix &m, int from, int to)
{
  getCol(to) = m.getCol(from);
}

void
TwoDMatrix::addColumn(double d, const ConstTwoDMatrix &m, int from, int to)
{
  getCol(to).add(d, m.getCol(from));
}

void
TwoDMatrix::save(const std::string &fname) const
{
  std::ofstream fd{fname, std::ios::out | std::ios::trunc};
  if (fd.fail())
    TL_RAISE("Cannot open file for writing in TwoDMatrix::save");

  fd << std::setprecision(std::numeric_limits<double>::digits10 + 1);

  for (int row = 0; row < nrows(); row++)
    {
      for (int col = 0; col < ncols(); col++)
        fd << ' ' << get(row, col);
      fd << '\n';
    }
  fd.close();
}
