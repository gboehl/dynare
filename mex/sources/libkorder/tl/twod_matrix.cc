/*
 * Copyright © 2004 Ondra Kamenik
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

#include "twod_matrix.hh"
#include "tl_exception.hh"

#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>

ConstTwoDMatrix::ConstTwoDMatrix(const TwoDMatrix& m) : ConstGeneralMatrix(m)
{
}

ConstTwoDMatrix::ConstTwoDMatrix(const TwoDMatrix& m, int first_col, int num) :
    ConstGeneralMatrix(m, 0, first_col, m.nrows(), num)
{
}

ConstTwoDMatrix::ConstTwoDMatrix(const ConstTwoDMatrix& m, int first_col, int num) :
    ConstGeneralMatrix(m, 0, first_col, m.nrows(), num)
{
}

ConstTwoDMatrix::ConstTwoDMatrix(int first_row, int num, const TwoDMatrix& m) :
    ConstGeneralMatrix(m, first_row, 0, num, m.ncols())
{
}

ConstTwoDMatrix::ConstTwoDMatrix(int first_row, int num, const ConstTwoDMatrix& m) :
    ConstGeneralMatrix(m, first_row, 0, num, m.ncols())
{
}

TwoDMatrix&
TwoDMatrix::operator=(const ConstTwoDMatrix& m)
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
TwoDMatrix::copyRow(const ConstTwoDMatrix& m, int from, int to)
{
  getRow(to) = m.getRow(from);
}

void
TwoDMatrix::addRow(double d, const ConstTwoDMatrix& m, int from, int to)
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
TwoDMatrix::copyColumn(const ConstTwoDMatrix& m, int from, int to)
{
  getCol(to) = m.getCol(from);
}

void
TwoDMatrix::addColumn(double d, const ConstTwoDMatrix& m, int from, int to)
{
  getCol(to).add(d, m.getCol(from));
}

void
TwoDMatrix::save(const std::string& fname) const
{
  std::ofstream fd {fname, std::ios::out | std::ios::trunc};
  if (fd.fail())
    TL_RAISE("Cannot open file for writing in TwoDMatrix::save");

  fd << std::setprecision(std::numeric_limits<double>::max_digits10);

  for (int row = 0; row < nrows(); row++)
    {
      for (int col = 0; col < ncols(); col++)
        fd << ' ' << get(row, col);
      fd << '\n';
    }
  fd.close();
}
