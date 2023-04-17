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

#ifndef MM_MATRIX_H
#define MM_MATRIX_H

#include "GeneralMatrix.hh"

#include <string>
#include <utility>
#include <vector>

class MMException
{
  std::string message;
public:
  MMException(std::string mes) : message(std::move(mes))
  {
  }
  std::string
  getMessage() const
  {
    return message;
  }
};

class MMMatrixIn
{
  std::vector<double> data;
  int rows;
  int cols;
public:
  MMMatrixIn(const std::string &fname);
  ~MMMatrixIn() = default;
  Vector
  getData()
  {
    return Vector{data.data(), size()};
  }
  int
  size() const
  {
    return rows*cols;
  }
  int
  row() const
  {
    return rows;
  }
  int
  col() const
  {
    return cols;
  }
};

class MMMatrixOut
{
public:
  static void write(const std::string &fname, const GeneralMatrix &m);
};

#endif /* MM_MATRIX_H */
