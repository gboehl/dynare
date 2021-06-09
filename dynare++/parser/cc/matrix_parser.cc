/*
 * Copyright © 2006 Ondra Kamenik
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

#include "parser_exception.hh"
#include "matrix_parser.hh"
#include "location.hh"
#include "matrix_tab.hh"

#include <algorithm>

using namespace ogp;

/** A global symbol for passing info to the MatrixParser from
 * matrix_parse(). */
MatrixParser *mparser;

/** The declaration of functions defined in matrix_ll.cc and
 * matrix_tab.cc generated from matrix.lex and matrix.y. */
void *matrix__scan_string(const char *);
void matrix__destroy_buffer(void *);
int matrix_parse();
extern ogp::location_type matrix_lloc;

void
MatrixParser::parse(const string &stream)
{
  // reinitialize the object
  data.clear();
  row_lengths.clear();
  nc = 0;
  // allocate temporary buffer and parse
  matrix_lloc.off = 0;
  matrix_lloc.ll = 0;
  void *p = matrix__scan_string(stream.c_str());
  mparser = this;
  matrix_parse();
  matrix__destroy_buffer(p);
}

void
MatrixParser::add_item(double v)
{
  data.push_back(v);
  if (row_lengths.size() == 0)
    row_lengths.push_back(0);
  (row_lengths.back())++;
  nc = std::max(nc, row_lengths.back());
}

void
MatrixParser::start_row()
{
  row_lengths.push_back(0);
}

void
MatrixParser::error(string mes) const
{
  throw ParserException(std::move(mes), matrix_lloc.off);
}

int
MatrixParser::find_first_non_empty_row(int start) const
{
  int r = start;
  while (r < static_cast<int>(row_lengths.size()) && row_lengths[r] == 0)
    r++;
  return r;
}

MPIterator
MatrixParser::begin() const
{
  MPIterator it(*this);
  return it;
}

MPIterator
MatrixParser::end() const
{
  MPIterator it(*this, "end");
  return it;
}

MPIterator::MPIterator(const MatrixParser &mp)
  : p(&mp), i(0), r(mp.find_first_non_empty_row())
{
}

MPIterator::MPIterator(const MatrixParser &mp, const char *dummy)
  : p(&mp), i(mp.data.size()), r(mp.row_lengths.size())
{
}

MPIterator &
MPIterator::operator++()
{
  i++;
  c++;
  if (p->row_lengths[r] <= c)
    {
      c = 0;
      r = p->find_first_non_empty_row(r+1);
    }
  return *this;
}
