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
 * along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "parser_exception.hh"

using namespace ogp;

ParserException::ParserException(string m, int offset)
  : mes(std::move(m)), off(offset),
    aux_i1(-1), aux_i2(-1), aux_i3(-1)
{
}

ParserException::ParserException(string m, const char *dum, int i1)
  : mes(std::move(m)), off(0),
    aux_i1(i1), aux_i2(-1), aux_i3(-1)
{
}

ParserException::ParserException(string m, const char *dum, int i1, int i2)
  : mes(std::move(m)), off(0),
    aux_i1(i1), aux_i2(i2), aux_i3(-1)
{
}

ParserException::ParserException(string m, const char *dum, int i1, int i2, int i3)
  : mes(std::move(m)), off(0),
    aux_i1(i1), aux_i2(i2), aux_i3(i3)
{
}

ParserException::ParserException(const ParserException &m, int plus_offset)
  : aux_i1(-1), aux_i2(-1), aux_i3(-1)
{
  copy(m);
  off += plus_offset;
}

ParserException::ParserException(const ParserException &m, const char *dum, int i)
  : aux_i1(-1), aux_i2(-1), aux_i3(-1)
{
  copy(m);
  aux_i3 = m.aux_i2;
  aux_i2 = m.aux_i1;
  aux_i1 = i;
}

ParserException::ParserException(const ParserException &m, const char *dum, int i1, int i2)
  : aux_i1(-1), aux_i2(-1), aux_i3(-1)
{
  copy(m);
  aux_i3 = m.aux_i1;
  aux_i2 = i2;
  aux_i1 = i1;
}

ParserException::ParserException(const ParserException &m, const char *dum, int i1, int i2, int i3)
  : aux_i1(-1), aux_i2(-1), aux_i3(-1)
{
  copy(m);
  aux_i3 = i3;
  aux_i2 = i2;
  aux_i1 = i1;
}

void
ParserException::copy(const ParserException &e)
{
  mes = e.mes;
  off = e.off;
  aux_i1 = e.aux_i1;
  aux_i2 = e.aux_i2;
  aux_i3 = e.aux_i3;
}
