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

#ifndef DYNARE_EXCEPTION_H
#define DYNARE_EXCEPTION_H

#include <string>
#include <utility>

class DynareException
{
  std::string mes;
public:
  DynareException(const std::string &m, const std::string &fname, int line, int col)
    : mes{"Parse error at " + fname + ", line " + std::to_string(line) + ", column "
          + std::to_string(col) + ": " + m}
  {
  }
  DynareException(const std::string &fname, int line, const std::string &m)
    : mes{fname + ':' + std::to_string(line) + ": " + m}
  {
  }
  DynareException(const std::string &m, int offset)
    : mes{"Parse error in provided string at offset " + std::to_string(offset) + ": " + m}
  {
  }
  const std::string &
  message() const
  {
    return mes;
  }
};

#endif
