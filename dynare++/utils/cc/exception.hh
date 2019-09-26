/*
 * Copyright © 2005 Ondra Kamenik
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

#ifndef OGU_EXCEPTION_H
#define OGU_EXCEPTION_H

#include <string>
#include <iostream>
#include <utility>

namespace ogu
{
  /** A primitive exception. */
  class Exception
  {
  protected:
    const std::string file;
    const int line;
    const std::string mes;
  public:
    Exception(std::string file_arg, int line_arg, std::string mes_arg)
      : file{std::move(file_arg)},
        line{line_arg},
        mes{std::move(mes_arg)}
    {
    }
    virtual ~Exception() = default;

    void
    print(std::ostream &out) const
    {
      out << file << ':' << line << ": " << mes << std::endl;
    }

    void
    print() const
    {
      print(std::cout);
    }

    std::string
    message() const
    {
      return mes;
    }
  };
};

#endif
