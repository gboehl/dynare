/*
 * Copyright © 2005 Ondra Kamenik
 * Copyright © 2019-2023 Dynare Team
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

// Exception

/* This is a simple code defining an exception and two convenience macros. */

#include <iostream>
#include <string>

#ifndef KORD_EXCEPTION_H
# define KORD_EXCEPTION_H

# define KORD_RAISE(mes) throw KordException(__FILE__, __LINE__, mes);

# define KORD_RAISE_IF(expr, mes)                                                                  \
  if (expr)                                                                                        \
   throw KordException(__FILE__, __LINE__, mes);

# define KORD_RAISE_X(mes, c) throw KordException(__FILE__, __LINE__, mes, c);

# define KORD_RAISE_IF_X(expr, mes, c)                                                             \
  if (expr)                                                                                        \
   throw KordException(__FILE__, __LINE__, mes, c);

class KordException
{
protected:
  std::string fname;
  int lnum;
  std::string message;
  int cd;

public:
  KordException(std::string f, int l, std::string mes, int c = 255) :
      fname {std::move(f)}, lnum {l}, message {std::move(mes)}, cd {c}
  {
  }
  virtual ~KordException() = default;
  virtual void
  print() const
  {
    std::cout << "At " << fname << ':' << lnum << ":(" << cd << "):" << message << '\n';
  }
  [[nodiscard]] virtual int
  code() const
  {
    return cd;
  }
  [[nodiscard]] const std::string&
  get_message() const
  {
    return message;
  }
};

// KordException error code definitions
constexpr int KORD_FP_NOT_CONV = 254;
constexpr int KORD_FP_NOT_FINITE = 253;
constexpr int KORD_MD_NOT_STABLE = 252;

#endif
