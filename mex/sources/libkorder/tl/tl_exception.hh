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

// Exception.

/* Within the code we often check some state of variables, typically
   preconditions or postconditions. If the state is not as required, it
   is worthless to continue, since this means some fatal error in
   algorithms. In this case we raise an exception which can be caught at
   some higher level. This header file defines a simple infrastructure
   for this. */

#ifndef TL_EXCEPTION_H
#define TL_EXCEPTION_H

#include <iostream>
#include <utility>
#include <string>

/* The basic idea of raising an exception if some condition fails is that the
   conditions is checked only if required. We define global TL_DEBUG macro
   which is integer and says, how many debug messages the programm has to emit.
   We also define TL_DEBUG_EXCEPTION which says, for what values of TL_DEBUG we
   will check for conditions of the exceptions. If the TL_DEBUG is equal or
   higher than TL_DEBUG_EXCEPTION, the exception conditions are checked.

   We define TL_RAISE, and TL_RAISE_IF macros which throw an instance of
   TLException (only if TL_DEBUG >= TL_DEBUG_EXCEPTION for the latter). The
   first is unconditional throw, the second is conditioned by a given
   expression. Note that if TL_DEBUG < TL_DEBUG_EXCEPTION then the code is
   compiled but evaluation of the condition is passed. If code is optimized,
   the optimizer also passes evaluation of TL_DEBUG and TL_DEBUG_EXCEPTION
   comparison (I hope).

   We provide default values for TL_DEBUG and TL_DEBUG_EXCEPTION. */

#ifndef TL_DEBUG_EXCEPTION
# define TL_DEBUG_EXCEPTION 1
#endif

#ifndef TL_DEBUG
# define TL_DEBUG 0
#endif

#define TL_RAISE(mes)                           \
  throw TLException(__FILE__, __LINE__, mes)

#define TL_RAISE_IF(expr, mes)                                          \
  if (TL_DEBUG >= TL_DEBUG_EXCEPTION && (expr)) throw TLException(__FILE__, __LINE__, mes);

/* Primitive exception class containing file name, line number and message. */

class TLException
{
  const std::string fname;
  int lnum;
public:
  const std::string message;
  TLException(std::string fname_arg, int lnum_arg, std::string message_arg)
    : fname{std::move(fname_arg)},
      lnum{lnum_arg},
      message{std::move(message_arg)}
  {
  }
  virtual ~TLException() = default;
  virtual void
  print() const
  {
    std::cout << "At " << fname << ':' << lnum << ':' << message << std::endl;
  }
};

#endif
