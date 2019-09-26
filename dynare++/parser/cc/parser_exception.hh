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

#ifndef OG_FORMULA_PARSER_H
#define OG_FORMULA_PARSER_H

#include <string>

namespace ogp
{
  using std::string;

  /** This is an easy exception, which, besides the message, stores
   * also an offset of the parse error. Since we might need to track
   * the argument number and for example the filed in the argument
   * which caused the error, we add three integers, which have no
   * semantics here. They should be documented in the function which
   * throws an exception and sets them. Their default value is -1,
   * which means they have not been set. */
  class ParserException
  {
  protected:
    string mes;
    int off;
    int aux_i1;
    int aux_i2;
    int aux_i3;
  public:
    ParserException(string m, int offset);
    ParserException(string m, const char *dum, int i1);
    ParserException(string m, const char *dum, int i1, int i2);
    ParserException(string m, const char *dum, int i1, int i2, int i3);
    ParserException(const ParserException &e, int plus_offset);
    /** Makes a copy and pushes given integer to aux_i1 shuffling
     * others and forgetting the last. */
    ParserException(const ParserException &e, const char *dum, int i);
    /** Makes a copy and pushes given two integers to aux_i1 and aux_i2  shuffling
     * others and forgetting the last two. */
    ParserException(const ParserException &e, const char *dum, int i1, int i2);
    /** Makes a copy and pushes given three integers to aux_i1, aux_i2, aus_i3 shuffling
     * others and forgetting the last three. */
    ParserException(const ParserException &e, const char *dum, int i1, int i2, int i3);
    ParserException(const ParserException &e) = default;
    virtual ~ParserException() = default;
    const string &
    message() const
    {
      return mes;
    }
    int
    offset() const
    {
      return off;
    }
    const int &
    i1() const
    {
      return aux_i1;
    }
    int &
    i1()
    {
      return aux_i1;
    }
    const int &
    i2() const
    {
      return aux_i2;
    }
    int &
    i2()
    {
      return aux_i2;
    }
    const int &
    i3() const
    {
      return aux_i3;
    }
    int &
    i3()
    {
      return aux_i3;
    }
  protected:
    void copy(const ParserException &e);
  };
};

#endif
