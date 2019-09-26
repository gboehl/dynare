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

#ifndef PASCAL_TRIANGLE_H
#define PASCAL_TRIANGLE_H

#include <vector>

class PascalRow : public std::vector<int>
{
  int k{1};
public:
  PascalRow() : std::vector<int>{}
  {
    push_back(2);
  }
  void setFromPrevious(const PascalRow &prev);
  void prolong(const PascalRow &prev);
  void prolongFirst(int n);
  void print() const;
};

namespace PascalTriangle
{
  void ensure(int n, int k);
  /*                              ⎛n⎞
    Computes binomial coefficient ⎝k⎠, hence the function name (“n over k”).
  */
  int noverk(int n, int k);
  void print();
};

#endif
