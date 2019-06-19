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
 * along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "tl_static.hh"
#include "pascal_triangle.hh"
#include "tl_exception.hh"

#include <mutex>
#include <limits>
#include <cmath>

/* Note that we allow for repeated calls of init(). This is not normal
   and the only purpose of allowing this is the test suite. */

namespace TLStatic
{
  EquivalenceBundle ebundle(1);
  PermutationBundle pbundle(1);
  std::mutex mut;

  const EquivalenceSet &
  getEquiv(int n)
  {
    return ebundle.get(n);
  }

  const PermutationSet &
  getPerm(int n)
  {
    return pbundle.get(n);
  }

  void
  init(int dim, int nvar)
  {
    // Check that tensor indices will not overflow (they are stored as signed int, hence on 31 bits)
    if (std::log2(nvar)*dim > std::numeric_limits<int>::digits)
      throw TLException(__FILE__, __LINE__, "Problem too large, you should decrease the approximation order");

    std::lock_guard<std::mutex>{mut};
    ebundle.generateUpTo(dim);
    pbundle.generateUpTo(dim);

    PascalTriangle::ensure(nvar, dim);
  }
}
