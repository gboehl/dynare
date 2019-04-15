// Copyright 2004, Ondra Kamenik

#include "tl_static.hh"
#include "pascal_triangle.hh"
#include "tl_exception.hh"

#include <mutex>
#include <limits>
#include <cmath>

/* Note that we allow for repeated calls of |init|. This is not normal
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
