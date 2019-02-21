// Copyright 2004, Ondra Kamenik

#include "tl_static.hh"
#include "pascal_triangle.hh"

#include <mutex>

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
    std::lock_guard<std::mutex>{mut};
    ebundle.generateUpTo(dim);
    pbundle.generateUpTo(dim);

    PascalTriangle::ensure(nvar, dim);
  }
}
