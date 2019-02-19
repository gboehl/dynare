// Copyright 2004, Ondra Kamenik

#include "tl_static.hh"

TLStatic tls;

/* Note that we allow for repeated calls of |init|. This is not normal
   and the only purpose of allowing this is the test suite. */

TLStatic::TLStatic()
{
  ebundle = nullptr;
  pbundle = nullptr;
}

TLStatic::~TLStatic()
{
  if (ebundle)
    delete ebundle;
  if (pbundle)
    delete pbundle;
}

void
TLStatic::init(int dim, int nvar)
{
  if (ebundle)
    ebundle->generateUpTo(dim);
  else
    ebundle = new EquivalenceBundle(dim);

  if (pbundle)
    pbundle->generateUpTo(dim);
  else
    pbundle = new PermutationBundle(dim);
}
