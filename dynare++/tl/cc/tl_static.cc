// Copyright 2004, Ondra Kamenik

#include "tl_static.hh"
#include "tl_exception.hh"

TLStatic tls;

/* Note that we allow for repeated calls of |init|. This is not normal
   and the only purpose of allowing this is the test suite. */

TLStatic::TLStatic()
{
  ebundle = NULL;
  pbundle = NULL;
  ptriang = NULL;
}

TLStatic::~TLStatic()
{
  if (ebundle)
    delete ebundle;
  if (pbundle)
    delete pbundle;
  if (ptriang)
    delete ptriang;
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

  if (ptriang)
    delete ptriang;
  ptriang = new PascalTriangle(nvar, dim);
}

/* The coefficients are stored in |data| row by row where a row are
   coeffs with the same $k$.

   We first initialize the first row with ones. Then for each other row
   we initialize the first item to one, and other items are a sum of
   coefficients of $n-1$ which is in code |i+j-1|. */

PascalTriangle::PascalTriangle(int n, int k)
  : data(new int[(n+1)*(k+1)]), kmax(k), nmax(n)
{
  for (int i = 0; i <= n; i++)
    data[i] = 1;
  for (int j = 1; j <= k; j++)
    {
      data[j*(nmax+1)] = 1;
      for (int i = 1; i <= n; i++)
        data[j*(nmax+1)+i] = noverk(i+j-1, j) + noverk(i+j-1, j-1);
    }
}

/* Clear. Recall, that there are |nmax+1| items in a row. */

int
PascalTriangle::noverk(int n, int k) const
{
  TL_RAISE_IF(k > n || n < 0,
              "Wrong arguments for PascalTriangle::noverk");

  if (k <= kmax && n-k <= nmax)
    return data[k*(nmax+1)+n-k];

  if (n-k <= kmax && k <= nmax)
    return data[(n-k)*(nmax+1)+k];

  TL_RAISE("n or k out of range in PascalTriangle::noverk");
  return 0;
}
