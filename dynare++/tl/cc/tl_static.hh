// Copyright 2004, Ondra Kamenik

// Tensor library static data.

/* The purpose of this file is to make a unique static variable which
   would contain all other static variables and be responsible for their
   correct initialization and destruction. The variables include an
   equivalence bundle and a Pascal triangle for binomial
   coefficients. Both depend on dimension of the problem, and maximum
   number of variables.

   So we declare static |tls| variable of type |TLStatic| encapsulating
   the variables. The |tls| must be initialized at the beginning of
   the program, as dimension and number of variables is known.

   Also we define a class for Pascal triangle. */

#ifndef TL_STATIC_H
#define TL_STATIC_H

#include "equivalence.hh"
#include "permutation.hh"

/* Pascal triangle is a storage for binomial coefficients. We store in
   |data| array the coefficients of rectangle starting at $\pmatrix{0\cr
   0}$, and ending $\pmatrix{nmax+kmax\cr kmax}$. */

class PascalTriangle
{
  int *data;
  int kmax;
  int nmax;
public:
  PascalTriangle(int n, int k);
  ~PascalTriangle()
  {
    delete [] data;
  }
  int noverk(int n, int k) const;
};

struct TLStatic
{
  EquivalenceBundle *ebundle;
  PermutationBundle *pbundle;
  PascalTriangle *ptriang;

  TLStatic();
  ~TLStatic();
  void init(int dim, int nvar);
};

extern TLStatic tls;

#endif
