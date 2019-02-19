// Copyright 2004, Ondra Kamenik

// Tensor library static data.

/* The purpose of this file is to make a unique static variable which
   would contain all other static variables and be responsible for their
   correct initialization and destruction. The variables include an
   equivalence bundle and a permutation bundle. Both depend on dimension of the
   problem, and maximum number of variables.

   So we declare static |tls| variable of type |TLStatic| encapsulating
   the variables. The |tls| must be initialized at the beginning of
   the program, as dimension and number of variables is known. */

#ifndef TL_STATIC_H
#define TL_STATIC_H

#include "equivalence.hh"
#include "permutation.hh"

struct TLStatic
{
  EquivalenceBundle *ebundle;
  PermutationBundle *pbundle;

  TLStatic();
  ~TLStatic();
  void init(int dim, int nvar);
};

extern TLStatic tls;

#endif
