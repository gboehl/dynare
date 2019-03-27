/* $Header: /var/lib/cvs/dynare_cpp/sylv/cc/KronUtils.h,v 1.1.1.1 2004/06/04 13:00:31 kamenik Exp $ */

/* Tag $Name:  $ */

#ifndef KRON_UTILS_H
#define KRON_UTILS_H

#include "KronVector.hh"
#include "QuasiTriangular.hh"

class KronUtils
{
public:
  /* Computes x = (Iₘ⊗…⊗Iₘ⊗T⊗Iₘ⊗…⊗Iₘ⊗Iₙ)·x, where x is n×mᵈ.
     T must be m×m, number of ⊗ is d, level is the number of Iₘ’s
     between T and Iₙ plus 1. If level=0, then we multiply by Iₘ⊗…⊗Iₘ⊗T,
     T must be n×n. */
  static void multAtLevel(int level, const QuasiTriangular &t,
                          KronVector &x);
  static void multAtLevelTrans(int level, const QuasiTriangular &t,
                               KronVector &x);

  // Computes x=(Fᵀ⊗Fᵀ⊗…⊗K)·x
  static void multKron(const QuasiTriangular &f, const QuasiTriangular &k,
                       KronVector &x);
};

#endif /* KRON_UTILS_H */
