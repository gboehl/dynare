/* $Header: /var/lib/cvs/dynare_cpp/sylv/cc/IterativeSylvester.h,v 1.1.1.1 2004/06/04 13:00:20 kamenik Exp $ */

/* Tag $Name:  $ */

#ifndef ITERATIVE_SYLVESTER_H
#define ITERATIVE_SYLVESTER_H

#include "SylvesterSolver.hh"
#include "KronVector.hh"
#include "QuasiTriangular.hh"
#include "SimilarityDecomp.hh"

class IterativeSylvester : public SylvesterSolver
{
public:
  IterativeSylvester(const QuasiTriangular &k, const QuasiTriangular &f)
    : SylvesterSolver(k, f)
  {
  }
  IterativeSylvester(const SchurDecompZero &kdecomp, const SchurDecomp &fdecomp)
    : SylvesterSolver(kdecomp, fdecomp)
  {
  }
  IterativeSylvester(const SchurDecompZero &kdecomp, const SimilarityDecomp &fdecomp)
    : SylvesterSolver(kdecomp, fdecomp)
  {
  }
  void solve(SylvParams &pars, KronVector &x) const override;
private:
  double performFirstStep(KronVector &x) const;
  static double performStep(const QuasiTriangular &k, const QuasiTriangular &f,
                            KronVector &x);
};

#endif /* ITERATIVE_SYLVESTER_H */
