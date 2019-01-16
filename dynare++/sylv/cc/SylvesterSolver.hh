/* $Header: /var/lib/cvs/dynare_cpp/sylv/cc/SylvesterSolver.h,v 1.1.1.1 2004/06/04 13:00:54 kamenik Exp $ */

/* Tag $Name:  $ */

#ifndef SYLVESTER_SOLVER_H
#define SYLVESTER_SOLVER_H

#include "KronVector.hh"
#include "QuasiTriangular.hh"
#include "QuasiTriangularZero.hh"
#include "SimilarityDecomp.hh"
#include "SylvParams.hh"
#include "SchurDecomp.hh"

#include <memory>

class SylvesterSolver
{
protected:
  const std::unique_ptr<const QuasiTriangular> matrixK;
  const std::unique_ptr<const QuasiTriangular> matrixF;
private:
  /* return true when it is more efficient to use QuasiTriangular
   * than QuasiTriangularZero */
  static bool
  zeroPad(const SchurDecompZero &kdecomp)
  {
    return ((kdecomp.getZeroCols()*3 < kdecomp.getDim()*2)
            || (kdecomp.getZeroCols() < 10));
  }
public:
  SylvesterSolver(const QuasiTriangular &k, const QuasiTriangular &f)
    : matrixK(std::make_unique<QuasiTriangular>(k)),
      matrixF(std::make_unique<QuasiTriangular>(f))
  {
  }
  SylvesterSolver(const SchurDecompZero &kdecomp, const SchurDecomp &fdecomp)
    : matrixK((zeroPad(kdecomp)) ?
              std::make_unique<QuasiTriangular>(kdecomp) :
              std::make_unique<QuasiTriangularZero>(kdecomp)),
      matrixF(std::make_unique<QuasiTriangular>(fdecomp))
  {
  }
  SylvesterSolver(const SchurDecompZero &kdecomp, const SimilarityDecomp &fdecomp)
    : matrixK((zeroPad(kdecomp)) ?
              std::make_unique<QuasiTriangular>(kdecomp) :
              std::make_unique<QuasiTriangularZero>(kdecomp)),
      matrixF(std::make_unique<BlockDiagonal>(fdecomp.getB()))
  {
  }
  virtual ~SylvesterSolver() = default;
  virtual void solve(SylvParams &pars, KronVector &x) const = 0;
};

#endif /* SYLVESTER_SOLVER_H */
