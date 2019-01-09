/* $Header: /var/lib/cvs/dynare_cpp/sylv/cc/SchurDecompEig.h,v 1.1.1.1 2004/06/04 13:00:44 kamenik Exp $ */

/* Tag $Name:  $ */

// contains algorithms for eigenvalue reordering

#ifndef SCHUR_DECOMP_EIG_H
#define SCHUR_DECOMP_EIG_H

#include "SchurDecomp.hh"
#include "QuasiTriangular.hh"

class SchurDecompEig : public SchurDecomp
{
public:
  using diag_iter = QuasiTriangular::diag_iter;
  SchurDecompEig(const SqSylvMatrix &m) : SchurDecomp(m)
  {
  }
  SchurDecompEig(const QuasiTriangular &tr) : SchurDecomp(tr)
  {
  };
  SchurDecompEig(QuasiTriangular &tr) : SchurDecomp(tr)
  {
  }
  diag_iter bubbleEigen(diag_iter from, diag_iter to);
  void orderEigen();
protected:
  bool tryToSwap(diag_iter &it, diag_iter &itadd);
};

#endif /* SCHUR_DECOMP_EIG_H */

// Local Variables:
// mode:C++
// End:
